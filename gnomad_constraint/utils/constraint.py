"""Script containing utility functions used in the constraint pipeline."""
import logging
from typing import Dict, Optional, Tuple

import hail as hl
import numpy as np
from gnomad.utils.constraint import (
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    collapse_strand,
    compute_expected_variants,
    compute_pli,
    count_variants_by_group,
    get_constraint_flags,
    get_downsampling_freq_indices,
    oe_aggregation_expr,
    oe_confidence_interval,
    trimer_from_heptamer,
)
from gnomad.utils.filtering import (
    add_filters_expr,
    filter_by_numeric_expr_range,
    filter_for_mu,
    filter_to_autosomes,
)
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import (
    COVERAGE_CUTOFF,
    get_checkpoint_path,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def add_vep_context_annotations(
    ht: hl.Table, annotated_context_ht: hl.Table
) -> hl.Table:
    """
    Add annotations from VEP context Table to gnomAD data.

    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp

    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from
    gnomAD data.

    .. note::
        Function expects that multiallelic variants in the VEP context Table have been
        split.

    :param ht: gnomAD exomes or genomes public Hail Table.
    :param annotated_context_ht: VEP context Table.
    :return: Table with annotations.
    """
    context_ht = annotated_context_ht.drop("a_index", "was_split")
    context_ht = context_ht.annotate(vep=context_ht.vep.drop("colocated_variants"))
    return ht.annotate(**context_ht[ht.key])


def prepare_ht_for_constraint_calculations(
    ht: hl.Table, require_exome_coverage: bool = True
) -> hl.Table:
    """
    Filter input Table and add annotations used in constraint calculations.

    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - exome_coverage
        - pass_filters - Whether the variant passed all variant filters
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and
          `add_most_severe_csq_to_tc_within_vep_root()`

    :param ht: Input Table to be annotated.
    :param require_exome_coverage: Filter to sites where exome coverage is defined.
        Default is True.
    :return: Table with annotations.
    """
    ht = trimer_from_heptamer(ht)

    if "filters" in ht.row_value.keys():
        ht = ht.annotate(pass_filters=hl.len(ht.filters) == 0)

    # Add annotations for 'ref' and 'alt'.
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])

    # Filter to SNPs and context fields where the bases are either A, T, C, or G.
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt) & ht.context.matches(f"[ATCG]{{{3}}}"))

    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and
    # collapse strands to deduplicate the context.
    ht = annotate_mutation_type(collapse_strand(ht))

    # Obtain field name for median exome coverage.
    # TODO: Edit coverage field once decide what to use for v4.
    exome_median_cov_field = (
        "median_approx" if "median_approx" in ht.coverage.exomes else "median"
    )

    # Define methylation level cutoffs based on fields present in the 'methylation'
    # annotation.
    if "MEAN" in ht.methylation:
        # The GRCh37 methylation resource provides a MEAN score ranging from 0-1.
        methylation_expr = ht.methylation.MEAN
        methylation_cutoffs = (0.6, 0.2)
    elif "methylation_level" in ht.methylation:
        # The GRCh38 methylation resource provides a score ranging from 0-15 for autosomes. The
        # determination of this score is described in Chen et al:
        # https://www.biorxiv.org/content/10.1101/2022.03.20.485034v2.full
        # For chrX, methylation scores reange from 0-12, but these scores are not directly comparable
        # to the autosome scores (chrX and autosomes were analyzed separately and levels are relative).
        # Cutoffs to translate these scores to the 0-2 methylation level were determined by
        # correlating these scores with the GRCh37 liftover scores. Proposed cutoffs are:
        # 0, 1-5, 6+ for autosomes, and 0, 1-3, 4+ for chrX.
        methylation_expr = ht.methylation.methylation_level
        methylation_cutoffs = hl.if_else(ht.locus.contig != "chrX", (5, 0), (3, 0))
    else:
        raise ValueError(
            "No 'methylation_level' or 'MEAN' found in 'methylation' annotation."
        )

    # Add annotations for methylation level and median exome coverage.
    ht = ht.annotate(
        methylation_level=(
            hl.case()
            .when(ht.cpg & (methylation_expr > methylation_cutoffs[0]), 2)
            .when(ht.cpg & (methylation_expr > methylation_cutoffs[1]), 1)
            .default(0)
        ),
        exome_coverage=ht.coverage.exomes[exome_median_cov_field],
    )

    # Add most_severe_consequence annotation to 'transcript_consequences' within the
    # vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)

    if require_exome_coverage:
        # Filter out locus with undefined exome coverage.
        ht = ht.filter(hl.is_defined(ht.exome_coverage))

    return ht


def create_observed_and_possible_ht(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    pops: Tuple[str] = (),
    grouping: Tuple[str] = ("exome_coverage",),
    partition_hint: int = 100,
    filter_coverage_over_0: bool = False,
    low_coverage_filter: int = None,
    transcript_for_synonymous_filter: str = None,
    global_annotation: Optional[str] = None,
) -> hl.Table:
    """
    Count the observed variants and possible variants by substitution, context, methylation level, and additional `grouping`.

    Prior to computing variant counts the following variants are removed:
        - Variants not observed by any samples in the dataset: `(freq_expr.AC > 0)`
        - Low-quality variants: `exome_ht.pass_filters`
        - Variants with allele frequency above `max_af` cutoff: `(freq_expr.AF <=
          max_af)`
        - Variants that are not synonymous or in the canonical/MANE Select transcript if specified

    For each substitution, context, methylation level, and exome coverage, the rest of
    variants in `exome_ht` are counted and annotated as `observed_variants`, and the
    rest of variants in `context_ht` are counted and annotated as `possible_variants`.
    The final Table is the outer-join of the filtered `exome_ht` and `context_ht` with
    the `observed_variants` and `possible_variants` annotations.

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation_level
        - observed_variants - observed variant counts in `exome_ht`
        - possible_variants - possible variant counts in `context_ht`
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param grouping: Annotations other than 'context', 'ref', 'alt', and
        `methylation_level` to group by when counting variants. Default is
        ('exome_coverage',).
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :param filter_coverage_over_0: Whether to filter the exome Table and context Table
        to variants with coverage larger than 0. Default is False.
    :param low_coverage_filter: Lower median coverage cutoff for coverage filter. Sites
        with coverage below this cutoff will be removed from the `exome_ht` and
        'context_ht'.
    :param transcript_for_synonymous_filter: Transcript to use when filtering to
        synonymous variants. Choices: ["mane_select", "canonical", None]. If "canonical", will
        filter to variants with a synonymous consequence in Ensembl canonical
        transcripts. If "mane_select", will filter to variants with a synonymous consequence
        in MANE Select transcripts. If None, no transcript/synonymous filter will be
        applied. Default is None.
    :param global_annotation: The annotation name to use as a global StructExpression
        annotation containing input parameter values. If no value is supplied, this
        global annotation will not be added. Default is None.
    :return: Table with observed variant and possible variant count.
    """
    if low_coverage_filter is not None:
        context_ht = context_ht.filter(context_ht.exome_coverage >= low_coverage_filter)
        exome_ht = exome_ht.filter(exome_ht.exome_coverage >= low_coverage_filter)

    # Allele frequency information for high-quality genotypes (GQ >= 20; DP >= 10; and
    # AB >= 0.2 for heterozygous calls) in all release samples in gnomAD.
    freq_expr = exome_ht.freq[0]

    # Set up the criteria to exclude variants not observed in the dataset, low-quality
    # variants, variants with allele frequency above the `max_af` cutoff, and variants
    # with exome coverage larger than 0 if requested.
    keep_criteria = (
        (freq_expr.AC > 0) & exome_ht.pass_filters & (freq_expr.AF <= max_af)
    )
    if filter_coverage_over_0:
        keep_criteria &= exome_ht.coverage > 0

    keep_annotations += grouping

    # Keep variants that satisfy the criteria above.
    filtered_exome_ht = exome_ht.filter(keep_criteria)

    # Filter context ht to sites with defined exome coverage.
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    # If requested keep only variants that are synonymous in either MANE Select or
    # canonical transcripts.
    if transcript_for_synonymous_filter is not None:
        if transcript_for_synonymous_filter == "canonical":
            canonical, mane_select = True, False
        elif transcript_for_synonymous_filter == "mane_select":
            canonical, mane_select = False, True
        else:
            raise ValueError(
                "If transcript_for_synonymous_filter is not None, must be either"
                " 'canonical' or 'mane_select'"
            )
        filtered_exome_ht = filter_vep_transcript_csqs(
            exome_ht.filter(keep_criteria), canonical=canonical, mane_select=mane_select
        )
        context_ht = filter_vep_transcript_csqs(
            context_ht, canonical=canonical, mane_select=mane_select
        )
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `keep_annotations`.
    observed_ht = count_variants_by_group(
        filtered_exome_ht.select(*list(keep_annotations) + ["freq"]),
        additional_grouping=grouping,
        partition_hint=partition_hint,
        count_downsamplings=pops,
        use_table_group_by=True,
    )

    # TODO: Remove repartition once partition_hint bugs are resolved.
    observed_ht = observed_ht.repartition(partition_hint)
    observed_ht = observed_ht.transmute(observed_variants=observed_ht.variant_count)

    # Filter the `exome_ht` to rows that don’t match the criteria above.
    # Anti join the `context_ht` with filtered `exome_ht`, so that `context_ht` only
    # has rows that match the criteria above in the `exome_ht` or are never in
    # the `exome_ht`.
    context_ht = context_ht.select(*keep_annotations).anti_join(
        exome_ht.filter(keep_criteria, keep=False)
    )

    # Count the possible variants in the context Table grouped by `keep_annotations`.
    possible_ht = count_variants_by_group(
        context_ht,
        additional_grouping=grouping,
        partition_hint=partition_hint,
        use_table_group_by=True,
    )
    possible_ht = annotate_with_mu(possible_ht, mutation_ht)
    possible_ht = possible_ht.transmute(possible_variants=possible_ht.variant_count)

    # Outer join the Tables with possible variant counts and observed variant counts.
    ht = observed_ht.join(possible_ht, "outer")
    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    # Annotate the Table with 'cpg' and 'mutation_type' (one of "CpG", "non-CpG
    # transition", or "transversion").
    ht = annotate_mutation_type(ht)

    if global_annotation:
        ht = ht.annotate_globals(
            **{global_annotation: hl.struct(max_af=max_af, pops=pops)}
        )

    return ht


def apply_models(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Optional[Tuple[float, float]] = None,
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    pops: Tuple[str] = (),
    obs_pos_count_partition_hint: int = 2000,
    expected_variant_partition_hint: int = 1000,
    custom_vep_annotation: str = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    low_coverage_filter: int = None,
    use_mane_select_instead_of_canonical: bool = False,
) -> hl.Table:
    """
    Compute the expected number of variants and observed:expected ratio using plateau models and coverage model.

    This function sums the number of possible variants times the mutation rate for all
    variants, and applies the calibration model separately for CpG transitions and
    other sites. For sites with coverage lower than the coverage cutoff, the value
    obtained from the previous step is multiplied by the coverage correction factor.
    These values are summed across the set of variants of interest to obtain the
    expected number of variants.

    A brief view of how to get the expected number of variants:
        mu_agg = the number of possible variants * the mutation rate (all variants)
        predicted_proportion_observed = sum(plateau model slop * mu_agg + plateau model intercept) (separately for CpG transitions and other sites)
        if 0 < coverage < coverage cutoff:
            coverage_correction = coverage_model slope * log10(coverage) + coverage_model intercept
            expected_variants = sum(predicted_proportion_observed * coverage_correction)
        else:
            expected_variants = sum(predicted_proportion_observed)
        The expected_variants are summed across the set of variants of interest to
        obtain the final expected number of variants.

    Function adds the following annotations all grouped by groupings (output of
        `annotate_exploded_vep_for_constraint_groupings()`):
        - observed_variants - observed variant counts annotated by `count_variants`
          function
        - predicted_proportion_observed (including those for each population) - the sum
          of mutation rate adjusted by plateau models and possible variant counts
        - possible_variants (including those for each population if `pops` is
          specified) - the sum of possible variant counts derived from the context
          Table
        - expected_variants (including those for each population if `pops` is
          specified) - the sum of expected variant counts
        - mu - sum(mu_snp * possible_variant * coverage_correction)
        - obs_exp - observed:expected ratio
        - annotations annotated by `annotate_exploded_vep_for_constraint_groupings()`

    :param exome_ht: Exome sites Table (output of `prepare_ht_for_constraint_calculations
        ()`) filtered to autosomes and pseudoautosomal regions.
    :param context_ht: Context Table (output of `prepare_ht_for_constraint_calculations
        ()`) filtered to autosomes and pseudoautosomal regions.
    :param mutation_ht: Mutation rate Table with 'mu_snp' field.
    :param plateau_models: Linear models (output of `build_models()` in
        gnomad_methods`), with the values of the dictionary formatted as a
        StrucExpression of intercept and slope, that calibrates mutation rate to
        proportion observed for high coverage exome. It includes models for CpG sites,
        non-CpG sites, and each population in `POPS`.
    :param coverage_model: A linear model (output of `build_models()` in
        gnomad_methods), formatted as a Tuple of intercept and slope, that calibrates a
        given coverage level to observed:expected ratio. It's a correction factor for
        low coverage sites.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table and exome Table.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param obs_pos_count_partition_hint: Target number of partitions for
        aggregation when counting variants. Default is 2000.
    :param expected_variant_partition_hint: Target number of partitions for sum
        aggregators when computation is done. Default is 1000.
    :param custom_vep_annotation: The customized model (one of
        "transcript_consequences" or "worst_csq_by_gene"), Default is None.
    :param high_cov_definition: Median coverage cutoff. Sites with coverage above this cutoff
        are considered well covered and was used to build plateau models. Sites
        below this cutoff have low coverage and was used to build coverage models.
        Default is `COVERAGE_CUTOFF`.
    :param low_coverage_filter: Lower median coverage cutoff for coverage filter.
        Sites with coverage below this cutoff will be removed from`exome_ht` and
        'context_ht'.
    :param use_mane_select_instead_of_canonical: Use MANE Select rather than canonical
        grouping. Only used when `custom_vep_annotation` is set to
        'transcript_consequences'.


    :return: Table with `expected_variants` (expected variant counts) and `obs_exp`
        (observed:expected ratio) annotations.
    """
    # Filter context ht to sites with defined exome coverage.
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    if low_coverage_filter is not None:
        context_ht = context_ht.filter(context_ht.exome_coverage >= low_coverage_filter)
        exome_ht = exome_ht.filter(exome_ht.exome_coverage >= low_coverage_filter)

    # Add necessary constraint annotations for grouping.
    if custom_vep_annotation == "worst_csq_by_gene":
        vep_annotation = "worst_csq_by_gene"
    else:
        vep_annotation = "transcript_consequences"
        if use_mane_select_instead_of_canonical:
            include_canonical_group, include_mane_select_group = False, True
        else:
            include_canonical_group, include_mane_select_group = True, False

    context_ht, _ = annotate_exploded_vep_for_constraint_groupings(
        ht=context_ht,
        vep_annotation=vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )
    exome_ht, grouping = annotate_exploded_vep_for_constraint_groupings(
        ht=exome_ht,
        vep_annotation=vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )

    # Compute observed and possible variant counts.
    ht = create_observed_and_possible_ht(
        exome_ht,
        context_ht,
        mutation_ht,
        max_af,
        keep_annotations,
        pops,
        grouping,
        obs_pos_count_partition_hint,
        filter_coverage_over_0=True,
        transcript_for_synonymous_filter=None,
    )

    mu_expr = ht.mu_snp * ht.possible_variants
    # Determine coverage correction to use based on coverage value. If no
    # coverage model is provided, set to 1 as long as coverage > 0.
    cov_corr_expr = (
        hl.case()
        .when(ht.coverage == 0, 0)
        .when(ht.coverage >= high_cov_definition, 1)
        .default(
            (coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0])
            if coverage_model is not None
            else 1
        )
    )
    # Generate sum aggregators for 'mu' on the entire dataset.
    agg_expr = {"mu": hl.agg.sum(mu_expr * cov_corr_expr)}
    agg_expr.update(
        compute_expected_variants(ht, plateau_models, mu_expr, cov_corr_expr, ht.cpg)
    )
    for pop in pops:
        agg_expr.update(
            compute_expected_variants(
                ht, plateau_models, mu_expr, cov_corr_expr, ht.cpg, pop
            )
        )

    grouping = list(grouping)
    grouping.remove("coverage")

    # Aggregate the sum aggregators grouped by `grouping`.
    ht = (
        ht.group_by(*grouping)
        .partition_hint(expected_variant_partition_hint)
        .aggregate(**agg_expr)
    )

    # TODO: Remove repartition once partition_hint bugs are resolved.
    ht = ht.repartition(expected_variant_partition_hint)

    # Annotate global annotations.
    coverage_model_global = coverage_model if coverage_model else "None"
    ht = ht.annotate_globals(
        apply_model_params=hl.struct(
            max_af=max_af,
            pops=pops,
            plateau_models=plateau_models,
            coverage_model=coverage_model_global,
        )
    )
    # Compute the observed:expected ratio.
    return ht.annotate(obs_exp=ht.observed_variants / ht.expected_variants)


def calculate_mu_by_downsampling(
    genome_ht: hl.Table,
    context_ht: hl.Table,
    recalculate_all_possible_summary: bool = True,
    omit_methylation: bool = False,
    count_singletons: bool = False,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    ac_cutoff: int = 5,
    downsampling_level: int = 1000,
    total_mu: float = 1.2e-08,
    pops: Tuple[str] = (),
    min_cov: int = 15,
    max_cov: int = 60,
    gerp_lower_cutoff: float = -3.9885,
    gerp_upper_cutoff: float = 2.6607,
) -> hl.Table:
    """
    Calculate mutation rate using the downsampling with size specified by `downsampling_level` in genome sites Table.

    Prior to computing mutation rate, only the following variants are kept:
        - variants with the mean coverage in the gnomAD genomes between `min_cov` and
          `max_cov`.
        - variants where the most severe consequence was 'intron_variant' or
          'intergenic_variant'.
        - variants with the GERP score between `gerp_lower_cutoff` and
          `gerp_upper_cutoff` (these default to -3.9885 and 2.6607, respectively -
          these values were precalculated on the GRCh37 context Table and define the
          5th and 95th percentiles).
        - high-quality variants: `genome_ht.pass_filters`.
        - variants with allele count below `ac_cutoff`: `(freq_expr.AC <= ac_cutoff)`.

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context.
        - ref - the reference allele.
        - alt - the alternate base.
        - methylation_level - methylation_level.
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`.
        - mu_snp - SNP mutation rate.
        - annotations added by `annotate_mutation_type`.

    :param genome_ht: Genome sites Table for autosome/pseudoautosomal regions.
    :param context_ht: Context Table for autosome/pseudoautosomal regions.
    :param recalculate_all_possible_summary: Whether to calculate possible
        variants using context Table with locus that is only on an autosome or
        in a pseudoautosomal region. Default is True.
    :param omit_methylation: Whether to omit 'methylation_level' from the
        grouping when counting variants. Default is False.
    :param count_singletons: Whether to count singletons. Default is False.
    :param keep_annotations: Annotations to keep in the context Table and genome
        sites Table.
    :param ac_cutoff: The cutoff of allele count when filtering context Table
        and genome sites Table.
    :param downsampling_level: The size of downsamplings will be used to count
        variants. Default is 1000.
    :param total_mu: The per-generation mutation rate. Default is 1.2e-08.
    :param pops: List of populations to use for downsampling counts. If empty
        Tuple is supplied, will default to '['global']'.
    :param min_cov: Minimum coverage required to keep a site when calculating
        the mutation rate. Default is 15.
    :param max_cov: Maximum coverage required to keep a site when calculating
        the mutation rate. Default is 60.
    :param gerp_lower_cutoff: Minimum GERP score for variant to be included
        when calculating the mutation rate. Default is -3.9885.
    :param gerp_upper_cutoff: Maximum GERP score for variant to be included
        when calculating the mutation rate. Default is 2.6607.
    :return: Mutation rate Table.
    """
    if not pops:
        pops = ["global"]

    # Filter to autosomal sites (remove pseudoautosomal regions) between
    # min_cov and max_cov.
    context_ht = filter_to_autosomes(
        filter_by_numeric_expr_range(
            context_ht, context_ht.coverage.genomes.mean, (min_cov, max_cov)
        )
    )
    genome_ht = filter_to_autosomes(
        filter_by_numeric_expr_range(
            genome_ht, genome_ht.coverage.genomes.mean, (min_cov, max_cov)
        )
    )

    # Filter the Table so that the most severe annotation is 'intron_variant' or
    # 'intergenic_variant', and that the GERP score is between 'gerp_lower_cutoff' and
    # 'gerp_upper_cutoff' (ideally these values will define the 5th and 95th
    # percentile of the genome-wide distribution).
    context_ht = filter_for_mu(context_ht, gerp_lower_cutoff, gerp_upper_cutoff)
    genome_ht = filter_for_mu(genome_ht, gerp_lower_cutoff, gerp_upper_cutoff)

    context_ht = context_ht.select(*keep_annotations)
    genome_ht = genome_ht.select(*list(keep_annotations) + ["freq", "pass_filters"])

    # Get the frequency index of downsampling with size of `downsampling_level`.
    downsampling_meta = get_downsampling_freq_indices(genome_ht.freq_meta)
    downsampling_idx = hl.eval(
        downsampling_meta.filter(
            lambda x: x[1]["downsampling"] == str(downsampling_level)
        )[0][0]
    )
    freq_expr = genome_ht.freq[downsampling_idx]

    # Set up the criteria to filter out low-quality sites, and sites found in greater
    # than 'ac_cutoff' copies in the downsampled set.
    keep_criteria = (freq_expr.AC <= ac_cutoff) & genome_ht.pass_filters

    # Count the observed variants in the genome sites Table.
    observed_ht = count_variants_by_group(
        genome_ht.filter(keep_criteria).select(*list(keep_annotations) + ["freq"]),
        count_downsamplings=pops,
        count_singletons=count_singletons,
        omit_methylation=omit_methylation,
        use_table_group_by=True,
    )

    # Count possible variants in context Table, only keeping variants not in the genome
    # dataset, or with AC <= 'ac_cutoff' and passing filters.
    all_possible_ht = count_variants_by_group(
        context_ht.anti_join(genome_ht.filter(keep_criteria, keep=False)).select(
            *keep_annotations
        ),
        omit_methylation=omit_methylation,
        use_table_group_by=True,
    )
    all_possible_ht = all_possible_ht.checkpoint(
        get_checkpoint_path("all_possible_summary"),
        _read_if_exists=not recalculate_all_possible_summary,
        overwrite=recalculate_all_possible_summary,
    )

    ht = observed_ht.annotate(
        possible_variants=all_possible_ht[observed_ht.key].variant_count
    )

    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants)) // 3
    logger.info(
        "Total bases to use when calculating correction_factors: %f", total_bases
    )

    # Get the index of dowsampling with size of `downsampling_level`.
    downsampling_idx = hl.eval(
        downsampling_meta.map(lambda x: hl.int(x[1]["downsampling"])).index(
            downsampling_level
        )
    )

    # Compute the proportion observed, which represents the relative mutability of each
    # variant class.
    ann_expr = {
        "proportion_observed": ht.variant_count / ht.possible_variants,
        f"proportion_observed_{downsampling_level}": (
            ht.downsampling_counts_global[downsampling_idx] / ht.possible_variants
        ),
        "downsamplings_frac_observed": (
            ht.downsampling_counts_global / ht.possible_variants
        ),
    }

    for pop in pops:
        pop_counts_expr = ht[f"downsampling_counts_{pop}"]
        correction_factors = ht.aggregate(
            total_mu / (hl.agg.array_sum(pop_counts_expr) / total_bases),
            _localize=False,
        )
        downsamplings_mu_expr = (
            correction_factors * pop_counts_expr / ht.possible_variants
        )
        ann_expr[f"downsamplings_mu_{'snp' if pop == 'global' else pop}"] = (
            downsamplings_mu_expr
        )
        ann_expr[f"mu_snp{'' if pop == 'global' else f'_{pop}'}"] = (
            downsamplings_mu_expr[downsampling_idx]
        )

    ht = ht.annotate(**ann_expr).checkpoint(
        new_temp_file(prefix="calculate_mu_by_downsampling", extension="ht")
    )

    ht = ht.annotate_globals(
        ac_cutoff=ac_cutoff,
        downsampling_level=downsampling_level,
        total_mu=total_mu,
        min_cov=min_cov,
        max_cov=max_cov,
        gerp_lower_cutoff=gerp_lower_cutoff,
        gerp_upper_cutoff=gerp_upper_cutoff,
    )

    return annotate_mutation_type(ht)


def compute_constraint_metrics(
    ht: hl.Table,
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Tuple[str] = (
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ),
    pops: Tuple[str] = (),
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
    raw_z_outlier_threshold_lof: float = -5.0,
    raw_z_outlier_threshold_missense: float = -5.0,
    raw_z_outlier_threshold_lower_syn: float = -5.0,
    raw_z_outlier_threshold_upper_syn: float = 5.0,
    include_os: bool = False,
) -> hl.Table:
    """
    Compute the pLI scores, observed:expected ratio, 90% confidence interval around the observed:expected ratio, and z scores for synonymous variants, missense variants, and predicted loss-of-function (pLoF) variants.

    .. note::
        The following annotations should be present in `ht`:
            - modifier
            - annotation
            - observed_variants
            - mu
            - possible_variants
            - expected_variants
            - expected_variants_{pop} (if `pops` is specified)
            - downsampling_counts_{pop} (if `pops` is specified)

    :param ht: Input Table with the number of expected variants (output of
        `get_proportion_observed()`).
    :param keys: The keys of the output Table, defaults to ('gene', 'transcript',
        'canonical').
    :param classic_lof_annotations: Classic LoF Annotations used to filter the input
        Table. Default is {"stop_gained", "splice_donor_variant",
        "splice_acceptor_variant"}.
    :param pops: List of populations used to compute constraint metrics. Default is ().
    :param expected_values: Dictionary containing the expected values for 'Null',
        'Rec', and 'LI' to use as starting values.
    :param min_diff_convergence: Minimum iteration change in LI to consider the EM
        model convergence criteria as met. Default is 0.001.
    :param raw_z_outlier_threshold_lof: Value at which the raw z-score is considered an outlier for lof variants. Values below this threshold will be considered outliers. Default is -5.0.
    :param raw_z_outlier_threshold_missense: Value at which the raw z-score is considered an outlier for missense variants. Values below this threshold will be considered outliers. Default is -5.0.
    :param raw_z_outlier_threshold_lower_syn: Lower value at which the raw z-score is considered an outlier for synonymous variants. Values below this threshold will be considered outliers. Default is -5.0.
    :param raw_z_outlier_threshold_upper_syn: Upper value at which the raw z-score is considered an outlier for synonymous variants. Values above this threshold will be considered outliers. Default is  5.0.
    :param include_os: Whether or not to include OS (other splice) as a grouping when
        stratifying calculations by lof HC.
    :return: Table with pLI scores, observed:expected ratio, confidence interval of the
        observed:expected ratio, and z scores.
    """
    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.463, "LI": 0.089}
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X.
    # `annotation_dict` stats the rule of filtration for each annotation.
    annotation_dict = {
        # Filter to classic LoF annotations with LOFTEE HC or LC.
        "lof_hc_lc": hl.literal(set(classic_lof_annotations)).contains(
            ht.annotation
        ) & ((ht.modifier == "HC") | (ht.modifier == "LC")),
        # Filter to LoF annotations with LOFTEE HC.
        "lof": ht.modifier == "HC",
        # Filter to missense variants.
        "mis": ht.annotation == "missense_variant",
        # Filter to probably damaging missense variants predicted by PolyPen-2.
        "mis_pphen": ht.modifier == "probably_damaging",
        # Filter to synonymous variants.
        "syn": ht.annotation == "synonymous_variant",
    }

    # Define two lists of 'annotation_dict' keys that require different computations.
    # The 90% CI around obs:exp and z-scores are only computed for lof, mis, and syn.
    oe_ann = ["lof", "mis", "syn"]
    # pLI scores are only computed for LoF variants.
    lof_ann = ["lof_hc_lc", "lof"]

    # Create dictionary with outlier z-score thresholds with annotation as key and list of thresholds [lower, upper] as values.
    z_score_outlier_dict = {
        "lof": [raw_z_outlier_threshold_lof, None],
        "mis": [raw_z_outlier_threshold_missense, None],
        "syn": [raw_z_outlier_threshold_lower_syn, raw_z_outlier_threshold_upper_syn],
    }

    if include_os:
        # Filter to LoF annotations with LOFTEE HC or OS.
        annotation_dict.update(
            {"lof_hc_os": (ht.modifier == "HC") | (ht.modifier == "OS")}
        )
        lof_ann.append("lof_hc_os")

    # Compute the observed:expected ratio. Will not compute per pop for "mis_pphen".
    ht = ht.group_by(*keys).aggregate(
        **{
            ann: oe_aggregation_expr(
                ht,
                filter_expr,
                pops=() if ann == "mis_pphen" else pops,
                exclude_mu_sum=True if ann == "mis_pphen" else False,
            )
            for ann, filter_expr in annotation_dict.items()
        }
    )
    # Filter to only rows with at least 1 obs or exp across all keys in annotation_dict.
    ht = ht.filter(
        hl.sum(
            [
                hl.or_else(ht[ann].obs, 0) + hl.or_else(ht[ann].exp, 0)
                for ann in annotation_dict
            ]
        )
        > 0
    )
    ht = ht.checkpoint(
        new_temp_file(prefix="compute_constraint_metrics", extension="ht")
    )

    # Compute the pLI scores for LoF variants.
    ann_expr = {
        ann: ht[ann].annotate(
            **compute_pli(
                ht,
                obs_expr=ht[ann].obs,
                exp_expr=ht[ann].exp,
                expected_values=expected_values,
                min_diff_convergence=min_diff_convergence,
            )
        )
        for ann in lof_ann
    }

    # Add a 'no_variants' flag indicating that there are zero observed variants summed
    # across pLoF, missense, and synonymous variants.
    constraint_flags_expr = {
        "no_variants": hl.sum([hl.or_else(ht[ann].obs, 0) for ann in oe_ann]) == 0
    }
    constraint_flags = {}
    for ann in oe_ann:
        obs_expr = ht[ann].obs
        exp_expr = ht[ann].exp
        # Compute the 90% confidence interval around the observed:expected ratio.
        oe_ci_expr = oe_confidence_interval(obs_expr, exp_expr)
        # Compute raw z-scores.
        raw_z_expr = calculate_raw_z_score(obs_expr, exp_expr)
        # Add flags that define why constraint will not be calculated.
        ann_constraint_flags_expr = get_constraint_flags(
            exp_expr=exp_expr,
            raw_z_expr=raw_z_expr,
            raw_z_lower_threshold=z_score_outlier_dict[ann][0],
            raw_z_upper_threshold=z_score_outlier_dict[ann][1],
            flag_postfix=ann,
        )
        constraint_flags_expr.update(ann_constraint_flags_expr)
        # The constraint_flags dict is used to filter the final ht.constraint_flags
        # annotation to the flags that should be considered in the z-score 'sd'
        # computation of the specified ann.
        constraint_flags[ann] = hl.set(
            ann_constraint_flags_expr.keys() | {"no_variants"}
        )
        # Add initial ann to ann_expr if it isn't present.
        # The ann_expr dict will already have all ann in lof_ann.
        if ann not in ann_expr:
            ann_expr[ann] = ht[ann]

        ann_expr[ann] = ann_expr[ann].annotate(
            oe_ci=oe_ci_expr,
            z_raw=raw_z_expr,
        )

    ann_expr["constraint_flags"] = add_filters_expr(filters=constraint_flags_expr)
    ht = ht.annotate(**ann_expr)
    ht = ht.checkpoint(
        new_temp_file(prefix="compute_constraint_metrics", extension="ht")
    )

    # Add z-score 'sd' annotation to globals.
    ht = ht.annotate_globals(
        sd_raw_z=ht.aggregate(
            hl.struct(
                **{
                    ann: calculate_raw_z_score_sd(
                        raw_z_expr=ht[ann].z_raw,
                        flag_expr=ht.constraint_flags.intersection(
                            constraint_flags[ann]
                        ),
                        mirror_neg_raw_z=(ann != "syn"),
                    )
                    for ann in oe_ann
                }
            )
        )
    )

    # Compute z-score from raw z-score and standard deviations.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(z_score=ht[ann].z_raw / ht.sd_raw_z[ann])
            for ann in oe_ann
        }
    )

    return ht


def calculate_gerp_cutoffs(ht: hl.Table) -> Tuple[float, float]:
    """
    Find GERP cutoffs determined by the 5% and 95% percentiles.

    :param ht: Input Table.
    :return: Tuple containing values determining the 5-95th percentile of the GERP score.
    """
    # Aggregate histogram of GERP values from -12.3 to 6.17 (-12.3 to 6.17 is the range
    # of GERP values where 6.17 is the most conserved).
    summary_hist = ht.aggregate(hl.struct(gerp=hl.agg.hist(ht.gerp, -12.3, 6.17, 100)))

    # Get cumulative sum of the hist array and add value of n_smaller to every value in
    # the cumulative sum array.
    cumulative_data = (
        np.cumsum(summary_hist.gerp.bin_freq) + summary_hist.gerp.n_smaller
    )

    # Append final value to the cumulative sum array (value added is last value of the
    # array plus n_larger).
    np.append(cumulative_data, [cumulative_data[-1] + summary_hist.gerp.n_larger])

    # Get zip of (bin_edge, value in cumulative sum array divided by max value in
    # cumulative sum array).
    zipped = zip(summary_hist.gerp.bin_edges, cumulative_data / max(cumulative_data))

    # Define lower and upper GERP cutoffs based on 5th and 95th percentiles.
    cutoff_lower = list(filter(lambda i: i[1] > 0.05, zipped))[0][0]

    zipped = zip(summary_hist.gerp.bin_edges, cumulative_data / max(cumulative_data))
    cutoff_upper = list(filter(lambda i: i[1] < 0.95, zipped))[-1][0]

    return cutoff_lower, cutoff_upper


def annotate_context_ht(
    ht: hl.Table,
    coverage_hts: Dict[str, hl.Table],
    methylation_ht: hl.Table,
    gerp_ht: hl.Table,
) -> hl.Table:
    """
    Split multiallelic sites if needed and add 'methylation', 'coverage', and 'gerp' annotation to context Table with VEP annotation.

    .. note::
        Checks for 'was_split' annotation in Table. If not present, splits
        multiallelic sites.

    :param ht: Input context Table with VEP annotation.
    :param coverage_hts: A Dictionary with key as one of 'exomes' or 'genomes' and
        values as corresponding coverage Tables.
    :param methylation_ht: Methylation Table.
    :param gerp_ht: Table with GERP annotation.
    :return: Table with sites split and necessary annotations.
    """
    # Check if context Table is split, and if not, split multiallelic sites.
    if "was_split" not in list(ht.row):
        ht = hl.split_multi_hts(ht)

    # Filter Table to only contigs 1-22, X, Y.
    ref = get_reference_genome(ht.locus)
    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(c, ref.name) for c in ref.contigs[:24]]
    )

    # If neccessary, pull out first element of coverage statistics (which includes all samples). Relevant to v4, where
    # coverage stats include additional elements to stratify by ukb subset and
    # platforms.
    if "coverage_stats" in coverage_hts["exomes"].row:
        coverage_hts["exomes"] = coverage_hts["exomes"].transmute(
            **coverage_hts["exomes"].coverage_stats[0]
        )

    # Add 'methylation', 'coverage', and 'gerp' annotation.
    ht = ht.annotate(
        methylation=methylation_ht[ht.locus],
        coverage=hl.struct(
            **{loc: coverage_ht[ht.locus] for loc, coverage_ht in coverage_hts.items()}
        ),
        gerp=gerp_ht[ht.locus].S,
    )

    ht = ht.annotate(gerp=hl.if_else(hl.is_missing(ht.gerp), 0, ht.gerp))

    return ht
