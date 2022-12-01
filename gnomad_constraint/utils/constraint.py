"""Script containing utility functions used in the constraint pipeline."""
# cSpell: disable
import pickle
from typing import Dict, List, Optional, Set, Tuple, Union

import hail as hl
from gnomad.utils.constraint import (
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    calculate_all_z_scores,
    collapse_strand,
    compute_all_pLI_scores,
    compute_expected_variants,
    compute_oe_per_transcript,
    count_variants_by_group,
    get_downsamplings,
    oe_confidence_interval,
    trimer_from_heptamer,
)
from gnomad.utils.filtering import (
    filter_for_mu,
    filter_to_autosomes,
    remove_coverage_outliers,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import COVERAGE_CUTOFF, mutation_rate_ht


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


def prepare_ht_for_constraint_calculations(ht: hl.Table) -> hl.Table:
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
    :return: Table with annotations.
    """
    ht = trimer_from_heptamer(ht)

    if "filters" in ht.row_value.keys():
        ht = ht.annotate(pass_filters=hl.len(ht.filters) == 0)

    # Add annotations for 'ref' and 'alt'
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
    # Filter to SNPs and context fields where the bases are either A, T, C, or G
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt) & ht.context.matches(f"[ATCG]{{{3}}}"))
    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and
    # collapse strands to deduplicate the context
    ht = annotate_mutation_type(collapse_strand(ht))
    # Add annotation for the methylation level
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    # Add annotation for the median exome coverage
    annotation["exome_coverage"] = ht.coverage.exomes.median
    ht = ht.annotate(**annotation)

    # Add most_severe_consequence annotation to 'transcript_consequences' within the
    # vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)

    # Filter out locus with undefined exome coverage
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
    filter_to_canonical_synonymous: bool = False,
    global_annotation: Optional[str] = None,
) -> hl.Table:
    """
    Count the observed variants and possible variants by substitution, context, methylation level, and additional `grouping`.

    Prior to computing variant counts the following variants are removed:
        - Variants not observed by any samples in the dataset: `(freq_expr.AC > 0)`
        - Low-quality variants: `exome_ht.pass_filters`
        - Variants with allele frequency above `max_af` cutoff: `(freq_expr.AF <=
          max_af)`
        - Variants that are not synonymous or in the canonical transcript

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
    :param filter_to_canonical_synonymous: Whether to keep only canonical synonymous
        variants in the exome Table. Default is False.
    :param global_annotation: The annotation name to use as a global StructExpression annotation containing
        input parameter values. If no value is supplied, this global annotation will not be added. Default is None.
    :return: Table with observed variant and possible variant count.
    """
    # Allele frequency information for high-quality genotypes (GQ >= 20; DP >= 10; and
    # AB >= 0.2 for heterozygous calls) in all release samples in gnomAD.
    freq_expr = exome_ht.freq[0]

    # Set up the criteria to exclude variants not observed in the dataset, low-quality
    # variants, variants with allele frequency above the `max_af` cutoff, and variants
    # with exome coverage larger then 0 if requested.
    keep_criteria = (
        (freq_expr.AC > 0) & exome_ht.pass_filters & (freq_expr.AF <= max_af)
    )
    if filter_coverage_over_0:
        keep_criteria &= exome_ht.coverage > 0

    keep_annotations += grouping

    # Keep variants that satisfy the criteria above.
    filtered_exome_ht = exome_ht.filter(keep_criteria)

    # If requested keep only variants that are synonymous, canonical
    if filter_to_canonical_synonymous:
        filtered_exome_ht = filter_vep_transcript_csqs(exome_ht.filter(keep_criteria))
        context_ht = filter_vep_transcript_csqs(context_ht)
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `keep_annotations`.
    observed_ht = count_variants_by_group(
        filtered_exome_ht.select(*list(keep_annotations) + ["freq"]),
        additional_grouping=grouping,
        partition_hint=partition_hint,
        count_downsamplings=pops,
        use_table_group_by=True,
    )
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
    coverage_model: Tuple[float, float],
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    pops: Tuple[str] = (),
    partition_hint_for_counting_variants: int = 2000,
    partition_hint_for_aggregation: int = 1000,
    custom_vep_annotation: str = None,
    cov_cutoff: int = COVERAGE_CUTOFF,
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
    :param partition_hint_for_counting_variants: Target number of partitions for
        aggregation when counting variants. Default is 2000.
    :param partition_hint_for_aggregation: Target number of partitions for sum
        aggregators when computation is done. Default is 1000.
    :param custom_vep_annotation: The customized model (one of
        "transcript_consequences" or "worst_csq_by_gene"), Default is None.
    :param cov_cutoff: Median coverage cutoff. Sites with coverage above this cutoff
        are considered well covered and was used to build plateau models. Sites
        below this cutoff have low coverage and was used to build coverage models.
        Default is `COVERAGE_CUTOFF`.
    :return: Table with `expected_variants` (expected variant counts) and `obs_exp`
        (observed:expected ratio) annotations.
    """
    # Add necessary constraint annotations for grouping
    if custom_vep_annotation == "worst_csq_by_gene":
        vep_annotation = "worst_csq_by_gene"
    else:
        vep_annotation = "transcript_consequences"

    context_ht, _ = annotate_exploded_vep_for_constraint_groupings(
        context_ht, vep_annotation
    )
    exome_ht, grouping = annotate_exploded_vep_for_constraint_groupings(
        exome_ht, vep_annotation
    )

    # Compute observed and possible variant counts
    ht = create_observed_and_possible_ht(
        exome_ht,
        context_ht,
        mutation_ht,
        max_af,
        keep_annotations,
        pops,
        grouping,
        partition_hint_for_counting_variants,
        filter_coverage_over_0=True,
        filter_to_canonical_synonymous=False,
    )

    mu_expr = ht.mu_snp * ht.possible_variants
    # Determine coverage correction to use based on coverage value.
    cov_corr_expr = (
        hl.case()
        .when(ht.coverage == 0, 0)
        .when(ht.coverage >= cov_cutoff, 1)
        .default(coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0])
    )
    # Generate sum aggregators for 'mu' on the entire dataset.
    agg_expr = {"mu": hl.agg.sum(mu_expr * cov_corr_expr)}
    agg_expr.update(
        compute_expected_variants(ht, plateau_models, mu_expr, cov_corr_expr)
    )
    for pop in pops:
        agg_expr.update(
            compute_expected_variants(ht, plateau_models, mu_expr, cov_corr_expr, pop)
        )

    grouping = list(grouping)
    grouping.remove("coverage")

    # Aggregate the sum aggregators grouped by `grouping`.
    ht = (
        ht.group_by(*grouping)
        .partition_hint(partition_hint_for_aggregation)
        .aggregate(**agg_expr)
    )

    # Annotate global annotations.
    ht = ht.annotate_globals(
        apply_model_params=hl.struct(
            max_af=max_af,
            pops=pops,
            plateau_models=plateau_models,
            coverage_model=coverage_model,
        )
    )
    # Compute the observed:expected ratio.
    return ht.annotate(obs_exp=ht.observed_variants / ht.expected_variants)


def calculate_mu_by_downsampling(
    genome_ht: hl.Table,
    raw_context_ht: hl.Table,
    recalculate_all_possible_summary: bool = True,
    recalculate_all_possible_summary_unfiltered: bool = False,
    omit_methylation: bool = False,
    count_singletons: bool = False,
    summary_file: str = None,
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
) -> hl.Table:
    """
    Calculate mutation rate using the downsampling with size specified by `downsampling_level` in genome sites Table.

    Prior to computing mutation rate the only following variants are kept:
        - variants with the mean coverage in the gnomAD genomes was between 15X and 60X.
        - variants whoes the most severe annotation was intron_variant or
            intergenic_variant
        - variants with the GERP score between the 5th and 95th percentile of the
            genomewide distribution.
        - high-quality variants: `exome_ht.pass_filters`
        - Variants with allele count below `ac_cutoff`: `(freq_expr.AC <= ac_cutoff)`

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation_level
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param genome_ht: Genome sites Table.
    :param raw_context_ht: Context Table with locus that is on an autosome or in a
        pseudoautosomal region and sex chromosomes.
    :param recalculate_all_possible_summary: Whether to calculate possible variants
        using context Table with locus that is only on an autosome or in a
        pseudoautosomal region. Default is True.
    :param recalculate_all_possible_summary_unfiltered: Whether to calculate possible
        variants using raw context Table. Default is True.
    :param omit_methylation: Whether to omit 'methylation_level' from the grouping when
        counting variants. Default is False.
    :param count_singletons: Whether to count singletons. Default is False.
    :param summary_file: _description_, defaults to None
    :param keep_annotations: Annotations to keep in the context Table and genome sites
        Table.
    :param ac_cutoff: The cutoff of allele count when filtering context Table and genome sites Table.
    :param downsampling_level: The size of downsamplings will be used to count variants. Default is 1000.
    :param total_mu: The per-generation mutation rate. Default is 1.2e-08.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :return: Mutation rate Table.
    """
    # keep only loci where the mean coverage in the gnomAD genomes was between 15X and
    # 60X.
    context_ht = filter_to_autosomes(remove_coverage_outliers(raw_context_ht))
    genome_ht = filter_to_autosomes(remove_coverage_outliers(genome_ht))

    # filter the Table so that the most severe annotation was intron_variant or
    # intergenic_variant, and that the GERP score was between the 5th and 95th
    # percentile of the genomewide distribution.
    context_ht = filter_for_mu(context_ht)
    genome_ht = filter_for_mu(genome_ht)

    context_ht = context_ht.select(*keep_annotations)
    genome_ht = genome_ht.select(*list(keep_annotations) + ["freq", "pass_filters"])

    # downsampled the dataset to 1,000 genomes
    freq_index = hl.eval(
        hl.enumerate(genome_ht.freq_meta).find(
            lambda f: (f[1].get("downsampling") == str(downsampling_level))
            & (f[1].get("pop") == "global")
            & (f[1].get("group") == "adj")
            & (f[1].size() == 3)
        )
    )[0]
    freq_expr = genome_ht.freq[freq_index]
    # Set up the criteria to filtered out  low quality sites, and sites found in
    # greater than 5 copies in the downsampled set.
    keep_criteria = (freq_expr.AC <= ac_cutoff) & genome_ht.pass_filters

    # Only keep variants with AC <= ac_cutoff and passing filters in genome site
    # Table.
    genome_ht = genome_ht.filter(keep_criteria)

    # In the denominator, only keep variants not in the genome dataset, or with AC
    # <= ac_cutoff and passing filters.
    context_ht = context_ht.anti_join(genome_ht.filter(keep_criteria, keep=False))

    if not summary_file:
        summary_file = new_temp_file(prefix="constraint", extension="he")

    all_possible_dtype = count_variants_by_group(
        context_ht,
        omit_methylation=omit_methylation,
        return_type_only=True,
    ).variant_count.dtype

    # Count possible variants in context Table.
    if recalculate_all_possible_summary:
        all_possible = count_variants_by_group(
            context_ht,
            omit_methylation=omit_methylation,
        ).variant_count
        # with hl.hadoop_open(summary_file, "wb") as f:
        #     pickle.dump(all_possible, f)
        hl.experimental.write_expression(all_possible, summary_file)
    # with hl.hadoop_open(summary_file, "rb") as f:
    #     all_possible = pickle.load(f)
    all_possible = hl.experimental.read_expression(summary_file)

    if recalculate_all_possible_summary_unfiltered:
        all_possible_unfiltered = count_variants_by_group(
            raw_context_ht.rows(), omit_methylation=omit_methylation
        ).variant_count
        all_possible_unfiltered = {
            x: y
            for x, y in list(all_possible_unfiltered.items())
            if x.context is not None
        }
        hl.experimental.write_expression(
            all_possible_unfiltered, new_temp_file(prefix="constraint", extension="he")
        )
        # with hl.hadoop_open(new_temp_file(prefix="constraint", extension="he"), "wb") as f:
        #     pickle.dump(all_possible_unfiltered, f)
    # Uncomment this line and next few commented lines to back-calculate total_mu from the old mutation rate dataset
    # all_possible_unfiltered = load_all_possible_summary(filtered=False)

    # Count the observed variants in the genome sites Table.
    genome_ht = count_variants_by_group(
        genome_ht,
        count_downsamplings=pops,
        count_singletons=count_singletons,
        omit_methylation=omit_methylation,
    )

    old_mu_data = mutation_rate_ht.ht()
    ht = genome_ht.annotate(
        possible_variants=hl.literal(all_possible, dtype=all_possible_dtype)[
            genome_ht.key
        ],
        # possible_variants_unfiltered=hl.literal(all_possible_unfiltered, dtype=all_possible_dtype)[genome_ht.key],
        old_mu_snp=old_mu_data[
            hl.struct(context=genome_ht.context, ref=genome_ht.ref, alt=genome_ht.alt)
        ].mu_snp,
    )

    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants)) // 3
    # total_bases_unfiltered = ht.aggregate(hl.agg.sum(ht.possible_variants_unfiltered)) // 3
    # total_mu = ht.aggregate(hl.agg.sum(ht.old_mu_snp * ht.possible_variants_unfiltered) / total_bases_unfiltered)

    # Get the index of dowsamplings with size of 1000 genomes
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht.freq_meta)))
    downsampling_idx = downsamplings.index(1000)

    for pop in pops:
        correction_factors = ht.aggregate(
            total_mu
            / (hl.agg.array_sum(ht[f"downsampling_counts_{pop}"]) / total_bases)
        )
        ht = ht.annotate(
            **{
                f"downsamplings_mu_{pop}": hl.literal(correction_factors)
                * ht[f"downsampling_counts_{pop}"]
                / ht.possible_variants
            }
        )
        if pop == "global":
            ht = ht.annotate(
                **{"mu_snp": ht[f"downsamplings_mu_{pop}"][downsampling_idx]}
            )
        else:
            ht = ht.annotate(
                **{f"mu_snp_{pop}": ht[f"downsamplings_mu_{pop}"][downsampling_idx]}
            )

    ht = annotate_mutation_type(
        ht.annotate(
            downsamplings_frac_observed=ht.downsampling_counts_global
            / ht.possible_variants
        )
    )

    # Compute the proportion observed, which represents the relative mutability of each
    # variant class
    return ht.annotate(
        proportion_observed_1kg=ht.downsampling_counts_global[downsampling_idx]
        / ht.possible_variants,
        proportion_observed=ht.variant_count / ht.possible_variants,
    )


def compute_constraint_metrics(
    ht: hl.Table,
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Set[str] = {
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    },
    pops: Tuple[str] = (),
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
    :param classic_lof_annotations: Classic LoF Annotations used to filter Input Table. Default is {}.
    :param pops: List of populations used to compute constraint metrics. Default is ().
    :return: Table with pLI scores, observed:expected ratio, confidence interval of the
        observed:expected ratio, and z scores.
    """
    classic_lof_annotations = hl.literal(classic_lof_annotations)
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X.
    # `annotation_dict` stats the rule of filtration for each annotation.
    ht_annotations = ["lof_hc_lc", "lof_hc_os", "lof", "mis", "mis_pphen", "syn"]
    annotation_dict = {
        # Filter to classic LoF annotations with LOFTEE HC or LC.
        "lof_hc_lc": classic_lof_annotations.contains(ht.annotation)
        & ((ht.modifier == "HC") | (ht.modifier == "LC")),
        # Filter to LoF annotations with LOFTEE HC or OS.
        "lof_hc_os": (ht.modifier == "HC") | (ht.modifier == "OS"),
        # Filter to LoF annotations with LOFTEE HC.
        "lof": ht.modifier == "HC",
        # Filter to missense variants.
        "mis": ht.annotation == "missense_variant",
        # Filter to probably damaging missense variants predicted by PolyPen-2.
        "mis_pphen": ht.modifier == "probably_damaging",
        # Filter to synonymous variants.
        "syn": ht.annotation == "synonymous_variant",
    }

    all_ht = None
    for annotation_name in ht_annotations:
        # Compute the observed: expected ratio.
        if annotation_name != "mis_pphen":
            ht = compute_oe_per_transcript(
                ht, annotation_name, annotation_dict[annotation_name], keys, pops
            )
        else:
            ht = compute_oe_per_transcript(
                ht, annotation_name, annotation_dict[annotation_name], keys
            )

        # Compute the pLI scores for LoF variants.
        if "lof" in annotation_name:
            ht = compute_all_pLI_scores(ht, keys, annotation_name)

        # Compute the 90% confidence interval around the observed:expected ratio and z
        # scores.
        if annotation_name in ["lof", "mis", "syn"]:
            oe_ci = oe_confidence_interval(
                ht,
                ht[f"obs_{annotation_name}"],
                ht[f"exp_{annotation_name}"],
                prefix=f"oe_{annotation_name}",
            )
            ht = ht.annotate(**oe_ci[ht.key])
            ht = calculate_all_z_scores(ht)

        # Combine all the metrics annotations.
        if not all_ht:
            all_ht = ht
        else:
            all_ht = all_ht.annotate(**ht[all_ht.key])

    return all_ht


def split_context_ht(
    vep_context_ht: str,
    coverage_hts: Dict[str, hl.Table],
    methylation_ht: str,
    gerp_ht: str,
) -> hl.Table:
    """
    Split multiallelic sites and add 'methylation', 'coverage', and 'gerp' annotation to context Table with VEP annotation.

    :param vep_context_ht: Input Table with VEP annotation.
    :param coverage_hts: A Dictionary with key as one of 'exomes' or 'genomes' and
        values as corresponding coverage Tables.
    :param methylation_ht: Methylation Table.
    :param gerp_ht: Table with gerp annotation.
    :return: Table with splitted sites and necessary annotations.
    """
    # Split multiallelic sites in context Table.
    split_context_ht = hl.split_multi_hts(vep_context_ht)
    # raw_context_ht = raw_context_ht.checkpoint(f'{raw_context_ht_path}.temp_split.ht', overwrite)

    # Drop unnecessary annotations in coverage Table.
    coverage_hts = {
        loc: coverage_ht.drop("#chrom", "pos")
        if "#chrom" in list(coverage_ht.row)
        else coverage_ht
        for loc, coverage_ht in coverage_hts.items()
    }

    # Add 'methylation', 'coverage', and 'gerp' annotation.
    annotated_context_ht = split_context_ht.annotate(
        methylation=methylation_ht[split_context_ht.locus],
        coverage=hl.struct(
            **{
                loc: coverage_ht[split_context_ht.locus]
                for loc, coverage_ht in coverage_hts.items()
            }
        ),
        gerp=gerp_ht[split_context_ht.locus].S,
    )
    return split_context_ht, annotated_context_ht
