"""Script containing utility functions used in the constraint pipeline."""

from functools import reduce
from typing import Dict, Optional, Set, Tuple, Union

import hail as hl
from gnomad.utils.constraint import (
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    calculate_z_score,
    collapse_strand,
    compute_expected_variants,
    compute_pli,
    count_variants_by_group,
    oe_aggregation_expr,
    oe_confidence_interval,
    trimer_from_heptamer,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import COVERAGE_CUTOFF


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
    :param global_annotation: The annotation name to use as a global StructExpression
        annotation containing input parameter values. If no value is supplied, this
        global annotation will not be added. Default is None.
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
    obs_pos_count_partition_hint: int = 2000,
    expected_variant_partition_hint: int = 1000,
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
    :param obs_pos_count_partition_hint: Target number of partitions for
        aggregation when counting variants. Default is 2000.
    :param expected_variant_partition_hint: Target number of partitions for sum
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
        obs_pos_count_partition_hint,
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


def compute_constraint_metrics(
    ht: hl.Table,
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Set[str] = {
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    },
    pops: Tuple[str] = (),
    expected_values: Dict[str, float] = {"Null": 1.0, "Rec": 0.463, "LI": 0.089},
    min_diff_convergence: float = 0.001,
    raw_z_outlier_threshold: int = 5,
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
    :param raw_z_outlier_threshold: Value at which the raw z-score is considered an outlier. Values below the negative of '--raw-z-outlier-threshold' will be considered outliers for lof and missense varaint counts (indicating too many variants), whereas values either above '--raw-z-outlier-threshold' or below the negative of '--raw-z-outlier-threshold' will be considered outliers for synonymous varaint counts (indicating too few or too many variants)."
    :return: Table with pLI scores, observed:expected ratio, confidence interval of the
        observed:expected ratio, and z scores.
    """
    classic_lof_annotations = hl.literal(classic_lof_annotations)
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X.
    # `annotation_dict` stats the rule of filtration for each annotation.
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
    ht = ht.checkpoint(
        new_temp_file(prefix="compute_constraint_metrics", extension="ht")
    )

    # Compute the pLI scores for LoF variants.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(**compute_pli(ht, ht[ann].obs, ht[ann].exp))
            for ann in ["lof_hc_lc", "lof_hc_os", "lof"]
        }
    )

    # Compute the 90% confidence interval around the observed:expected ratio and z
    # scores.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(
                oe_ci=oe_confidence_interval(ht[ann].obs, ht[ann].exp),
                z_raw=calculate_raw_z_score(ht[ann].obs, ht[ann].exp),
            )
            for ann in ["lof", "mis", "syn"]
        }
    )

    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(
                constraint_flags=add_constraint_flags(
                    exp_expr=ht[ann].exp,
                    outlier_expr=ht[ann].z_raw < -raw_z_outlier_threshold
                    if ann != "syn"
                    else hl.abs(ht[ann].z_raw) > raw_z_outlier_threshold,
                )
            )
            for ann in ["lof", "mis", "syn"]
        }
    )

    # Add 'no_variants' to the constraint flags if there are a total of 0 observed
    # variants summed across lof, mis, and syn.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(
                constraint_flags=hl.if_else(
                    hl.or_else(ht.lof["obs"], 0)
                    + hl.or_else(ht.mis["obs"], 0)
                    + hl.or_else(ht.syn["obs"], 0)
                    == 0,
                    ht[ann].constraint_flags.add("no_variants"),
                    ht[ann].constraint_flags,
                )
            )
            for ann in ["lof", "mis", "syn"]
        }
    )

    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(
                **calculate_z_score(
                    raw_z_expr=ht[ann].z_raw,
                    flag_expr=ht[ann].constraint_flags,
                    additional_requirements_expr=ht[ann].z_raw < 0
                    if ann != "syn"
                    else True,
                    both=True if ann != "syn" else False,
                )
            )
            for ann in ["lof", "mis", "syn"]
        }
    )

    flag_lists = [
        ht[ann].constraint_flags.map(
            lambda x: hl.if_else(x != "no_variants", (x + "_" + ann), x)
        )
        for ann in ["lof_hc_lc", "new_lof_hc_lc"]
    ]
    ht = ht.annotate(all_constraint_flags=flag_lists[0].union(*flag_lists[1:]))

    # TODO: Move standard deviations into globals.

    return ht
