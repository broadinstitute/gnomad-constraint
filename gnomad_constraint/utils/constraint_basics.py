# noqa: D100
# cSpell: disable

from typing import Optional, Tuple

import hail as hl

from gnomad.utils.constraint import (
    trimer_from_heptamer,
    annotate_mutation_type,
    collapse_strand,
    annotate_with_mu,
    count_variants_by_group,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)

from .generic import (
    fast_filter_vep,
    count_variants,
    annotate_with_mu,
    remove_unnecessary_variants,
    get_all_pop_lengths,
)

POPS = ("global", "afr", "amr", "eas", "nfe", "sas")
HIGH_COVERAGE_CUTOFF = 40


def add_vep_context_annotations(ht: hl.Table, annotated_context_ht: str) -> hl.Table:
    """
    Add annotations from VEP context Table to gnomAD data.

    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp
        - pass_filters

    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from gnomAD data.

    .. note::
        Function expects that multiallelic variants in the VEP context Table have been split.

    :param ht: gnomAD exomes or genomes public Hail Table.
    :param annotated_context_ht: VEP context Table.
    """
    context_ht = annotated_context_ht.drop("a_index", "was_split")
    context_ht = context_ht.annotate(vep=context_ht.vep.drop("colocated_variants"))
    return ht.annotate(**context_ht[ht.key])


def prepare_ht_for_constraint_calculations(ht: hl.Table):
    """
    Filter input Table and add annotations used in constraint calculations.

    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and `add_most_severe_csq_to_tc_within_vep_root()`
        - exome_coverage (optionally added)

    :param ht: Input Table to be annotated.
    :param annotate_coverage: Whether to annotate the coverage of exome. Defaults to True.
    :return: Table with annotations.
    """
    ht = trimer_from_heptamer(ht)

    if "filters" in ht.row_value.keys():
        ht = ht.annotate(pass_filters=hl.len(ht.filters) == 0)

    # Add annotations for 'ref' and 'alt'
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
    # Filter to SNPs and context fields where the bases are either A, T, C, or G
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt) & ht.context.matches(f"[ATCG]{{{3}}}"))
    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and collapse strands to deduplicate the context
    ht = annotate_mutation_type(collapse_strand(ht))
    # Add annotation for the methylation level
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    # Add annotation for the median exome coverage if requested
    annotation["exome_coverage"] = ht.coverage.exomes.median
    ht = ht.annotate(**annotation)
    # Add most_severe_consequence annotation to 'transcript_consequences' within the vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)

    # Filter out locus with undefined exome coverage
    ht = ht.filter(hl.is_defined(ht.exome_coverage))
    return ht


def create_constraint_training_dataset(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
    ),
    pops: Optional[Tuple[str]] = (),
    grouping: Tuple[str] = ("exome_coverage",),
    partition_hint: int = 100,
) -> hl.Table:
    """
    Count the observed variants and possible variants by substitution, context, methylation level, and additional `grouping`.

    Variants where there was no coverage, a low-quality variant was observed, or a variant
    above 0.1% frequency was observed  in `exome_ht` and `context_ht` was removed. For each substitution,
    context, and methylation level, observed variants were counted using `exome_ht`, and possible variants
    was counted using `context_ht` that anti joins `exome_ht`.

    The returned Table includes following annotations:
        - context - trinucleotide genomic context
        - ref - the middle base of `context`
        - alt - the alternate base
        - methylation_level - methylation_level
        - observed_variants - observed variant counts in `exome_ht`
        - possible_variants - possible variant counts in `context_ht`
        - downsampling_counts_{pop} - variant counts in downsaplings for populations in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param max_af: Maximum allele frequency for a variant to be included in returned counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table.
    :param pops: List of populations for choosing downsamplings when counting variants. Default is ().
    :param grouping: Annotations other than 'context', 'ref', 'alt' to group by when counting variants. Default is ('exome_coverage',).
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :return: Table with observed variant and possible variant count.
    """
    # It's adjusted allele frequency information for all release samples in gnomAD
    freq_expr = exome_ht.freq[0]

    keep_criteria = (
        (freq_expr.AC > 0) & exome_ht.pass_filters & (freq_expr.AF <= max_af)
    )

    ht = count_variants_by_group(
        filter_vep_transcript_csqs(exome_ht.filter(keep_criteria)).select(
            *list(keep_annotations) + ["freq"]
        ),
        additional_grouping=grouping,
        partition_hint=partition_hint,
        count_downsamplings=pops,
        use_table_group_by=True,
    )
    ht = ht.transmute(observed_variants=ht.variant_count)

    context_ht = filter_vep_transcript_csqs(context_ht).select(*keep_annotations)
    context_ht = context_ht.anti_join(exome_ht.filter(keep_criteria, keep=False))
    possible_ht = count_variants_by_group(
        context_ht, additional_grouping=grouping, use_table_group_by=True
    )
    possible_ht = annotate_with_mu(possible_ht, mutation_ht)
    possible_ht = possible_ht.transmute(possible_variants=possible_ht.variant_count)

    ht = ht.join(possible_ht, "outer")
    ht = annotate_mutation_type(ht)
    return ht

def get_proportion_observed_by_coverage(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    recompute_possible: bool = False,
    dataset: str = "gnomad",
    impose_high_af_cutoff_upfront: bool = True,
) -> hl.Table:
    """
    Count the observed variants and possible variants by exome coverage.

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param recompute_possible: Whether to use context Table to recompute the number of possible variants instead of using a precomputed intermediate Table if it exists. Defaults to False.
    :param dataset: Dataset to use when computing frequency index. Defaults to 'gnomad'.
    :param impose_high_af_cutoff_upfront: Whether to remove high frequency alleles. Defaults to True.
    :return: Table with observed variant and possible variant count.
    """
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    context_ht = fast_filter_vep(context_ht).select(
        "context", "ref", "alt", "methylation_level", "exome_coverage"
    )
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    exome_ht = fast_filter_vep(exome_ht).select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
        "freq",
        "pass_filters",
    )

    grouping = ("exome_coverage",)
    af_cutoff = 0.001

    exome_join = exome_ht[context_ht.key]
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters
        if impose_high_af_cutoff_upfront:
            crit &= ht.freq[freq_index].AF <= af_cutoff
        return crit

    context_ht = context_ht.filter(
        hl.is_missing(exome_join) | keep_criteria(exome_join)
    )

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    possible_file = f"{root}/tmp/possible_coverage.ht"
    if recompute_possible:
        possible_ht = count_variants(
            context_ht, additional_grouping=grouping, force_grouping=True
        )
        possible_ht = annotate_with_mu(possible_ht, mutation_ht)
        possible_ht.transmute(possible_variants=possible_ht.variant_count).write(
            possible_file, True
        )

    possible_ht = hl.read_table(possible_file)
    ht = count_variants(
        exome_ht,
        additional_grouping=grouping,
        partition_hint=100,
        count_downsamplings=POPS,
        impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
    )
    ht = ht.join(possible_ht, "outer")
    ht = annotate_variant_types(ht)
    return ht

def get_proportion_observed_by_coverage(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    possible_file_path: str,
    recompute_possible: bool = False,
    impose_high_af_cutoff_upfront: bool = True,
    dataset: str = "gnomad",
    af_cutoff=0.001,
    kept_context_annotations: List[str] = [
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
    ],
    kept_exome_annotations: List[str] = [
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
        "freq",
        "pass_filters",
    ],
    pops: Tuple[str] = POPS,
    grouping: Tuple[str] = ("exome_coverage"),
) -> hl.Table:
    """
    Count the observed variants and possible variants by exome coverage.

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param possible_file_path: Path to save table with possible variants.
    :param recompute_possible: Whether to use context Table to recompute the number of possible variants instead of using a precomputed intermediate Table if it exists. Defaults to False.
    :param impose_high_af_cutoff_upfront: Whether to remove high frequency alleles, defaults to True.
    :param dataset: Dataset to use when computing frequency index. Defaults to 'gnomad'.
    :param af_cutoff: Variants with AF above than AF cutoff will be removed, defaults to 0.001.
    :param kept_context_annotations: Annotations to keep in the context Table.
    :param kept_exome_annotations: Annotations to keep in the exome Table.
    :param pops: List of populations. Defaults to `POPS`.
    :param grouping: Annotations other than 'context', 'ref', 'alt' to group by when counting variants, defaults to 'exome_coverage'.
    :return: Table with observed variant and possible variant count.
    """
    context_ht = fast_filter_vep(context_ht).select(kept_context_annotations)
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))
    exome_ht = fast_filter_vep(exome_ht).select(kept_exome_annotations)

    context_ht, exome_ht = remove_unnecessary_variants(
        context_ht, exome_ht, dataset, af_cutoff, impose_high_af_cutoff_upfront
    )

    if recompute_possible:
        possible_ht = count_variants(
            context_ht, additional_grouping=grouping, force_grouping=True
        )
        possible_ht = annotate_with_mu(possible_ht, mutation_ht)
        possible_ht.transmute(possible_variants=possible_ht.variant_count).write(
            possible_file_path, True
        )

    possible_ht = hl.read_table(possible_file_path)
    ht = count_variants(
        exome_ht,
        additional_grouping=grouping,
        partition_hint=100,
        count_downsamplings=pops,
        impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
    )
    ht = ht.join(possible_ht, "outer")
    ht = annotate_variant_types(ht)
    return ht


def build_models(
    coverage_ht: hl.Table,
    trimers: bool = False,
    weighted: bool = False,
    pop_anal: bool = True,
    pops: Tuple[str] = POPS,
    cov_cutoff: int = HIGH_COVERAGE_CUTOFF,
) -> Tuple[Tuple[float, float], Dict[str, Tuple[float, float]]]:
    """
    Build coverage model and plateau models.

    This function builds plateau models to calibrate mutation rate estimates against the proportion observed
    of each substitution, context, and methylation level in `coverage_ht` considering only high coverage sites,
    or sites above a median coverage of `HIGH_COVERAGE_CUTOFF` (or half of `HIGH_COVERAGE_CUTOFF` if `half_cutoff`
    is True). If using the output of `get_proportion_observed_by_coverage` as `coverage_ht`, the proportion observed
    will be high-quality variants below 0.1% frequency at synonymous sites. Two plateau models are fit, one for CpG
    transitions and one for the remainder of sites (transversions and non-CpG transitions).

    The x and y of the plateau models:
        x: `mu_snp` - mutation rate
        y: proportion observed ('variant_count' or 'observed_{pop}' / 'possible_variants')

    For low coverage sites, or sites below `HIGH_COVERAGE_CUTOFF` (or half of `HIGH_COVERAGE_CUTOFF` if `half_cutoff` is
    True), this function performs a base-level resolution rather than exon-level to compute a coverage correction factor
    to reduce the inaccuracy of expected variant counts caused by low coverage on each base. The coverage models are built
    by first defining a metric that is derived by dividing the number of observed variants with the total number of
    possible variants times the mutation rate summed across all substitutions, contexts, and methylation level. If using
    the output of `get_proportion_observed_by_coverage` as `coverage_ht`, the number of observed variants and possible
    variants will be at synonymous sites. The function computes this metric for high coverage sites as a global scaling
    factor, and divides this metric at low coverage sites by this scaling factor to create an observed:expected ratio. Then the
    coverage model is built of log10(coverage) to this scaled ratio as a correction factor for low coverage sites.

    The x and y of the coverage model:
        x: log10('exome_coverage') at low coverage site
        y: sum('variant_count')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage site
            where `high_coverage_scale_factor` = sum('variant_count') / sum('possible_variants' * 'mu_snp') at high coverage site

    .. note::
        This function expects that the input `coverage_ht` is the output of `get_proportion_observed_by_coverage`, and
        therefore the following fields should be present in `coverage_ht`:
            - context - trinucleotide genomic context
            - ref - the middle base of `context`
            - alt - the alternate base
            - methylation_level - methylation level
            - exome_coverage - median exome coverage at integer values between 1-100
            - variant_count - the number of observed variants in the dataset for each substitution (`alt`), `context`,
            `methylation_level`, and median exome coverage (`exome_coverage`)
            - downsampling_counts_{pop} (pop defaults to `POPS`) - array of observed variant counts per population after downsampling
            - mu_snp - mutation rate
            - possible_variants - the number of possible variants in the dataset for each substitution (`alt`), `context`,
            `methylation_level`, and median exome coverage (`exome_coverage`)

    :param coverage_ht: Input coverage Table.
    :param trimers: Whether the contexts were trimmed or not. Defaults to False.
    :param weighted: Whether to weight the high coverage model (a linear regression model) by 'possible_variants'. Defaults to False.
    :param pop_anal: Whether to build models for each population. Defaults to True.
    :param pops: List of populations. Defaults to `POPS`.
    :param cov_cutoff: A coverage cutoff. Coverage of sites above the cutoff will be used to build plateau models. Otherwise they will be used to build coverage models. defaults to `HIGH_COVERAGE_CUTOFF`.
    :return: Coverage model and plateau models.
    """
    keys = ["context", "ref", "alt", "methylation_level", "mu_snp"]

    all_high_coverage_ht = coverage_ht.filter(coverage_ht.exome_coverage >= cov_cutoff)
    agg_expr = {
        "observed_variants": hl.agg.sum(all_high_coverage_ht.variant_count),
        "possible_variants": hl.agg.sum(all_high_coverage_ht.possible_variants),
    }
    if pop_anal:
        for pop in pops:
            agg_expr[f"observed_{pop}"] = hl.agg.array_sum(
                all_high_coverage_ht[f"downsampling_counts_{pop}"]
            )
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(**agg_expr)

    high_coverage_ht = annotate_variant_types(high_coverage_ht, not trimers)
    plateau_models = build_plateau_models(high_coverage_ht, weighted=weighted)

    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.variant_count)
        / hl.agg.sum(
            all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp
        )
    )

    all_low_coverage_ht = coverage_ht.filter(
        (coverage_ht.exome_coverage < cov_cutoff) & (coverage_ht.exome_coverage > 0)
    )

    low_coverage_ht = all_low_coverage_ht.group_by(
        log_coverage=hl.log10(all_low_coverage_ht.exome_coverage)
    ).aggregate(
        low_coverage_obs_exp=hl.agg.sum(all_low_coverage_ht.variant_count)
        / (
            high_coverage_scale_factor
            * hl.agg.sum(
                all_low_coverage_ht.possible_variants * all_low_coverage_ht.mu_snp
            )
        )
    )
    coverage_model = build_coverage_model(low_coverage_ht)
    # TODO: consider weighting here as well

    return coverage_model, plateau_models


def build_plateau_models(
    ht: hl.Table,
    weighted: bool = False,
    pop_anal: bool = True,
    pops: Tuple[str] = POPS,
) -> Dict[str, Tuple[float, float]]:
    """
    Calibrate high coverage model (returns intercept and slope).

    The function fits two models, one for CpG transitions and one for the remainder of sites, to calibrate
    from the mutation rate to proportion observed in total and in each population.

    .. note::
        The following annotations should be present in `ht`:
            - observed_variants - observed variant counts for each combination of 'context', 'ref', 'alt', 'methylation_level', 'mu_snp
            - observed_variants_{pop} (where pop is each population) - observed variant counts for each population
            - possible_variants - possible variant counts for each combination of 'context', 'ref', 'alt', 'methylation_level', 'mu_snp
            - cpg - whether it's a cpg site or not
            - mu_snp - mutation rate

    :param ht: High coverage Table.
    :param pop_anal: Whether to build models for each population. Defaults to True.
    :param pops: List of populations. Defaults to `POPS`.
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants'.
    :return: A Dictionary of intercepts and slopes for the full plateau models.
    """
    pop_lengths = get_all_pop_lengths(ht, pops)

    plateau_models = build_plateau_models_total(ht, weighted=weighted)

    if pop_anal:
        arg_expr = {}
        for length, pop in pop_lengths:
            plateau_model = build_plateau_models_pop(ht, pop, length, weighted=weighted)
            arg_expr[pop] = plateau_model[pop]
        plateau_models = plateau_models.annotate(**arg_expr)

    return plateau_models


def build_plateau_models_pop(
    ht: hl.Table, pop, length, weighted: bool = False
) -> Dict[str, Tuple[float, float]]:
    """
    Calibrate high coverage model for each population in `POPS`.

    :param ht: High coverage Table.
    :param pop: List of populations. Defaults to `POPS`.
    :param length: The indices of each population.
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants', defaults to False.
    :return: A Dictionary of intercepts and slopes for plateau models of each population.
    """
    agg_expr = {
        pop: [
            hl.agg.group_by(
                ht.cpg,
                hl.agg.linreg(
                    ht[f"observed_{pop}"][i] / ht.possible_variants,
                    [1, ht.mu_snp],
                    weight=ht.possible_variants if weighted else None,
                ).beta,
            )
            for i in range(length)
        ]
    }
    return ht.aggregate(hl.struct(**agg_expr))


def build_plateau_models_total(
    ht: hl.Table, weighted: bool = False
) -> Dict[str, Tuple[float, float]]:
    """
    Calibrate high coverage model in total.

    :param ht: High coverage Table.
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants', defaults to False.
    :return: A Dictionary of intercepts and slopes for a plateau model.
    """
    agg_expr = {
        "total": hl.agg.group_by(
            ht.cpg,
            hl.agg.linreg(
                ht.observed_variants / ht.possible_variants,
                [1, ht.mu_snp],
                weight=ht.possible_variants if weighted else None,
            ).beta,
        )
    }
    return ht.aggregate(hl.struct(**agg_expr))


def build_coverage_model(coverage_ht: hl.Table) -> (float, float):
    """
    Calibrate coverage model.

    This function uses linear regression to build a model of log10(coverage) to this scaled ratio as a correction
    factor for low coverage sites.

    .. note::
        The following annotations should be present in `coverage_ht`:
            - low_coverage_obs_exp - an observed:expected ratio for a given coverage level
            - log_coverage - log10 coverage

    :param coverage_ht: Low coverage Table.
    :return: Tuple with intercept and slope of the model.
    """
    return tuple(
        coverage_ht.aggregate(
            hl.agg.linreg(
                coverage_ht.low_coverage_obs_exp, [1, coverage_ht.log_coverage]
            )
        ).beta
    )
