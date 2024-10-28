"""Script containing utility functions used in the constraint pipeline."""

import logging
from typing import Dict, List, Optional, Tuple

import hail as hl
import numpy as np
from gnomad.utils.constraint import (
    add_gencode_transcript_annotations,
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    compute_expected_variants,
    compute_pli,
    count_observed_and_possible_by_group,
    count_variants_by_group,
    get_annotation,
    get_constraint_flags,
    get_counts_agg_expr,
    get_pop_freq_indices,
    get_single_variant_count_expr,
    oe_aggregation_expr,
    oe_confidence_interval,
)
from gnomad.utils.filtering import (
    add_filters_expr,
    filter_by_numeric_expr_range,
    filter_for_mu,
    filter_to_autosomes,
)
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import filter_vep_transcript_csqs
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import (
    COVERAGE_CUTOFF,
    GENOMIC_REGIONS,
    get_checkpoint_path,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def prepare_ht_for_constraint_calculations(
    ht: hl.Table,
    exome_coverage_metric: str = "median",
) -> hl.Table:
    """
    Prepare Table for constraint calculations.

    Perform the following steps:

        - Add exome coverage annotation based on `exome_coverage_metric`.
        - Add genomic region annotation based on the locus. The genomic regions are:
            - 'autosome_or_par': autosomes and pseudoautosomal regions.
            - 'chrx_nonpar': Chromosome X non-pseudoautosomal regions.
            - 'chry_nonpar': Chromosome Y non-pseudoautosomal regions.

    :param ht: Annotated context Table.
    :param exome_coverage_metric: Metric to use for exome coverage. One of ["median",
        "AN", "AN_percent"]. Default is "median".
    :return: Table prepared for constraint calculations.
    """
    if exome_coverage_metric == "median":
        # Obtain field name for median exome coverage.
        exome_coverage_metric = (
            "median_approx" if "median_approx" in ht.coverage.exomes else "median"
        )
        cov_expr = ht.coverage.exomes[exome_coverage_metric]
    elif exome_coverage_metric == "AN":
        cov_expr = ht.AN.exomes[0]
    elif exome_coverage_metric == "AN_percent":
        # Calculate total allele number from strata_sample_count and annotate
        # exomes_AN_percent (percent samples with AN).
        cov_expr = hl.int(
            (ht.AN.exomes[0] / (ht.an_globals.exomes.strata_sample_count[0] * 2)) * 100
        )
    else:
        raise ValueError(
            f"Exome coverage metric must be one of ['median', 'AN', 'AN_percent'], not {exome_coverage_metric}"
        )
    logger.info("Setting 'exome_coverage' to %s", exome_coverage_metric)

    # Add annotation for exome coverage and genomic region (autosome, X non-par,
    # Y non-par).
    ht = ht.annotate(
        exome_coverage=cov_expr,
        genomic_region=(
            hl.case()
            .when(ht.locus.in_autosome_or_par(), "autosome_or_par")
            .when(ht.locus.in_x_nonpar(), "chrx_nonpar")
            .when(ht.locus.in_y_nonpar(), "chry_nonpar")
            .or_missing()
        ),
    )

    # TODO: Remove this when we have X and Y methylation levels.
    ht = ht.filter(ht.locus.in_autosome())

    return ht


# TODO: for pops, should possible be calculated for each population's AF?
# TODO: handle pops and downsampling
def create_observed_and_possible_ht(
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    max_af: float = 0.001,
    additional_grouping: Tuple = (),
    pops: Tuple = (),
    downsamplings: Optional[List[int]] = None,
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

    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling counts.
    :param grouping: Annotations other than 'context', 'ref', 'alt', and
        `methylation_level` to group by when counting variants. Default is
        ('exome_coverage',).
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :param filter_coverage_over_0: Whether to filter the exome Table and context Table
        to variants with `coverage_metric` larger than 0. Default is False.
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
    # Set up the criteria to keep sites with defined exome coverage, high-quality
    # variants and variants with exome coverage larger than 0 if requested. For variants
    # with defined exome coverage, but undefined high-quality filters, the variant is
    # kept so that it can be counted in the possible variant count.
    low_coverage_filter = low_coverage_filter or (1 if filter_coverage_over_0 else 0)
    context_ht = context_ht.filter(
        hl.is_defined(context_ht.exome_coverage)
        & (context_ht.exome_coverage >= low_coverage_filter)
        & hl.or_else(hl.len(context_ht.filters.exomes) == 0, True)
    )

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
        context_ht = filter_vep_transcript_csqs(
            context_ht, canonical=canonical, mane_select=mane_select
        )

    # Allele frequency information for high-quality genotypes (GQ >= 20; DP >= 10; and
    # AB >= 0.2 for heterozygous calls) in all release samples in gnomAD.
    freq_expr = context_ht.freq.exomes[0]

    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `grouping`, context, ref, alt, and methylation_level.
    keys = ("context", "ref", "alt", "methylation_level") + additional_grouping
    additional_grouping = keys + ("cpg", "mutation_type")
    # TODO: Change from a struct of structs...?
    ht = count_observed_and_possible_by_group(
        context_ht,
        freq_expr=freq_expr,
        additional_grouping=additional_grouping,
        partition_hint=partition_hint,
        use_table_group_by=True,
        max_af=max_af,
    ).key_by(*keys)

    # Annotate with mutation rate.
    ht = annotate_with_mu(ht, mutation_ht)

    # TODO: Remove repartition once partition_hint bugs are resolved.
    ht = ht.repartition(partition_hint)

    if global_annotation:
        ht = ht.annotate_globals(
            **{global_annotation: hl.struct(max_af=max_af, pops=pops)}
        )

    return ht


def compute_variant_predicted_probability_observed(
    plateau_models_expr: hl.StructExpression,
    mu_expr: hl.Float64Expression,
    cpg_expr: hl.BooleanExpression,
    cov_corr_expr: hl.Float64Expression = None,
):
    """
    Compute the predicted probability observed for a variant.

    The predicted probability observed is computed as the mutation rate adjusted by the
    plateau model for CpG transitions and non-CpG transitions. If a coverage correction
    factor is provided, the predicted probability observed is multiplied by the coverage
    correction factor.

    :param plateau_models_expr: Linear models that calibrate mutation rate to proportion
        observed for high coverage exome. It includes models for CpG sites, non-CpG
        sites, and each population in `POPS`.
    :param mu_expr: Mutation rate.
    :param cpg_expr: Boolean expression indicating whether the variant is a CpG
        transition.
    :param cov_corr_expr: Coverage correction factor. Default is None.
    :return: Predicted probability observed for a variant.
    """
    plateau_model = hl.literal(plateau_models_expr)[cpg_expr]
    slope = plateau_model[1]
    intercept = plateau_model[0]

    ppo_expr = mu_expr * slope + intercept

    if cov_corr_expr is not None:
        ppo_expr = ppo_expr * cov_corr_expr

    return ppo_expr


def get_expected_variants_agg_expr(
    ht: hl.Table,
    plateau_models_expr: hl.StructExpression,
    mu_expr: hl.Float64Expression,
    cov_corr_expr: hl.Float64Expression,
    possible_variants_expr: hl.Int64Expression,
    cpg_expr: hl.BooleanExpression,
):
    """
    Compute the expected number of variants.

    The expected number of variants is computed as the sum of the predicted probability
    observed multiplied by the possible variant counts. If a coverage correction factor
    is provided, the expected number of variants is multiplied by the coverage correction
    factor.

    :param ht: Table with observed and possible variant counts.
    :param plateau_models_expr: Linear models that calibrate mutation rate to proportion
        observed for high coverage exome. It includes models for CpG sites, non-CpG
        sites, and each population in `POPS`.
    :param mu_expr: Mutation rate.
    :param cov_corr_expr: Coverage correction factor.
    :param possible_variants_expr: Possible variant counts.
    :param cpg_expr: Boolean expression indicating whether the variant is a CpG
        transition.
    :return: Aggregation expression for expected number of variants.
    """

    def _get_ppo(cov_corr_expr: hl.Float64Expression) -> hl.Float64Expression:
        """
        Get the predicted probability observed for a variant.

        :param cov_corr_expr: Coverage correction factor.
        :return: Predicted probability observed for a variant.
        """
        return compute_variant_predicted_probability_observed(
            plateau_models_expr=plateau_models_expr,
            mu_expr=mu_expr,
            cpg_expr=cpg_expr,
            cov_corr_expr=cov_corr_expr,
        )

    ann_expr = {
        "observed_variants": ht.observed_variants,
        "possible_variants": ht.possible_variants,
        "predicted_proportion_observed": _get_ppo(None),
        "expected_variants": _get_ppo(cov_corr_expr) * possible_variants_expr,
    }

    return {k: hl.agg.sum(v) for k, v in ann_expr.items()}


def apply_models(
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Optional[Tuple[float, float]] = None,
    log10_coverage: bool = True,
    max_af: float = 0.001,
    additional_grouping: Tuple = ("methylation_level",),
    pops: Tuple = (),
    downsamplings: Optional[List[int]] = None,
    obs_pos_count_partition_hint: int = 2000,
    expected_variant_partition_hint: int = 1000,
    custom_vep_annotation: str = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    low_coverage_filter: int = None,
    use_mane_select: bool = True,
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
        predicted_proportion_observed = sum(plateau model slope * mu_agg + plateau model intercept) (separately for CpG transitions and other sites)
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
    :param log10_coverage: Whether to convert coverage sites with log10 when building the coverage model. Default is True.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling counts.
    :param obs_pos_count_partition_hint: Target number of partitions for
        aggregation when counting variants. Default is 2000.
    :param expected_variant_partition_hint: Target number of partitions for sum
        aggregators when computation is done. Default is 1000.
    :param custom_vep_annotation: The customized model (one of
        "transcript_consequences" or "worst_csq_by_gene"). Default is None.
    :param high_cov_definition: Median coverage cutoff. Sites with coverage above this cutoff
        are considered well covered and was used to build plateau models. Sites
        below this cutoff have low coverage and was used to build coverage models.
        Default is `COVERAGE_CUTOFF`.
    :param low_coverage_filter: Lower median coverage cutoff for coverage filter.
        Sites with coverage below this cutoff will be removed from`exome_ht` and
        'context_ht'.
    :param use_mane_select: Use MANE Select transcripts in grouping.
        Only used when `custom_vep_annotation` is set to 'transcript_consequences'.
        Default is True.
    :return: Table with `expected_variants` (expected variant counts) and `obs_exp`
        (observed:expected ratio) annotations.
    """
    # Filter context ht to sites with defined exome coverage_metric.
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    if low_coverage_filter is not None:
        context_ht = context_ht.filter(context_ht.exome_coverage >= low_coverage_filter)

    # Add necessary constraint annotations for grouping.
    if custom_vep_annotation == "worst_csq_by_gene":
        vep_annotation = "worst_csq_by_gene"
        if use_mane_select:
            raise ValueError(
                "'mane_select' cannot be set to True when custom_vep_annotation is set"
                " to 'worst_csq_by_gene'."
            )
    else:
        vep_annotation = "transcript_consequences"
        include_canonical_group = True
        include_mane_select_group = use_mane_select

    context_ht, grouping = annotate_exploded_vep_for_constraint_groupings(
        ht=context_ht,
        coverage_expr=context_ht.exome_coverage,
        vep_annotation=vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )
    grouping = additional_grouping + grouping

    # Compute observed and possible variant counts.
    ht = create_observed_and_possible_ht(
        context_ht=context_ht,
        mutation_ht=mutation_ht,
        max_af=max_af,
        additional_grouping=grouping,
        pops=pops,
        downsamplings=downsamplings,
        partition_hint=obs_pos_count_partition_hint,
        filter_coverage_over_0=True,
        transcript_for_synonymous_filter=None,
    )

    # NOTE: In v2 ht.mu_snp was incorrectly multiplied here by possible_variants, but
    #  this multiplication has now been moved,
    # so that it is applied after the regression within compute_expected_variants.
    mu_expr = ht.mu_snp
    poss_expr = ht.possible_variants
    # Determine coverage correction to use based on coverage value. If no
    # coverage model is provided, set to 1 as long as coverage > 0.
    if log10_coverage:
        logger.info("Converting coverage sites by log10.")
        cov_value = hl.log10(ht.coverage)
    else:
        cov_value = ht.coverage

    cov_corr_expr = (
        hl.case()
        .when(ht.coverage == 0, 0)
        .when(ht.coverage >= high_cov_definition, 1)
        .default(
            (coverage_model[1] * cov_value + coverage_model[0])
            if coverage_model is not None
            else 1
        )
    )

    # Generate sum aggregators for 'mu' on the entire dataset.
    agg_expr = {"mu": hl.agg.sum(mu_expr * cov_corr_expr)}
    agg_expr.update(
        compute_expected_variants(
            ht=ht,
            plateau_models_expr=plateau_models,
            mu_expr=mu_expr,
            cov_corr_expr=cov_corr_expr,
            possible_variants_expr=poss_expr,
            cpg_expr=ht.cpg,
        )
    )
    downsampling_meta = {}
    for pop in pops:
        agg_expr.update(
            compute_expected_variants(
                ht=ht,
                plateau_models_expr=plateau_models,
                mu_expr=mu_expr,
                cov_corr_expr=cov_corr_expr,
                possible_variants_expr=poss_expr,
                cpg_expr=ht.cpg,
                pop=pop,
            )
        )

        # Store which downsamplings are obtained for each pop in a
        # downsampling_meta dictionary.
        ds = hl.eval(get_pop_freq_indices(ht.freq_meta, pop=pop))
        key_names = {key for _, meta_dict in ds for key in meta_dict.keys()}
        genetic_ancestry_label = "gen_anc" if "gen_anc" in key_names else "pop"
        downsampling_meta[pop] = [
            x[1].get("downsampling", "all")
            for x in ds
            if x[1][genetic_ancestry_label] == pop
            and (
                int(x[1].get("downsampling", 0)) in downsamplings
                if downsamplings is not None
                else True
            )
        ]

    # Remove coverage from grouping.
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
            genetic_ancestry_groups=pops,
            plateau_models=plateau_models,
            coverage_model=coverage_model_global,
            high_cov_definition=high_cov_definition,
            log10_coverage=log10_coverage,
            downsampling_meta=downsampling_meta if downsampling_meta else "None",
        )
    )
    # Compute the observed:expected ratio.
    return ht.annotate(obs_exp=ht.observed_variants / ht.expected_variants)


def calculate_mu_by_downsampling(
    context_ht: hl.Table,
    count_singletons: bool = False,
    additional_grouping: Tuple[str] = ("methylation_level",),
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

    :param context_ht: Context Table for autosome/pseudoautosomal regions.
    :param count_singletons: Whether to count singletons. Default is False.
    :param additional_grouping: Annotations other than 'context', 'ref', and 'alt'.
        Default is ('methylation_level',).
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

    # Filter to autosomal sites (remove pseudoautosomal regions).
    context_ht = filter_to_autosomes(context_ht)

    # Filter to sites with mean genome coverage between min_cov and max_cov.
    context_ht = filter_by_numeric_expr_range(
        context_ht, context_ht.coverage.genomes.mean, (min_cov, max_cov)
    )

    # Filter the Table so that the most severe annotation is 'intron_variant' or
    # 'intergenic_variant', and that the GERP score is between 'gerp_lower_cutoff' and
    # 'gerp_upper_cutoff' (ideally these values will define the 5th and 95th
    # percentile of the genome-wide distribution).
    context_ht = filter_for_mu(context_ht, gerp_lower_cutoff, gerp_upper_cutoff)

    # Get the frequency index of downsampling with size of `downsampling_level`.
    freq_meta = hl.eval(context_ht.freq_globals.genomes.freq_meta)
    downsampling_idx = freq_meta.index(
        {"group": "adj", "pop": "global", "downsampling": str(downsampling_level)}
    )
    freq_expr = context_ht.freq.genomes[downsampling_idx]

    # Set up the criteria to keep high-quality sites, and sites found in less than or
    # equal to 'ac_cutoff' copies in the downsampled set.
    # Count possible variants in context Table, only keeping variants not in the genome
    # dataset, or with AC <= 'ac_cutoff' and passing filters.
    context_ht = context_ht.filter(
        hl.or_else(
            (hl.len(context_ht.filters.genomes) == 0) & (freq_expr.AC <= ac_cutoff),
            True,
        )
    )

    meta_keep = [{"group": "adj"}]
    meta_keep += [{"group": "adj", "pop": pop} for pop in pops if pop != "global"]
    meta_keep += [
        {"group": "adj", "pop": pop, "downsampling": str(downsampling_level)}
        for pop in pops
    ]
    freq_expr = hl.array(
        [context_ht.freq.genomes[freq_meta.index(m)] for m in meta_keep]
    )

    # Count the observed variants in the entire Table and in each downsampling grouped
    # by context, ref, alt, and 'additional_grouping'.
    ht = count_observed_and_possible_by_group(
        context_ht,
        freq_expr=freq_expr,
        count_singletons=count_singletons,
        additional_grouping=additional_grouping,
        use_table_group_by=True,
    )

    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants[0])) // 3
    logger.info(
        "Total bases to use when calculating correction_factors: %f", total_bases
    )

    # Compute the proportion observed, which represents the relative mutability of each
    # variant class.
    downsampling_idx = meta_keep.index(
        {"group": "adj", "pop": "global", "downsampling": str(downsampling_level)}
    )
    po_expr = ht.observed_variants / ht.possible_variants[0]
    correction_factors = ht.aggregate(
        total_mu / (hl.agg.array_sum(ht.observed_variants) / total_bases),
        _localize=False,
    )
    mu_expr = correction_factors * ht.observed_variants / ht.possible_variants[0]
    ht = ht.annotate(
        proportion_observed=po_expr,
        mu=mu_expr,
        mu_snp=mu_expr[downsampling_idx],
    )
    ht = ht.annotate_globals(
        ac_cutoff=ac_cutoff,
        total_mu=total_mu,
        mu_meta=meta_keep,
        downsampling_level=downsampling_level,
        downsampling_idx=downsampling_idx,
        min_cov=min_cov,
        max_cov=max_cov,
        gerp_lower_cutoff=gerp_lower_cutoff,
        gerp_upper_cutoff=gerp_upper_cutoff,
    )

    return annotate_mutation_type(ht)


def add_oe_lof_upper_rank_and_bin(
    ht: hl.Table, use_mane_select_over_canonical: bool = True
) -> hl.Table:
    """
    Compute the rank and decile of the lof oe upper confidence interval for MANE Select or canonical ensembl transcripts.

    :param ht: Input Table with the value for the lof oe upper confidence interval stored in ht.lof.oe_ci.upper.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical transcripts for filtering the Table.
        If a gene does not have a MANE Select transcript, the canonical transcript (if available) will be used instead. Default is True.
    :return: Table with anntotations added for 'upper_rank', 'upper_bin_decile'.
    """
    # Filter to only ensembl transcripts of the specified transcript filter. If MANE
    # select is specified, and a gene does not have a MANE select transcript, use
    # canonical instead.
    if use_mane_select_over_canonical:
        genes = ht.group_by(ht.gene_id).aggregate(
            mane_present=hl.agg.any(ht.mane_select),
            canonical_present=hl.agg.any(ht.canonical),
        )

        genes = genes.annotate(
            only_canonical=~(genes.mane_present) & (genes.canonical_present)
        )

        ms_ht = ht.annotate(
            _only_canonical=genes[ht.gene_id].only_canonical,
            _mane_present=genes[ht.gene_id].mane_present,
        )
        total_count = ms_ht.count()
        ms_ht = ms_ht.filter(
            (ms_ht.transcript.startswith("ENST"))
            & (
                (ms_ht._mane_present & ms_ht.mane_select)
                | (ms_ht._only_canonical & ms_ht.canonical)
            )
        )
        filtered_count = ms_ht.count()
        logger.info(
            "Retaining %d out of %d transcripts to use for rank annotations.",
            filtered_count,
            total_count,
        )
    else:
        ms_ht = ht.filter((ht.canonical) & (ht.transcript.startswith("ENST")))

    # Rank lof.oe_ci.upper in ascending order.
    ms_ht = ms_ht.order_by(ms_ht.lof.oe_ci.upper).add_index(name="upper_rank")

    # Determine decile bins.
    n_transcripts = ms_ht.count()
    ms_ht = ms_ht.annotate(
        upper_bin_decile=hl.int(ms_ht.upper_rank * 10 / n_transcripts)
    )

    # Add rank and bin annotations back to original Table.
    ms_ht = ms_ht.key_by(*list(ht.key))
    ms_index = ms_ht[ht.key]
    ht = ht.annotate(
        lof=ht.lof.annotate(
            oe_ci=ht.lof.oe_ci.annotate(
                upper_rank=ms_index.upper_rank,
                upper_bin_decile=ms_index.upper_bin_decile,
            )
        )
    )

    return ht


def compute_constraint_metrics(
    ht: hl.Table,
    gencode_ht: hl.Table,
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Tuple[str] = (
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ),
    pops: Tuple[str] = (),
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
    raw_z_outlier_threshold_lower_lof: float = -8.0,
    raw_z_outlier_threshold_lower_missense: float = -8.0,
    raw_z_outlier_threshold_lower_syn: float = -8.0,
    raw_z_outlier_threshold_upper_syn: float = 8.0,
    include_os: bool = False,
    use_mane_select_over_canonical: bool = True,
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
    :param raw_z_outlier_threshold_lower_lof: Value at which the raw z-score is considered an outlier for lof variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_lower_missense: Value at which the raw z-score is considered an outlier for missense variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_lower_syn: Lower value at which the raw z-score is considered an outlier for synonymous variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_upper_syn: Upper value at which the raw z-score is considered an outlier for synonymous variants. Values above this threshold will be considered outliers. Default is  8.0.
    :param include_os: Whether or not to include OS (other splice) as a grouping when
        stratifying calculations by lof HC.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical transcripts for filtering the Table when determining ranks for the lof oe upper confidence interval.
        If a gene does not have a MANE Select transcript, the canonical transcript (if available) will be used instead. Default is True.
    :param gencode_ht: Table containing GENCODE annotations.
    :return: Table with pLI scores, observed:expected ratio, confidence interval of the
        observed:expected ratio, and z scores.
    """
    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X.
    # `annotation_dict` stats the rule of filtration for each annotation.
    annotation_dict = {
        # Filter to classic LoF annotations with LOFTEE HC or LC.
        "lof_hc_lc": hl.literal(set(classic_lof_annotations)).contains(ht.annotation)
        & ((ht.modifier == "HC") | (ht.modifier == "LC")),
        # Filter to LoF annotations with LOFTEE HC.
        "lof": ht.modifier == "HC",
        # Filter to LoF annotations with LOFTEE HC and flagged alpha missense variants.
        "lof_and_alphamissense": (ht.modifier == "HC") | ht.am,
        # Filter to missense variants.
        "mis": ht.annotation == "missense_variant",
        # Filter to probably damaging missense variants predicted by PolyPen-2.
        "mis_pphen": ht.modifier == "probably_damaging",
        # Filter to flagged alpha missense variants.
        "mis_alphamissense": ht.am,
        # Filter to synonymous variants.
        "syn": ht.annotation == "synonymous_variant",
    }

    # Define two lists of 'annotation_dict' keys that require different computations.
    # The 90% CI around obs:exp and z-scores are only computed for lof, mis, and syn.
    oe_ann = ["lof", "lof_and_alphamissense", "mis_alphamissense", "mis", "syn"]
    # pLI scores are only computed for LoF variants.
    lof_ann = ["lof_hc_lc", "lof"]

    # Create dictionary with outlier z-score thresholds with annotation as key
    # and list of thresholds [lower, upper] as values.
    z_score_outlier_dict = {
        "lof": [raw_z_outlier_threshold_lower_lof, None],
        "lof_and_alphamissense": [raw_z_outlier_threshold_lower_lof, None],
        "mis_alphamissense": [raw_z_outlier_threshold_lower_lof, None],
        "mis": [raw_z_outlier_threshold_lower_missense, None],
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

        gen_anc_lower_struct = {}
        gen_anc_upper_struct = {}
        gen_anc_z_raw_struct = {}

        # Calculate lower and upper cis, and raw z scores for each pop, excluding
        # downsamplings.
        for pop in pops:
            obs_expr = ht[ann]["gen_anc_obs"][pop][0]
            exp_expr = ht[ann]["gen_anc_exp"][pop][0]
            oe_ci_expr = oe_confidence_interval(obs_expr, exp_expr)
            raw_z_expr = calculate_raw_z_score(obs_expr, exp_expr)

            gen_anc_lower_struct[pop] = oe_ci_expr.lower
            gen_anc_upper_struct[pop] = oe_ci_expr.upper
            gen_anc_z_raw_struct[pop] = raw_z_expr

    # Annotate the table with the structs.
    ann_expr[ann] = ann_expr[ann].annotate(
        gen_anc_oe_ci=hl.struct(
            lower=hl.struct(**gen_anc_lower_struct),
            upper=hl.struct(**gen_anc_upper_struct),
        ),
        gen_anc_z_raw=hl.struct(**gen_anc_z_raw_struct),
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

    ht = ht.checkpoint(new_temp_file(prefix="z_scores", extension="ht"))

    # Compute the rank and decile of the lof oe upper confidence
    # interval for MANE Select or canonical ensembl transcripts.
    ht = add_oe_lof_upper_rank_and_bin(
        ht, use_mane_select_over_canonical=use_mane_select_over_canonical
    )

    # Add transcript annotations from GENCODE.
    ht = add_gencode_transcript_annotations(ht, gencode_ht)

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
