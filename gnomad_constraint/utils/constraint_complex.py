"""Script containing utility functions used in the constraint pipeline."""

import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
import numpy as np
from gnomad.assessment.summary_stats import generate_filter_combinations
from gnomad.utils.constraint import (
    add_gencode_transcript_annotations,
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    calculate_raw_z_score,
    get_annotation,
    get_counts_agg_expr,
    get_pop_freq_indices,
    get_single_variant_count_expr,
    obs_exp_aggregation_expr,
    oe_aggregation_expr,
    oe_confidence_interval,
)
from gnomad.utils.vep import filter_vep_transcript_csqs
from hail.utils.misc import divide_null, new_temp_file

from gnomad_constraint.resources.resource_utils import COVERAGE_CUTOFF

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def get_coverage_expr(ht, exome_coverage_metric):
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

    return cov_expr


def get_annotations_for_computing_mu(
    locus_expr: hl.expr.LocusExpression,
    genomes_filter_expr: hl.expr.StructExpression,
    genomes_freq_expr: hl.expr.StructExpression,
    genomes_freq_meta: List[Dict[str, str]],
    genomes_coverage_expr: hl.expr.Int32Expression,
    gerp_expr: hl.expr.Float64Expression,
    most_severe_consequence_expr: hl.expr.StringExpression,
    pops: Optional[List[str]] = None,
    downsampling_level: int = 1000,
    min_cov: int = 15,
    max_cov: int = 60,
    gerp_lower_cutoff: float = -3.9885,
    gerp_upper_cutoff: float = 2.6607,
    ac_cutoff: int = 5,
):
    """
    Filter to non-coding annotations and remove GERP outliers.

    .. note::

        Values for `gerp_lower_cutoff` and `gerp_upper_cutoff` default to -3.9885 and
        2.6607, respectively. These values were precalculated on the GRCh37 context
        table and define the 5th and 95th percentiles.

    :param ht: Input Table.
    :param gerp_lower_cutoff: Minimum GERP score for variant to be included. Default is -3.9885.
    :param gerp_upper_cutoff: Maximum GERP score for variant to be included. Default is 2.6607.
    :return: Table filtered to intron or intergenic variants with GERP outliers removed.
    """
    genomes_freq_expr, genomes_freq_meta = filter_freq_for_constraint(
        genomes_freq_expr,
        genomes_freq_meta,
        pops=None,
        downsamplings=[downsampling_level],
        downsampling_pops=pops,
        pop_label="pop",
    )
    downsampling_idx = genomes_freq_meta.index(
        {"group": "adj", "pop": "global", "downsampling": str(downsampling_level)}
    )

    # Filter to autosomal sites (remove pseudoautosomal regions).
    keep_expr = locus_expr.in_autosome()

    # Filter to sites with mean genome coverage between min_cov and max_cov.
    keep_expr &= (genomes_coverage_expr >= min_cov) & (genomes_coverage_expr <= max_cov)

    # Filter to sites where the GERP score is between 'gerp_lower_cutoff' and
    # 'gerp_upper_cutoff' (ideally these values will define the 5th and 95th
    # percentile of the genome-wide distribution).
    keep_expr &= (gerp_expr > gerp_lower_cutoff) & (gerp_expr < gerp_upper_cutoff)

    # Filter so that the most severe annotation is 'intron_variant' or
    # 'intergenic_variant'
    keep_expr &= (most_severe_consequence_expr == "intron_variant") | (
        most_severe_consequence_expr == "intergenic_variant"
    )

    # Set up the criteria to keep high-quality sites, and sites found in less than or
    # equal to 'ac_cutoff' copies in the downsampled set.
    # Count possible variants in context Table, only keeping variants not in the genome
    # dataset, or with AC <= 'ac_cutoff' and passing filters.
    genomes_filter_freq_expr = genomes_freq_expr[downsampling_idx]
    keep_expr &= hl.or_else(
        (hl.len(genomes_filter_expr) == 0) & (genomes_filter_freq_expr.AC <= ac_cutoff),
        True,
    )
    obs_pos_expr = hl.struct(
        observed=genomes_freq_expr.map(
            lambda x: get_single_variant_count_expr(freq_expr=x)
        ),
        possible=get_single_variant_count_expr(
            freq_expr=genomes_freq_expr[0], count_missing=True
        ),
    )
    obs_pos_expr = hl.struct(
        freq=genomes_freq_expr,
        used_in_mu=keep_expr,
        **hl.or_missing(keep_expr, obs_pos_expr),
    )
    obs_pos_globals = hl.struct(
        freq_meta=genomes_freq_meta,
        ac_cutoff=ac_cutoff,
        min_cov=min_cov,
        max_cov=max_cov,
        gerp_lower_cutoff=gerp_lower_cutoff,
        gerp_upper_cutoff=gerp_upper_cutoff,
        downsampling_level=downsampling_level,
        downsampling_idx=downsampling_idx,
        most_severe_consequence=["intron_variant", "intergenic_variant"],
    )

    return obs_pos_expr, obs_pos_globals


def get_exomes_observed_and_possible(
    exomes_filter_expr: hl.expr.SetExpression,
    exomes_freq_expr: hl.expr.ArrayExpression,
    exomes_freq_meta: List[Dict[str, str]],
    exomes_coverage_expr: hl.expr.Int32Expression,
    pops: Optional[List[str]] = None,
    downsamplings: Optional[List[int]] = None,
    max_af: float = 0.001,
):
    # Filter frequency array for computing the observed expression on all requested
    # populations and downsamplings.
    exomes_freq_expr, exomes_freq_meta = filter_freq_for_constraint(
        exomes_freq_expr,
        exomes_freq_meta,
        pops=pops,
        downsamplings=downsamplings,
        downsampling_pops=pops if downsamplings is not None else None,
    )
    keep_expr = hl.is_defined(exomes_coverage_expr) & hl.or_else(
        hl.len(exomes_filter_expr) == 0, True
    )

    obs_pos_expr = hl.struct(
        observed=exomes_freq_expr.map(
            lambda x: get_single_variant_count_expr(freq_expr=x, max_af=max_af)
        ),
        possible=get_single_variant_count_expr(
            freq_expr=exomes_freq_expr[0], max_af=max_af, count_missing=True
        ),
    )
    obs_pos_expr = hl.struct(
        exomes_freq=exomes_freq_expr,
        used_in_apply=keep_expr,
        **hl.or_missing(keep_expr, obs_pos_expr),
    )
    obs_pos_globals = hl.struct(
        exomes_freq_meta=exomes_freq_meta,
        genetic_ancestry_groups=pops,
        downsamplings=downsamplings or "None",
    )

    return obs_pos_expr, obs_pos_globals


def prepare_ht_for_constraint_calculations(
    ht: hl.Table,
    exome_coverage_metric: str = "median",
    pops: Optional[List[str]] = None,
    downsamplings: Optional[List[int]] = None,
    mu_downsampling_level: int = 1000,
    min_cov: int = 15,
    max_cov: int = 60,
    gerp_lower_cutoff: float = -3.9885,
    gerp_upper_cutoff: float = 2.6607,
    ac_cutoff: int = 5,
    max_af: float = 0.001,
    build_model_low_cov_cutoff: int = None,
    build_model_high_cov_cutoff: int = COVERAGE_CUTOFF,
    build_model_upper_cov_cutoff: int = None,
    apply_model_low_cov_cutoff: int = COVERAGE_CUTOFF,
    apply_model_high_cov_cutoff: int = COVERAGE_CUTOFF,
    skip_coverage_model: bool = False,
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
    # Add annotation for exome coverage and genomic region (autosome, X non-par,
    # Y non-par).
    genomic_region_expr = (
        hl.case()
        .when(ht.locus.in_autosome_or_par(), "autosome_or_par")
        .when(ht.locus.in_x_nonpar(), "chrx_nonpar")
        .when(ht.locus.in_y_nonpar(), "chry_nonpar")
        .or_missing()
    )

    mu_expr, mu_globals = get_annotations_for_computing_mu(
        ht.locus,
        ht.filters.genomes,
        ht.freq.genomes,
        ht.freq_globals.genomes.freq_meta,
        ht.coverage.genomes.mean,
        ht.gerp,
        ht.vep.most_severe_consequence,
        pops=pops,
        downsampling_level=mu_downsampling_level,
        min_cov=min_cov,
        max_cov=max_cov,
        gerp_lower_cutoff=gerp_lower_cutoff,
        gerp_upper_cutoff=gerp_upper_cutoff,
        ac_cutoff=ac_cutoff,
    )

    exomes_coverage_expr = get_coverage_expr(ht, exome_coverage_metric)
    exomes_obs_pos_expr, exomes_obs_pos_globals = get_exomes_observed_and_possible(
        ht.filters.exomes,
        ht.freq.exomes,
        ht.freq_globals.exomes.freq_meta,
        exomes_coverage_expr,
        pops=pops,
        downsamplings=downsamplings,
        max_af=max_af,
    )

    ann_expr = {}
    if build_model_high_cov_cutoff is not None:
        build_model_high_cov_expr = exomes_coverage_expr >= build_model_high_cov_cutoff

        # Filter to sites with coverage_expr equal to or below `upper_cov_cutoff` if
        # specified.
        if build_model_upper_cov_cutoff is not None:
            build_model_high_cov_expr &= (
                exomes_coverage_expr <= build_model_upper_cov_cutoff
            )

        if build_model_low_cov_cutoff is not None:
            build_model_low_cov_expr = (
                exomes_coverage_expr >= build_model_low_cov_cutoff
            )
        else:
            build_model_low_cov_expr = hl.literal(False)

        if skip_coverage_model:
            add_low_coverage_model = hl.literal(False)
        else:
            add_low_coverage_model = exomes_coverage_expr > 0

        ann_expr["build_model"] = (
            hl.case()
            .when(build_model_high_cov_expr, ("high_cov", hl.missing(hl.tint)))
            .when(
                build_model_low_cov_expr & add_low_coverage_model,
                ("low_cov", exomes_coverage_expr),
            )
            .or_missing()
        )

    if apply_model_high_cov_cutoff is not None:
        ann_expr["apply_model"] = (
            hl.case()
            .when(exomes_coverage_expr >= apply_model_high_cov_cutoff, "high_cov")
            .when(exomes_coverage_expr > apply_model_low_cov_cutoff, "low_cov")
            .or_missing()
        )

    ht = ht.annotate(
        genomic_region=genomic_region_expr,
        mu=mu_expr,
        exomes_coverage=exomes_coverage_expr,
        **exomes_obs_pos_expr,
        **hl.or_missing(exomes_obs_pos_expr.used_in_apply, hl.struct(**ann_expr)),
    )
    ht = ht.select_globals(
        calculate_mu_globals=mu_globals,
        build_models_globals=hl.struct(
            plateau_models=hl.struct(),
            coverage_model=hl.struct(),
            build_model_high_cov_cutoff=build_model_high_cov_cutoff,
            build_model_upper_cov_cutoff=build_model_upper_cov_cutoff,
        ),
        apply_models_globals=hl.struct(
            apply_model_high_cov_cutoff=apply_model_high_cov_cutoff,
        ),
        **exomes_obs_pos_globals,
    )

    # TODO: Remove this when we have X and Y methylation levels.
    ht = ht.filter(ht.locus.in_autosome())

    return ht


def filter_context_for_observed_and_possible_counts(
    context_ht: hl.Table,
    synonymous_transcript_filter_field: str = None,
) -> hl.Table:
    """
    Filter the fully annotated context Table for observed and possible variant counts.

    The following variants are kept:

        - High-quality exome variants based on `context_ht.filters.exomes`.
        - Variants with defined exome coverage based on `context_ht.exome_coverage`.

    And if requested:

        - High exome coverage sites: `context_ht.exome_coverage` greeater than or equal
            to `low_coverage_filter`.
        - Variants that are synonymous in either MANE Select or canonical transcripts,
            if specified (`synonymous_transcript_filter_field`).

    :param context_ht: Preprocessed context Table.
    :param low_coverage_filter: Lower median coverage cutoff for coverage filter. Sites
        with coverage below this cutoff will be removed from the `exome_ht` and
        'context_ht'.
    :param synonymous_transcript_filter_field: Transcript to use when filtering to
        synonymous variants. Choices: ["mane_select", "canonical", None]. If
        "canonical", will filter to variants with a synonymous consequence in Ensembl
        canonical transcripts. If "mane_select", will filter to variants with a
        synonymous consequence in MANE Select transcripts. If None, no
        transcript/synonymous filter will be applied. Default is None.
    :return: Filtered context Table.
    """
    # If requested keep only variants that are synonymous in either MANE Select or
    # canonical transcripts.
    if synonymous_transcript_filter_field is not None:
        if synonymous_transcript_filter_field == "canonical":
            canonical, mane_select = True, False
        elif synonymous_transcript_filter_field == "mane_select":
            canonical, mane_select = False, True
        else:
            raise ValueError(
                "If synonymous_transcript_filter_field is not None, must be either"
                " 'canonical' or 'mane_select'"
            )
        context_ht = filter_vep_transcript_csqs(
            context_ht, canonical=canonical, mane_select=mane_select
        )

    return context_ht


def filter_freq_for_constraint(
    freq_expr,
    freq_meta_expr,
    pops: Optional[List[str]] = None,
    downsamplings: Optional[List[int]] = None,
    downsampling_pops: Optional[List[str]] = None,
    pop_label="gen_anc",
) -> Tuple[hl.ArrayExpression, List[Dict[str, str]]]:
    """
    Filter the frequency array for constraint calculations.

    The frequency array is filtered to include only adj frequencies for
    the populations in `pops` and the downsamplings in `downsamplings`, for the
    populations in `downsampling_pops`.

    No matter the input, the frequency array is always filtered to include the
    "adj" frequency for the full dataset.

    If `downsamplings` is None, no downsamplings are included. If `downsamplings` is
    provided, and `downsampling_pops` is None, only the "global" downsampling is
    included. If `downsampling_pops` is provided, the downsamplings for the populations
    in `downsampling_pops` are included as well as the "global" downsampling.

    :param freq_expr: Frequency array.
    :param freq_meta_expr: Frequency metadata array.
    :param pops: Optional list of populations to include in the frequency array. Default
        is None.
    :param downsamplings: Optional list of downsampling indices to include in the
        frequency array. Default is None.
    :param downsampling_pops: Optional list of populations to include in the frequency
        array for downsamplings. Default is None.
    :return: Filtered frequency array and metadata.
    """
    freq_meta = hl.eval(freq_meta_expr)
    meta_keep = [{"group": "adj"}]

    if pops is not None:
        meta_keep += [{"group": "adj", pop_label: pop} for pop in pops]

    if downsamplings is not None:
        downsampling_pops = ["global"] + (downsampling_pops or [])
        meta_keep += [
            {"group": "adj", pop_label: pop, "downsampling": str(ds)}
            for pop in downsampling_pops
            for ds in downsamplings
        ]

    meta_keep = [m for m in meta_keep if m in freq_meta]
    freq_expr = hl.array([freq_expr[freq_meta.index(m)] for m in meta_keep])

    return freq_expr, meta_keep


def get_weighted_agg_sum_expr(
    expr: Union[hl.ArrayExpression, hl.NumericExpression],
    weight_expr: Union[hl.ArrayExpression, hl.NumericExpression],
) -> Union[hl.ArrayExpression, hl.NumericExpression]:
    """
    Get the expression for a weighted aggregate sum.

    The function will return the sum of the product of `expr` and `weight_expr`.

    The parameters `expr` and `weight_expr` can be either ArrayExpression or
    NumericExpression. If both are ArrayExpression, each element of the arrays will be
    multiplied together and summed. If one is ArrayNumericExpression and the other is
    NumericExpression, the NumericExpression will be multiplied by each element
    of the ArrayExpression and summed.

    :param expr: Expression to be weighted and summed.
    :param weight_expr: Expression to weight `expr` by.
    :return: Weighted aggregate sum expression.
    """
    expr_is_array = isinstance(expr, hl.expr.ArrayNumericExpression)
    weight_is_array = isinstance(weight_expr, hl.expr.ArrayNumericExpression)
    if expr_is_array and weight_is_array:
        return hl.agg.array_sum(hl.zip(expr, weight_expr).map(lambda x: x[0] * x[1]))
    elif not expr_is_array and not weight_is_array:
        return hl.agg.sum(expr * weight_expr)
    else:
        return hl.agg.array_sum(expr * weight_expr)


def count_observed_and_possible_by_group(
    ht: hl.Table,
    possible_expr: hl.expr.Int32Expression,
    observed_expr: hl.expr.ArrayExpression,
    additional_grouping: Union[List[str]] = ("methylation_level",),
    partition_hint: int = 100,
    weight_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
    additional_agg_sum_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
) -> Union[hl.Table, Any]:
    if isinstance(weight_exprs, list):
        weight_exprs = {k: ht[k] for k in weight_exprs}
    if isinstance(additional_agg_sum_exprs, list):
        additional_agg_sum_exprs = {k: ht[k] for k in additional_agg_sum_exprs}

    weight_exprs = weight_exprs or {}
    additional_agg_sum_exprs = additional_agg_sum_exprs or {}

    # Build the grouping struct for the variant count aggregation.
    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    grouping = grouping.annotate(
        **{g: ht[g] for g in additional_grouping if g not in grouping}
    )
    logger.info(
        "The following annotations will be used to group the input Table rows when"
        " counting variants: %s.",
        ", ".join(grouping.keys()),
    )

    agg_expr = {
        "observed_variants": hl.agg.array_sum(observed_expr),
        "possible_variants": hl.agg.sum(possible_expr),
    }

    # Update the possible variant count aggregation expression to include weighted sums
    # of possible variant counts.
    agg_expr.update(
        {
            k: get_weighted_agg_sum_expr(possible_expr, v)
            for k, v in weight_exprs.items()
        }
    )

    # Get sum aggregation expressions for requested fields.
    agg_expr.update(
        {
            k: (
                hl.agg.array_sum(v)
                if isinstance(v, hl.ArrayExpression)
                else hl.agg.sum(v)
            )
            for k, v in additional_agg_sum_exprs.items()
        }
    )

    # Apply each variant count aggregation in `agg_expr` to get counts for all
    # combinations of `grouping`.
    ht = ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**agg_expr)

    return ht


def create_observed_and_possible_ht(
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    additional_grouping: Tuple = (),
    partition_hint: int = 100,
    keys=("context", "ref", "alt", "methylation_level"),
    weight_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
    additional_agg_sum_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
) -> hl.Table:
    """
    Count the observed variants and possible variants by substitution, context, methylation level, and additional `grouping`.

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation_level
        - observed_variants - observed variant counts
        - possible_variants - possible variant counts
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling counts.
    :param additional_grouping: Annotations other than 'context', 'ref', 'alt', and
        `methylation_level` to group by when counting variants. Default is
        ('exome_coverage',).
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :param global_annotation: The annotation name to use as a global StructExpression
        annotation containing input parameter values. If no value is supplied, this
        global annotation will not be added. Default is None.
    :return: Table with observed variant and possible variant count.
    """
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `grouping`, context, ref, alt, and methylation_level.
    # keys = keys + additional_grouping
    # additional_grouping = keys + ("cpg", "mutation_type")
    ht = count_observed_and_possible_by_group(
        context_ht,
        context_ht.possible,
        context_ht.observed,
        additional_grouping=additional_grouping,
        partition_hint=partition_hint,
        weight_exprs=weight_exprs,
        additional_agg_sum_exprs=additional_agg_sum_exprs,
    )
    ht = ht.key_by(*keys)

    # TODO: Remove repartition once partition_hint bugs are resolved.
    ht = ht.repartition(partition_hint)

    return ht


def create_training_set(
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    low_coverage_filter=30,
    synonymous_transcript_filter_field="mane_select",
    partition_hint=100,
) -> hl.Table:
    context_ht = filter_context_for_observed_and_possible_counts(
        context_ht,
        synonymous_transcript_filter_field=synonymous_transcript_filter_field,
    )

    ht = create_observed_and_possible_ht(
        context_ht,
        mutation_ht,
        additional_grouping=("genomic_region", "build_model"),
        partition_hint=partition_hint,
    )

    # Annotate with mutation rate.
    ht = annotate_with_mu(ht, mutation_ht)

    ht = ht.annotate_globals(
        training_dataset_params=hl.struct(
            low_coverage_filter=low_coverage_filter,
            synonymous_transcript_filter_field=synonymous_transcript_filter_field,
        )
    )

    return ht


def run_build_models(
    ht: hl.Table,
    skip_coverage_model: bool = False,
    log10_coverage: bool = True,
    additional_grouping=(),
) -> Tuple[Optional[Tuple[float, float]], hl.expr.StructExpression]:
    """
    Build coverage and plateau models.

    This function builds models (plateau_models) using linear regression to calibrate
    mutation rate estimates against the proportion observed of each substitution,
    context, and methylation level in `coverage_ht`.

    Two plateau models are fit, one for CpG transitions, and one for the remainder of
    sites (transversions and non CpG transitions).

    The plateau models only consider high coverage sites, or sites above a median
    coverage of `high_cov_definition` and median coverage below `upper_cov_cutoff`.

    Plateau model: adjusts proportion of expected variation based on location in the
    genome and CpG status.
    The x and y of the plateau models:
    - x: `mu_snp` - mutation rate
    - y: proportion observed ('observed_variants' or 'observed_{pop}' / 'possible_variants')

    This function also builds models (coverage models) to calibrate the proportion of
    expected variation at low coverage sites (sites below `high_cov_definition`).

    The coverage models are built by creating a scaling factor across all high coverage
    sites, applying this ratio to the low coverage sites, and running a linear
    regression.

    Coverage model: corrects proportion of expected variation at low coverage sites.
    Low coverage sites are defined as sites with median coverage < `high_cov_definition`.

    The x and y of the coverage model:
    - x: log10 groupings of exome coverage at low coverage sites
    - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage sites

    `high_coverage_scale_factor` = sum('observed_variants') /
                        sum('possible_variants' * 'mu_snp') at high coverage sites

    .. note::

        This function expects that the input Table(`coverage_ht`) was created using
        `get_proportion_observed_by_coverage`, which means that `coverage_ht` should
        contain only high quality synonymous variants below 0.1% frequency.

        This function also expects that the following fields are present in
        `coverage_ht`:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate allele
        - methylation_level - methylation level
        - cpg - whether the site is CpG site
        - observed_variants - the number of observed variants in the dataset for each
        variant. Note that the term "variant" here refers to a specific substitution,
        context, methylation level, and coverage combination
        - downsampling_counts_{pop} (optional) - array of observed variant counts per
        population after downsampling. Used only when `pops` is specified.
        - mu_snp - mutation rate
        - possible_variants - the number of possible variants in the dataset for each
        variant

    :param coverage_ht: Input coverage Table.
    :param coverage_expr: Expression that defines the coverage metric.
    :param weighted: Whether to weight the plateau models (a linear regression
        model) by 'possible_variants'. Default is False.
    :param pops: List of populations used to build plateau models.
        Default is ().
    :param keys: Annotations used to group observed and possible variant counts.
        Default is ("context", "ref", "alt", "methylation_level", "mu_snp").
    :param high_cov_definition: Lower median coverage cutoff. Sites with coverage above this cutoff
        are considered well covered. Default is `COVERAGE_CUTOFF`.
    :param upper_cov_cutoff: Upper median coverage cutoff. Sites with coverage above this cutoff
        are excluded from the high coverage Table. Default is None.
    :param skip_coverage_model: Whether to skip generating the coverage model. If set to True,
        None is returned instead of the coverage model. Default is False.
    :param log10_coverage: Whether to convert coverage sites with log10 when building
        the coverage model. Default is True.
    :return: Coverage model and plateau models.
    """
    high_cov_ht = ht.filter(ht.build_model[0] == "high_cov")

    # Build plateau models.
    plateau_models_agg_expr = build_plateau_models(
        cpg_expr=high_cov_ht.cpg,
        mu_snp_expr=high_cov_ht.mu_snp,
        observed_variants_expr=high_cov_ht.observed_variants,
        possible_variants_expr=high_cov_ht.possible_variants,
        additional_grouping={g: high_cov_ht[g] for g in additional_grouping},
    )
    plateau_models = high_cov_ht.aggregate(plateau_models_agg_expr)

    if not skip_coverage_model:
        low_cov_ht = ht.filter(ht.build_model[0] == "low_cov")

        # Create a metric that represents the relative mutability of the exome calculated
        # on high coverage sites and will be used as scaling factor when building the
        # coverage model.
        high_coverage_scale_factor = high_cov_ht.aggregate(
            hl.agg.sum(high_cov_ht.observed_variants[0])
            / hl.agg.sum(high_cov_ht.possible_variants * high_cov_ht.mu_snp)
        )

        low_cov_group_ht = low_cov_ht.group_by(
            cov_value=low_cov_ht.build_model[1]
        ).aggregate(
            low_coverage_oe=hl.agg.sum(low_cov_ht.observed_variants[0])
            / (
                high_coverage_scale_factor
                * hl.agg.sum(low_cov_ht.possible_variants * low_cov_ht.mu_snp)
            )
        )

        # Generate a Table with all necessary annotations (x and y listed above)
        # for the coverage model.
        if log10_coverage:
            logger.info("Converting coverage sites by log10.")
            cov_value = hl.log10(low_cov_group_ht.cov_value)
        else:
            cov_value = low_cov_group_ht.cov_value

        # Build the coverage model.
        coverage_model_expr = build_coverage_model(
            low_coverage_oe_expr=low_cov_group_ht.low_coverage_oe,
            coverage_expr=cov_value,
        )
        coverage_model = tuple(low_cov_group_ht.aggregate(coverage_model_expr).beta)
    else:
        coverage_model = None

    return coverage_model, plateau_models


def build_plateau_models(
    cpg_expr: hl.expr.BooleanExpression,
    mu_snp_expr: hl.expr.Float64Expression,
    observed_variants_expr: hl.expr.Int64Expression,
    possible_variants_expr: hl.expr.Int64Expression,
    weighted: bool = False,
    additional_grouping: Dict[str, hl.expr.StringExpression] = {},
) -> Dict[str, Union[Dict[bool, hl.expr.ArrayExpression], hl.ArrayExpression]]:
    """
    Build plateau models to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:
    - x: `mu_snp_expr`
    - y: `observed_variants_expr` / `possible_variants_expr`
    or `pops_observed_variants_array_expr`[index] / `possible_variants_expr`
    if `pops` is specified

    :param cpg_expr: BooleanExpression noting whether a site is a CPG site.
    :param mu_snp_expr: Float64Expression of the mutation rate.
    :param observed_variants_expr: Int64Expression of the observed variant counts.
    :param possible_variants_expr: Int64Expression of the possible variant counts.
    :param pops_observed_variants_array_expr: Nested ArrayExpression with all observed
        variant counts ArrayNumericExpressions for specified populations. e.g., `[[1,1,
        1],[1,1,1]]`. Default is None.
    :param weighted: Whether to generalize the model to weighted least squares using
        'possible_variants'. Default is False.
    :return: A dictionary of intercepts and slopes of plateau models. The keys are
        'total' (for all sites) and 'pop' (optional; for populations). The values for
        'total' is a dictionary (e.g., <DictExpression of type dict<bool,
        array<float64>>>), and the value for 'pop' is a nested list of dictionaries (e.
        g., <ArrayExpression of type array<array<dict<bool, array<float64>>>>>). The
        key of the dictionary in the nested list is CpG status (BooleanExpression), and
        the value is an ArrayExpression containing intercept and slope values.
    """
    grouping = hl.struct(cpg=cpg_expr, **additional_grouping)

    # Build plateau models for all sites
    plateau_models_agg_expr = hl.agg.group_by(
        grouping,
        hl.agg.array_agg(
            lambda x: hl.agg.linreg(
                x / possible_variants_expr,
                [1, mu_snp_expr],
                weight=possible_variants_expr if weighted else None,
            ).beta,
            observed_variants_expr,
        ),
    )

    return plateau_models_agg_expr


def build_coverage_model(
    low_coverage_oe_expr: hl.expr.Float64Expression,
    coverage_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Build coverage model.

    This function uses linear regression to build a model of coverage to correct
    proportion of expected variation at low coverage sites.

    The x and y of the coverage model:
    - x: `coverage_expr`
    - y: `low_coverage_oe_expr`

    :param low_coverage_oe_expr: The Float64Expression of observed:expected ratio
        for a given coverage level.
    :param coverage_expr: The Float64Expression of the coverage expression.
    :return: StructExpression with intercept and slope of the model.
    """
    return hl.agg.linreg(low_coverage_oe_expr, [1, coverage_expr])


def get_coverage_correction_expr(
    coverage_expr: hl.Int64Expression,
    coverage_model: Tuple[float, float],
    model_type_expr: hl.StringExpression,
    log10_coverage: bool = False,
) -> hl.Float64Expression:
    """
    Get the coverage correction factor.

    The coverage correction factor is computed as the coverage value adjusted by the
    coverage model. If coverage is 0, the coverage correction factor is set to 0. If
    coverage is above the high coverage cutoff, the coverage correction factor is set to
    1.

    :param coverage_expr: Expression for exome coverage.
    :param coverage_model: A linear model (output of `build_models()` in
        gnomad_methods), formatted as a Tuple of intercept and slope, that calibrates a
        given coverage level to observed:expected ratio. Default is None.
    :param high_coverage_cutoff: Exome coverage cutoff. Sites with coverage above this
        cutoff are considered well covered and do not require a coverage correction, so
        the coverage correction factor is set to 1. Sites below this cutoff have low
        coverage and require the coverage correction defined by the `coverage_model`.
        Default is `COVERAGE_CUTOFF`.
    :param log10_coverage: Whether to convert coverage sites with log10 when applying
        the coverage model. Default is True.
    :return: Coverage correction expression.
    """
    if log10_coverage:
        cov_corr_expr = hl.log10(coverage_expr)
    else:
        cov_corr_expr = coverage_expr

    return (
        hl.case()
        .when(coverage_expr == 0, 0)
        .when(
            model_type_expr == "low_covergage",
            coverage_model[1] * cov_corr_expr + coverage_model[0],
        )
        .default(1)
    )


def apply_plateau_models_per_variant(
    mu_expr: hl.Float64Expression,
    cpg_expr: hl.BooleanExpression,
    genomic_region_expr: hl.StringExpression,
    plateau_models_expr: hl.StructExpression,
) -> hl.Float64Expression:
    """
    Compute the predicted probability observed for a variant.

    The predicted probability observed is computed as the mutation rate adjusted by the
    plateau model for CpG transitions and non-CpG transitions. If a coverage correction
    factor is provided, the predicted probability observed is multiplied by the coverage
    correction factor.

    :param mu_expr: Mutation rate.
    :param cpg_expr: Boolean expression indicating whether the variant is a CpG
        transition.
    :param plateau_models_expr: Linear models that calibrate mutation rate to proportion
        observed for high coverage exome. It includes models for CpG sites, non-CpG
        sites, and each population in `POPS`.
    :return: Predicted probability observed for a variant.
    """

    def _apply_model(plateau_model):
        slope = plateau_model[1]
        intercept = plateau_model[0]
        ppo_expr = mu_expr * slope + intercept

        return ppo_expr

    plateau_models_expr = plateau_models_expr.get(
        hl.Struct(cpg=cpg_expr, genomic_region=genomic_region_expr)
    )

    return plateau_models_expr.map(lambda x: _apply_model(x))


def apply_models_per_variant(
    mu_expr: hl.Float64Expression,
    cpg_expr: hl.BooleanExpression,
    coverage_expr: hl.Int64Expression,
    genomic_region_expr: hl.StringExpression,
    model_type_expr: hl.StringExpression,
    plateau_models_expr: hl.StructExpression,
    coverage_model_expr: Optional[Tuple[float, float]] = None,
    high_coverage_cutoff: int = COVERAGE_CUTOFF,
    log10_coverage: bool = False,
) -> hl.StructExpression:
    """
    Compute the predicted probability observed for a variant.

    The predicted probability observed is computed as the mutation rate adjusted by the
    plateau model for CpG transitions and non-CpG transitions. If a coverage correction
    factor is provided, the predicted probability observed is multiplied by the coverage
    correction factor.

    :param mu_expr: Mutation rate.
    :param cpg_expr: Boolean expression indicating whether the variant is a CpG
        transition.
    :param plateau_models_expr: Linear models that calibrate mutation rate to proportion
        observed for high coverage exome. It includes models for CpG sites, non-CpG
        sites, and each population in `POPS`.
    :param coverage_model_expr: Optional Coverage model to calibrate low coverage sites.
        Default is None.
    :param coverage_expr: Expression for exome coverage.
    :param high_coverage_cutoff: Exome coverage cutoff. Sites with coverage above this
        cutoff are considered well covered and do not require a coverage correction, so
        the coverage correction factor is set to 1. Sites below this cutoff have low
        coverage and require the coverage correction defined by the `coverage_model`.
        Default is `COVERAGE_CUTOFF`.
    :param log10_coverage: Whether to convert coverage sites with log10 when applying
        the coverage model. Default is True.
    :return: Predicted probability observed for a variant.
    """
    if coverage_model_expr is None:
        cov_corr_expr = hl.if_else(coverage_expr == 0, 0, 1)
    else:
        cov_corr_expr = get_coverage_correction_expr(
            coverage_expr,
            coverage_model_expr,
            model_type_expr,
            log10_coverage=log10_coverage,
        )

    ppo_expr = apply_plateau_models_per_variant(
        mu_expr, cpg_expr, genomic_region_expr, plateau_models_expr
    )
    ppo_with_cov_corr_expr = ppo_expr * cov_corr_expr

    exp_ppo_no_coverage_no_int_expr = 1 - (1.004728 * hl.exp(-67885262 * mu_expr))
    exp_ppo_no_coverage_expr = 1 - (
        0.02763048 + 0.9799107 * hl.exp(-71166173 * mu_expr)
    )

    return hl.struct(
        ppo_no_coverage=ppo_expr,
        coverage_correction=cov_corr_expr,
        predicted_probability_observed=ppo_with_cov_corr_expr,
        exp_ppo_no_coverage_no_int=exp_ppo_no_coverage_no_int_expr,
        exp_ppo_no_coverage=exp_ppo_no_coverage_expr,
        exp_ppo_no_int=exp_ppo_no_coverage_no_int_expr * cov_corr_expr,
        exp_ppo=exp_ppo_no_coverage_expr * cov_corr_expr,
    )


def create_per_variant_expected_ht(
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Optional[Tuple[float, float]] = None,
    log10_coverage: bool = True,
    custom_vep_annotation: str = None,
    use_mane_select: bool = True,
) -> hl.Table:
    # Add necessary constraint annotations for grouping.
    include_canonical_group = False
    include_mane_select_group = False
    if custom_vep_annotation == "worst_csq_by_gene":
        vep_annotation = "worst_csq_by_gene"
        if use_mane_select:
            raise ValueError(
                "'mane_select' cannot be set to True when custom_vep_annotation is set"
                " to 'worst_csq_by_gene'."
            )
    else:
        vep_annotation = custom_vep_annotation
        include_canonical_group = True
        include_mane_select_group = use_mane_select

    context_ht = context_ht.filter(hl.is_defined(context_ht.apply_model))
    context_ht, _ = annotate_exploded_vep_for_constraint_groupings(
        ht=context_ht,
        coverage_expr=context_ht.exomes_coverage,
        vep_annotation=vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )
    context_ht = annotate_with_mu(context_ht, mutation_ht)

    ppo_expr = apply_models_per_variant(
        context_ht.mu_snp,
        context_ht.cpg,
        context_ht.coverage,
        context_ht.genomic_region,
        context_ht.apply_model,
        plateau_models,
        coverage_model_expr=coverage_model,
        log10_coverage=log10_coverage,
    )

    context_ht = context_ht.annotate(**ppo_expr)
    coverage_model_global = coverage_model if coverage_model else "None"
    context_ht = context_ht.annotate_globals(
        apply_model_params=hl.struct(
            plateau_models=plateau_models,
            coverage_model=coverage_model_global,
            log10_coverage=log10_coverage,
        )
    )

    return context_ht


def aggregate_per_variant_expected_ht_old(
    ht,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Optional[Tuple[float, float]] = None,
    log10_coverage: bool = True,
    additional_grouping: Tuple[str] = (),
    include_genomic_constraint_adjustment: bool = False,
):
    # ht = ht.annotate(
    #    transcript=ht.original_transcript_consequences.transcript_id,
    # )
    # weighted_sum_exprs = {
    #    "base_mu_snp": ht.mu_snp,
    #    "cov_corr": ht.coverage_correction,
    #    "base_ppo": ht.ppo_no_coverage,
    #    "base_cov_corr_ppo": ht.predicted_probability_observed,
    # }

    second_grouping = (
        "annotation",
        "modifier",
        "lof_flags",
        "gene",
        "gene_id",
        "transcript",
    ) + additional_grouping
    # first_grouping = (
    #    "methylation_level",
    #    "apply_model",
    #    "cpg",
    #    "mutation_type",
    #    "coverage",
    # ) + second_grouping
    # ht = create_observed_and_possible_ht(
    #    ht,
    #    mutation_ht,
    #    additional_grouping=first_grouping,
    #    weight_exprs=weighted_sum_exprs,
    # )
    # ht.describe()
    # ht = ht.checkpoint(
    #    new_temp_file(prefix="temp_agg_of_per_variant_old", extension="ht")
    # )
    ht = hl.read_table(
        "gs://gnomad-tmp-4day/temp_agg_of_per_variant_old-Wo1d8mBvIrbnc4bY8ppRLd.ht",
        _n_partitions=5000,
    )
    # ht.show()

    ht = annotate_with_mu(ht, mutation_ht)
    ht.show()

    ppo_expr = apply_models_per_variant(
        ht.mu_snp,
        ht.cpg,
        ht.coverage,
        ht.genomic_region,
        ht.apply_model,
        plateau_models,
        coverage_model_expr=coverage_model,
        log10_coverage=log10_coverage,
    )

    ht = ht.annotate(
        mu_snp=ht.mu_snp * ht.possible_variants,
        ppo=ppo_expr.ppo_no_coverage * ht.possible_variants,
        cov_corr_ppo=ppo_expr.predicted_probability_observed * ht.possible_variants,
    )

    ht.describe()
    ht = ht.checkpoint(
        new_temp_file(prefix="temp_agg_of_per_variant2_old", extension="ht")
    )
    ht.show()

    additional_agg_sum_exprs = [
        "mu_snp",
        "ppo",
        "cov_corr_ppo",
        "base_mu_snp",
        "cov_corr",
        "base_ppo",
        "base_cov_corr_ppo",
    ]
    agg_expr = {
        "observed_variants": hl.agg.array_sum(ht.observed_variants),
        "possible_variants": hl.agg.sum(ht.possible_variants),
        **{
            k: (
                hl.agg.array_sum(ht[k])
                if isinstance(ht[k], hl.ArrayExpression)
                else hl.agg.sum(ht[k])
            )
            for k in additional_agg_sum_exprs
        },
    }

    # Apply each variant count aggregation in `agg_expr` to get counts for all
    # combinations of `grouping`.
    ht = ht.group_by(*second_grouping).partition_hint(1000).aggregate(**agg_expr)
    ht.describe()
    ht = ht.checkpoint(
        new_temp_file(prefix="temp_agg_of_per_variant3_old", extension="ht")
    )
    ht.show()

    # Compute the observed:expected ratio.
    # ht = ht.annotate(
    #    mu=[ht.mu_snp * x for x in mu_expr],
    #    expected_variants=exp_expr,
    # )
    # ht = ht.annotate_globals(
    #    mu_meta=hl.literal(mu_meta, dtype="array<array<str>>"),
    #    exp_meta=hl.literal(exp_meta, dtype="array<array<str>>"),
    # )

    return ht


def aggregate_per_variant_expected_ht(
    ht,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Optional[Tuple[float, float]] = None,
    log10_coverage: bool = True,
    additional_grouping: Tuple[str] = (),
    include_genomic_constraint_adjustment: bool = False,
):
    # weighted_sum_exprs = {
    #    "base_mu_snp": ht.mu_snp,
    #    "cov_corr": ht.coverage_correction,
    #    "base_ppo": ht.ppo_no_coverage,
    #    "base_cov_corr_ppo": ht.predicted_probability_observed,
    #    "base_exp_ppo_no_coverage_no_int": ht.exp_ppo_no_coverage_no_int,
    #    "base_exp_ppo_no_coverage": ht.exp_ppo_no_coverage,
    #    "base_exp_ppo_no_int": ht.exp_ppo_no_int,
    #    "base_exp_ppo": ht.exp_ppo,
    # }

    # if include_genomic_constraint_adjustment:
    #    adj_r_expr = hl.or_else(ht.adj_r, 1)
    #    weighted_sum_exprs.update(
    #        {
    #            "adj_r": ht.adj_r,
    #            "cov_corr_adj_r": ht.coverage_correction * adj_r_expr,
    #            "base_ppo_adj_r": ht.ppo_no_coverage * adj_r_expr,
    #            "base_cov_corr_ppo_adj_r": ht.predicted_probability_observed * adj_r_expr,
    #            "base_exp_ppo_no_coverage_no_int_adj_r": ht.exp_ppo_no_coverage_no_int * adj_r_expr,
    #            "base_exp_ppo_no_coverage_adj_r": ht.exp_ppo_no_coverage * adj_r_expr,
    #            "base_exp_ppo_no_int_adj_r": ht.exp_ppo_no_int * adj_r_expr,
    #            "base_exp_ppo_adj_r": ht.exp_ppo * adj_r_expr,
    #        }
    #    )

    second_grouping = (
        "annotation",
        "modifier",
        "lof_flags",
        "gene",
        "gene_id",
        "transcript",
        "canonical",
        "mane_select",
    ) + additional_grouping
    # first_grouping = ("methylation_level", "apply_model", "cpg", "mutation_type", "coverage") + second_grouping
    # ht = create_observed_and_possible_ht(
    #    ht,
    #    mutation_ht,
    #    additional_grouping=first_grouping,
    #    weight_exprs=weighted_sum_exprs,
    # )
    # ht.describe()
    # ht = ht.checkpoint(new_temp_file(prefix="temp_agg_of_per_variant", extension="ht"))
    ht = hl.read_table(
        "gs://gnomad-tmp-4day/temp_agg_of_per_variant-sBNtW9xtMHxu1W7fEJouis.ht",
        _n_partitions=5000,
    )
    ht.show()

    ht = annotate_with_mu(ht, mutation_ht)
    ht.show()

    ppo_expr = apply_models_per_variant(
        ht.mu_snp,
        ht.cpg,
        ht.coverage,
        ht.genomic_region,
        ht.apply_model,
        plateau_models,
        coverage_model_expr=coverage_model,
        log10_coverage=log10_coverage,
    )

    ht = ht.annotate(
        mu_snp=ht.mu_snp * ht.possible_variants,
        ppo=ppo_expr.ppo_no_coverage * ht.possible_variants,
        cov_corr_ppo=ppo_expr.predicted_probability_observed * ht.possible_variants,
        exp_ppo_no_coverage_no_int=ppo_expr.exp_ppo_no_coverage_no_int
        * ht.possible_variants,
        exp_ppo_no_coverage=ppo_expr.exp_ppo_no_coverage * ht.possible_variants,
        exp_ppo_no_int=ppo_expr.exp_ppo_no_int * ht.possible_variants,
        exp_ppo=ppo_expr.exp_ppo * ht.possible_variants,
    )

    ht.describe()
    ht = ht.checkpoint(new_temp_file(prefix="temp_agg_of_per_variant2", extension="ht"))
    ht.show()

    # ht = hl.read_table("gs://gnomad-tmp-4day/temp_agg_of_per_variant2-F8qM41isz2RhJX3LztNjTb.ht", _n_partitions=5000)

    additional_agg_sum_exprs = [
        "mu_snp",
        "ppo",
        "cov_corr_ppo",
        "exp_ppo_no_coverage_no_int",
        "exp_ppo_no_coverage",
        "exp_ppo_no_int",
        "exp_ppo",
        "base_mu_snp",
        "cov_corr",
        "base_ppo",
        "base_cov_corr_ppo",
        "base_exp_ppo_no_coverage_no_int",
        "base_exp_ppo_no_coverage",
        "base_exp_ppo_no_int",
        "base_exp_ppo",
        "adj_r",
        "cov_corr_adj_r",
        "base_ppo_adj_r",
        "base_cov_corr_ppo_adj_r",
        "base_exp_ppo_no_coverage_no_int_adj_r",
        "base_exp_ppo_no_coverage_adj_r",
        "base_exp_ppo_no_int_adj_r",
        "base_exp_ppo_adj_r",
    ]
    agg_expr = {
        "observed_variants": hl.agg.array_sum(ht.observed_variants),
        "possible_variants": hl.agg.sum(ht.possible_variants),
        **{
            k: (
                hl.agg.array_sum(ht[k])
                if isinstance(ht[k], hl.ArrayExpression)
                else hl.agg.sum(ht[k])
            )
            for k in additional_agg_sum_exprs
        },
    }

    # Apply each variant count aggregation in `agg_expr` to get counts for all
    # combinations of `grouping`.
    ht = ht.group_by(*second_grouping).partition_hint(1000).aggregate(**agg_expr)

    ht.describe()
    ht = ht.checkpoint(new_temp_file(prefix="temp_agg_of_per_variant3", extension="ht"))
    ht.show()

    # mu_meta = [[], ["cov_corr"], ["adj_r"], ["cov_corr", "adj_r"]]
    # mu_expr = [1, ht.cov_corr, ht.adj_r, ht.cov_corr_adj_r]
    # exp_meta = [["plateau", "per_base"] + x for x in mu_meta]
    # exp_meta += [["exp no int", "per_base"] + x for x in mu_meta]
    # exp_meta += [["exp with int", "per_base"] + x for x in mu_meta]
    # exp_expr = [ht.base_ppo, ht.base_cov_corr_ppo, ht.base_ppo_adj_r, ht.base_cov_corr_ppo_adj_r]
    # exp_expr += [ht.base_exp_ppo_no_coverage_no_int, ht.base_exp_ppo_no_int, ht.base_exp_ppo_no_coverage_no_int_adj_r, ht.base_exp_ppo_no_int_adj_r]
    # exp_expr += [ht.base_exp_ppo_no_coverage, ht.base_exp_ppo, ht.base_exp_ppo_no_coverage_adj_r, ht.base_exp_ppo_adj_r]

    adj_r_ht = hl.read_table(
        "gs://gnomad/v4.1/constraint/resources/ncc_adj_r_by_transcript_WG.ht"
    )
    adj_r_keyed = adj_r_ht[ht.transcript]
    # add_tx_adj_r = [["plateau"], ["plateau", "cov_corr"]]
    # add_tx_adj_r_idx = [exp_meta.index(x) for x in add_tx_adj_r]
    ht = ht.annotate(
        tx_adj_r=adj_r_keyed.adj_r,
    )

    # Compute the observed:expected ratio.
    # ht = ht.annotate(
    #    mu=[ht.mu_snp * x for x in mu_expr],
    #    expected_variants=exp_expr,
    # )
    # ht = ht.annotate_globals(
    #    mu_meta=hl.literal(mu_meta, dtype="array<array<str>>"),
    #    exp_meta=hl.literal(exp_meta, dtype="array<array<str>>"),
    # )

    return ht.naive_coalesce(1000)


def calculate_mu_by_downsampling(
    ht: hl.Table,
    additional_grouping: Tuple[str] = ("methylation_level",),
    total_mu: float = 1.2e-08,
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

    :param ht: Context Table for autosome/pseudoautosomal regions.
    :param additional_grouping: Annotations other than 'context', 'ref', and 'alt'.
        Default is ('methylation_level',).
    :param ac_cutoff: The cutoff of allele count when filtering context Table
        and genome sites Table.
    :param total_mu: The per-generation mutation rate. Default is 1.2e-08.
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
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by context, ref, alt, and 'additional_grouping'.
    ht = count_observed_and_possible_by_group(
        ht,
        ht.mu.possible,
        ht.mu.observed,
        additional_grouping=additional_grouping,
    )

    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants)) // 3
    logger.info(
        "Total bases to use when calculating correction_factors: %f", total_bases
    )

    # Compute the proportion observed, which represents the relative mutability of each
    # variant class.
    po_expr = ht.observed_variants / ht.possible_variants
    correction_factors = ht.aggregate(
        total_mu / (hl.agg.array_sum(ht.observed_variants) / total_bases),
        _localize=False,
    )
    mu_expr = correction_factors * ht.observed_variants / ht.possible_variants
    ht = ht.annotate(
        proportion_observed=po_expr,
        mu=mu_expr,
        mu_snp=mu_expr[ht.calculate_mu_globals.downsampling_idx],
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
    keys = list(ht.key)
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

    ms_ht = ms_ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    ms_ht.describe()
    lof_fields = [x for x in ht.row_value if x.startswith("lof")]

    # Determine decile bins.
    n_transcripts = ms_ht.count()

    for lof_field in lof_fields:
        # Rank lof.oe_ci.upper in ascending order.
        cov_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        cov_ms_ht = cov_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=cov_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(cov_ms_ht.upper_rank * 10 / n_transcripts),
                )
            }
        )
        cov_ms_ht = cov_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.1.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        no_cov_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci_no_cov_corr[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        no_cov_ms_ht = no_cov_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=no_cov_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(
                        no_cov_ms_ht.upper_rank * 10 / n_transcripts
                    ),
                )
            }
        )
        no_cov_ms_ht = no_cov_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.2.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        cov_adj_r_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci_adj_r[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        cov_adj_r_ms_ht = cov_adj_r_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=cov_adj_r_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(
                        cov_adj_r_ms_ht.upper_rank * 10 / n_transcripts
                    ),
                )
            }
        )
        cov_adj_r_ms_ht = cov_adj_r_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.3.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        no_cov_adj_r_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci_adj_r_no_cov_corr[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        no_cov_adj_r_ms_ht = no_cov_adj_r_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=no_cov_adj_r_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(
                        no_cov_adj_r_ms_ht.upper_rank * 10 / n_transcripts
                    ),
                )
            }
        )
        no_cov_adj_r_ms_ht = no_cov_adj_r_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.4.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        cov_tx_adj_r_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci_tx_adj_r[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        cov_tx_adj_r_ms_ht = cov_tx_adj_r_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=cov_tx_adj_r_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(
                        cov_tx_adj_r_ms_ht.upper_rank * 10 / n_transcripts
                    ),
                )
            }
        )
        cov_tx_adj_r_ms_ht = cov_tx_adj_r_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.5.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        no_cov_tx_adj_r_ms_ht = (
            ms_ht.order_by(ms_ht[lof_field].oe_ci_no_cov_corr_tx_adj_r[0].upper)
            .add_index(name="upper_rank")
            .key_by(*list(keys))
        )
        # Determine decile bins.
        no_cov_tx_adj_r_ms_ht = no_cov_tx_adj_r_ms_ht.select(
            **{
                f"oe_upper_rank_and_bin_{lof_field}": hl.struct(
                    upper_rank=no_cov_tx_adj_r_ms_ht.upper_rank,
                    upper_bin_decile=hl.int(
                        no_cov_tx_adj_r_ms_ht.upper_rank * 10 / n_transcripts
                    ),
                )
            }
        )
        no_cov_tx_adj_r_ms_ht = no_cov_tx_adj_r_ms_ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.6.ht",
            # _read_if_exists=True,
            overwrite=True,
        )

        # Add rank and bin annotations back to original Table.
        cov_ms_index = cov_ms_ht[ht.key]
        no_cov_ms_index = no_cov_ms_ht[ht.key]
        cov_adj_r_ms_index = cov_adj_r_ms_ht[ht.key]
        no_cov_adj_r_ms_index = no_cov_adj_r_ms_ht[ht.key]
        cov_tx_adj_r_ms_index = cov_tx_adj_r_ms_ht[ht.key]
        no_cov_tx_adj_r_ms_index = no_cov_tx_adj_r_ms_ht[ht.key]

        ht = ht.annotate(
            **{
                lof_field: ht[lof_field].annotate(
                    global_oe_ci=ht[lof_field]
                    .oe_ci[0]
                    .annotate(**cov_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]),
                    global_oe_ci_no_cov_corr=ht[lof_field]
                    .oe_ci_no_cov_corr[0]
                    .annotate(**no_cov_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]),
                    global_oe_ci_adj_r=ht[lof_field]
                    .oe_ci_adj_r[0]
                    .annotate(
                        **cov_adj_r_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]
                    ),
                    global_oe_ci_adj_r_no_cov_corr=ht[lof_field]
                    .oe_ci_adj_r_no_cov_corr[0]
                    .annotate(
                        **no_cov_adj_r_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]
                    ),
                    global_oe_ci_tx_adj_r=ht[lof_field]
                    .oe_ci_tx_adj_r[0]
                    .annotate(
                        **cov_tx_adj_r_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]
                    ),
                    global_oe_ci_no_cov_corr_tx_adj_r=ht[lof_field]
                    .oe_ci_no_cov_corr_tx_adj_r[0]
                    .annotate(
                        **no_cov_tx_adj_r_ms_index[f"oe_upper_rank_and_bin_{lof_field}"]
                    ),
                )
            }
        )
        ht = ht.checkpoint(
            f"gs://gnomad-tmp-4day/constraint/oe_lof_upper_rank_and_bin.{lof_field}.all.ht",
            # _read_if_exists=True,
            overwrite=True,
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
    mu_meta = [[], ["adj_r"], ["cov_corr", "adj_r"]]
    mu_expr = [1, ht.adj_r, ht.cov_corr_adj_r]

    exp_meta = [["exp no int"], ["exp no int", "cov_corr"]]
    exp_meta += [["exp with int"], ["exp with int", "cov_corr"]]

    exp_meta += [
        ["exp no int", "per_base"],
        ["exp no int", "per_base", "cov_corr"],
        ["exp no int", "per_base", "adj_r"],
        ["exp no int", "per_base", "cov_corr", "adj_r"],
    ]
    exp_meta += [
        ["exp with int", "per_base"],
        ["exp with int", "per_base", "cov_corr"],
        ["exp with int", "per_base", "adj_r"],
        ["exp with int", "per_base", "cov_corr", "adj_r"],
    ]

    exp_expr = [ht.exp_ppo_no_coverage_no_int, ht.exp_ppo_no_int]
    exp_expr += [ht.exp_ppo_no_coverage, ht.exp_ppo]

    exp_expr += [
        ht.base_exp_ppo_no_coverage_no_int,
        ht.base_exp_ppo_no_int,
        ht.base_exp_ppo_no_coverage_no_int_adj_r,
        ht.base_exp_ppo_no_int_adj_r,
    ]
    exp_expr += [
        ht.base_exp_ppo_no_coverage,
        ht.base_exp_ppo,
        ht.base_exp_ppo_no_coverage_adj_r,
        ht.base_exp_ppo_adj_r,
    ]

    array_exp_meta = [["plateau"], ["plateau", "cov_corr"]]
    array_exp_meta += [
        ["plateau", "per_base"],
        ["plateau", "per_base", "cov_corr"],
        ["plateau", "per_base", "adj_r"],
        ["plateau", "per_base", "cov_corr", "adj_r"],
    ]

    array_exp_expr = [ht.ppo, ht.cov_corr_ppo]
    array_exp_expr += [
        ht.base_ppo,
        ht.base_cov_corr_ppo,
        ht.base_ppo_adj_r,
        ht.base_cov_corr_ppo_adj_r,
    ]

    # Compute the observed:expected ratio.
    ht = ht.transmute(
        adj_r=hl.or_missing(ht.adj_r > 0, ht.adj_r),
        mu=[ht.mu_snp * x for x in mu_expr],
        expected_variants=exp_expr,
        expected_variants_arrays=array_exp_expr,
    )
    ht = ht.annotate_globals(
        mu_meta=hl.literal(mu_meta, dtype="array<array<str>>"),
        exp_meta=hl.literal(exp_meta, dtype="array<array<str>>"),
        array_exp_meta=hl.literal(array_exp_meta, dtype="array<array<str>>"),
    )

    # ht = ht.checkpoint(
    #    new_temp_file(prefix="compute_constraint_metrics", extension="ht")
    # )

    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}

    lof_classic_expr = hl.literal(set(classic_lof_annotations)).contains(ht.annotation)
    lof_hc_expr = ht.modifier == "HC"
    lof_hc_lc_expr = lof_hc_expr | (ht.modifier == "LC")
    mis_expr = ht.annotation == "missense_variant"
    annotation_dict = {
        "csq_set": {"syn": ht.annotation == "synonymous_variant", "mis": mis_expr},
        # Filter to missense variants that are in the alpha missense list.
        "alphamis": {
            "98_per": mis_expr & (hl.str(ht.am_per_98) == "true"),
            "0.999": mis_expr & (hl.str(ht.am_over_0_999) == "true"),
            "99_per": mis_expr & (hl.str(ht.am_per_99) == "true"),
        },
        "lof": {
            # Filter to classic LoF annotations.
            "classic": lof_classic_expr,
            # Filter to LOFTEE HC or LC.
            "hc_lc": lof_hc_lc_expr,
            # Filter to classic LoF annotations with LOFTEE HC or LC.
            "classic_hc_lc": lof_classic_expr & lof_hc_lc_expr,
            # Filter to LoF annotations with LOFTEE HC.
            "hc": lof_hc_expr,
        },
    }

    meta = generate_filter_combinations(
        [["csq_set"], ["alphamis"], ["lof"], ["lof", "alphamis"]],
        {k: list(v.keys()) for k, v in annotation_dict.items()},
    )

    meta_filter_expr = []
    for m in meta:
        if "csq_set" in m:
            meta_filter_expr.append(annotation_dict["csq_set"][m["csq_set"]])
        elif "lof" in m and "alphamis" in m:
            meta_filter_expr.append(
                annotation_dict["lof"][m["lof"]]
                | annotation_dict["alphamis"][m["alphamis"]]
            )
        elif "lof" in m:
            meta_filter_expr.append(annotation_dict["lof"][m["lof"]])
        elif "alphamis" in m:
            meta_filter_expr.append(annotation_dict["alphamis"][m["alphamis"]])

    ht = ht.annotate(constraint_group_filters=meta_filter_expr)
    ht = ht.annotate_globals(constraint_meta=meta)
    print("One")
    ht.describe()
    print()
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.constraint_group_filters.ht",
        _read_if_exists=True,
        # overwrite=True,
    )
    common_ann = [
        "adj_r",
        "possible_variants",
    ]
    agg_expr = {
        "mu": hl.agg.array_agg(
            lambda f: hl.agg.filter(f, hl.agg.array_sum(ht.mu)),
            ht.constraint_group_filters,
        ),
        **{
            ann: hl.agg.array_agg(
                lambda f: hl.agg.filter(f, hl.agg.sum(ht[ann])),
                ht.constraint_group_filters,
            )
            for ann in common_ann
        },
    }

    ht = ht.group_by(*keys).aggregate(
        **agg_expr,
        observed_variants=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, hl.agg.array_sum(ht.observed_variants)),
            ht.constraint_group_filters,
        ),
        expected_variants=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, hl.agg.array_sum(ht.expected_variants)),
            ht.constraint_group_filters,
        ),
        expected_variants_arrays=hl.agg.array_agg(
            lambda f: hl.agg.filter(
                f,
                hl.agg.array_agg(
                    lambda x: hl.agg.array_sum(x), ht.expected_variants_arrays
                ),
            ),
            ht.constraint_group_filters,
        ),
    )
    print("Two")
    ht.describe()
    print()
    # ht = ht.checkpoint(
    ht = hl.read_table(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.global_annotations.ht",
        _n_partitions=5000,
        #    _read_if_exists=True,
        # overwrite=True,
    )

    ht = ht.annotate(
        adj_r=ht.adj_r.map(lambda x: hl.or_missing(x > 0, x)),
        oe=hl.map(
            lambda obs, exp: exp.map(lambda e: divide_null(obs[0], e)),
            ht.observed_variants,
            ht.expected_variants,
        ),
        oe_ci=hl.map(
            lambda obs, exp: exp.map(lambda e: oe_confidence_interval(obs[0], e)),
            ht.observed_variants,
            ht.expected_variants,
        ),
        z_raw=hl.map(
            lambda obs, exp: exp.map(lambda e: calculate_raw_z_score(obs[0], e)),
            ht.observed_variants,
            ht.expected_variants,
        ),
        oe_arrays=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: divide_null(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants_arrays,
        ),
        oe_ci_arrays=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: oe_confidence_interval(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants_arrays,
        ),
        z_raw_arrays=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: calculate_raw_z_score(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants_arrays,
        ),
    )

    print("Three")
    ht.describe()
    print()
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.oe_agg.pop.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    # Filter to only rows with at least 1 obs or exp across all keys in annotation_dict.
    ht = ht.filter(
        hl.sum(
            ht.total.map(
                lambda x: hl.or_else(x.obs, 0)
                + hl.sum(x.exp.map(lambda y: hl.or_else(y, 0)))
            )
        )
        > 0
    )

    """
    # Create dictionary with outlier z-score thresholds with annotation as key
    # and list of thresholds [lower, upper] as values.
    lof_ann.extend(
        ["mis_am_per_98", "mis_am_0_999", "mis_am_per_99"]
        + [f"{k}_am_per_98" for k in lof_ann]
        + [f"{k}_am_0_999" for k in lof_ann]
        + [f"{k}_am_per_99" for k in lof_ann]
    )
    z_score_outlier_dict = {
        **{k: [raw_z_outlier_threshold_lower_lof, None] for k in lof_ann},
        "mis": [raw_z_outlier_threshold_lower_missense, None],
        "syn": [raw_z_outlier_threshold_lower_syn, raw_z_outlier_threshold_upper_syn],
    }

    # The constraint_flags dict is used to filter the final ht.constraint_flags
    # annotation to the flags that should be considered in the z-score 'sd'
    # computation of the specified ann.
    constraint_flags = {
        ann: get_constraint_flags(
            exp_expr=ht[ann].exp[0],
            raw_z_expr=ht[ann].z_raw[0],
            raw_z_lower_threshold=z_score_outlier_dict[ann][0],
            raw_z_upper_threshold=z_score_outlier_dict[ann][1],
            flag_postfix=ann,
        )
        for ann in annotation_dict
    }

    # Add a 'no_variants' flag indicating that there are zero observed variants summed
    # across pLoF, missense, and synonymous variants.
    constraint_flags_expr = {
        "no_variants": hl.sum(
            [hl.or_else(ht[ann].obs[0], 0) for ann in annotation_dict]
        )
        == 0,
    }
    for ann in annotation_dict:
        constraint_flags_expr.update(constraint_flags[ann])

    ht = ht.annotate(constraint_flags=add_filters_expr(filters=constraint_flags_expr))
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.constraint_flags.pop.ht",
        #_read_if_exists=True,
        overwrite=True,
    )

    # Add z-score 'sd' annotation to globals.
    ht = ht.annotate_globals(
        sd_raw_z=ht.aggregate(
            hl.struct(
                **{
                    ann: calculate_raw_z_score_sd(
                        raw_z_expr=ht[ann].z_raw[0],
                        flag_expr=ht.constraint_flags.intersection(
                            constraint_flags[ann].keys() | {"no_variants"}
                        ),
                        mirror_neg_raw_z=(ann != "syn"),
                    )
                    for ann in annotation_dict
                }
            )
        )
    )

    # Compute z-score from raw z-score and standard deviations.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(z_score=ht[ann].z_raw / ht.sd_raw_z[ann])
            for ann in annotation_dict
        }
    )

    # Compute the rank and decile of the lof oe upper confidence
    # interval for MANE Select or canonical ensembl transcripts.
    ht = add_oe_lof_upper_rank_and_bin(
        ht, use_mane_select_over_canonical=use_mane_select_over_canonical
    )
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.oe_lof_upper_rank_and_bin.pop.ht",
        #_read_if_exists=True,
        overwrite=True,
    )
    """
    # Add transcript annotations from GENCODE.
    ht = add_gencode_transcript_annotations(ht, gencode_ht)

    return ht.naive_coalesce(1000)


def compute_constraint_metrics_old(
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
    exp_meta = [["plateau"], ["plateau", "cov_corr"]]
    exp_meta += [["plateau", "per_base"], ["plateau", "per_base", "cov_corr"]]

    exp_expr = [ht.ppo, ht.cov_corr_ppo]
    exp_expr += [ht.base_ppo, ht.base_cov_corr_ppo]

    # Compute the observed:expected ratio.
    ht = ht.transmute(
        mu=ht.mu_snp,
        expected_variants=exp_expr,
    )
    ht = ht.annotate_globals(
        exp_meta=hl.literal(exp_meta, dtype="array<array<str>>"),
    )

    # ht = ht.checkpoint(
    #    new_temp_file(prefix="compute_constraint_metrics_old", extension="ht")
    # )

    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}

    lof_classic_expr = hl.literal(set(classic_lof_annotations)).contains(ht.annotation)
    lof_hc_expr = ht.modifier == "HC"
    lof_hc_lc_expr = lof_hc_expr | (ht.modifier == "LC")
    mis_expr = ht.annotation == "missense_variant"
    annotation_dict = {
        "csq_set": {"syn": ht.annotation == "synonymous_variant", "mis": mis_expr},
        # Filter to missense variants that are in the alpha missense list.
        "alphamis": {
            "98_per": mis_expr & (hl.str(ht.am_per_98) == "true"),
            "0.999": mis_expr & (hl.str(ht.am_over_0_999) == "true"),
            "99_per": mis_expr & (hl.str(ht.am_per_99) == "true"),
        },
        "lof": {
            # Filter to classic LoF annotations.
            "classic": lof_classic_expr,
            # Filter to LOFTEE HC or LC.
            "hc_lc": lof_hc_lc_expr,
            # Filter to classic LoF annotations with LOFTEE HC or LC.
            "classic_hc_lc": lof_classic_expr & lof_hc_lc_expr,
            # Filter to LoF annotations with LOFTEE HC.
            "hc": lof_hc_expr,
        },
    }

    meta = generate_filter_combinations(
        [["csq_set"], ["alphamis"], ["lof"], ["lof", "alphamis"]],
        {k: list(v.keys()) for k, v in annotation_dict.items()},
    )

    meta_filter_expr = []
    for m in meta:
        if "csq_set" in m:
            meta_filter_expr.append(annotation_dict["csq_set"][m["csq_set"]])
        elif "lof" in m and "alphamis" in m:
            meta_filter_expr.append(
                annotation_dict["lof"][m["lof"]]
                | annotation_dict["alphamis"][m["alphamis"]]
            )
        elif "lof" in m:
            meta_filter_expr.append(annotation_dict["lof"][m["lof"]])
        elif "alphamis" in m:
            meta_filter_expr.append(annotation_dict["alphamis"][m["alphamis"]])

    ht = ht.annotate(constraint_group_filters=meta_filter_expr)
    ht = ht.annotate_globals(constraint_meta=meta)
    print("One")
    ht.describe()
    print()
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.constraint_group_filters.old.ht",
        _read_if_exists=True,
        # overwrite=True,
    )
    common_ann = [
        "possible_variants",
    ]
    agg_expr = {
        "mu": hl.agg.array_agg(
            lambda f: hl.agg.filter(f, hl.agg.sum(ht.mu)),
            ht.constraint_group_filters,
        ),
        **{
            ann: hl.agg.array_agg(
                lambda f: hl.agg.filter(f, hl.agg.sum(ht[ann])),
                ht.constraint_group_filters,
            )
            for ann in common_ann
        },
    }

    ht = ht.group_by(*keys).aggregate(
        **agg_expr,
        observed_variants=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, hl.agg.array_sum(ht.observed_variants)),
            ht.constraint_group_filters,
        ),
        expected_variants=hl.agg.array_agg(
            lambda f: hl.agg.filter(
                f, hl.agg.array_agg(lambda x: hl.agg.array_sum(x), ht.expected_variants)
            ),
            ht.constraint_group_filters,
        ),
    )
    print("Two")
    ht.describe()
    print()
    # ht = ht.checkpoint(
    ht = hl.read_table(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.global_annotations.old.ht",
        _n_partitions=5000,
        #    _read_if_exists=True,
        #    overwrite=True,
    )

    ht = ht.annotate(
        oe=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: divide_null(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants,
        ),
        oe_ci=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: oe_confidence_interval(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants,
        ),
        z_raw=hl.map(
            lambda obs, exp_a: exp_a.map(
                lambda exp: hl.map(lambda o, e: calculate_raw_z_score(o, e), obs, exp)
            ),
            ht.observed_variants,
            ht.expected_variants,
        ),
    )

    print("Three")
    ht.describe()
    print()
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.oe_agg.pop.old.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    # Filter to only rows with at least 1 obs or exp across all keys in annotation_dict.
    ht = ht.filter(
        hl.sum(
            ht.total.map(
                lambda x: hl.or_else(x.obs, 0)
                + hl.sum(x.exp.map(lambda y: hl.or_else(y, 0)))
            )
        )
        > 0
    )

    """
    # Create dictionary with outlier z-score thresholds with annotation as key
    # and list of thresholds [lower, upper] as values.
    lof_ann.extend(
        ["mis_am_per_98", "mis_am_0_999", "mis_am_per_99"]
        + [f"{k}_am_per_98" for k in lof_ann]
        + [f"{k}_am_0_999" for k in lof_ann]
        + [f"{k}_am_per_99" for k in lof_ann]
    )
    z_score_outlier_dict = {
        **{k: [raw_z_outlier_threshold_lower_lof, None] for k in lof_ann},
        "mis": [raw_z_outlier_threshold_lower_missense, None],
        "syn": [raw_z_outlier_threshold_lower_syn, raw_z_outlier_threshold_upper_syn],
    }

    # The constraint_flags dict is used to filter the final ht.constraint_flags
    # annotation to the flags that should be considered in the z-score 'sd'
    # computation of the specified ann.
    constraint_flags = {
        ann: get_constraint_flags(
            exp_expr=ht[ann].exp[0],
            raw_z_expr=ht[ann].z_raw[0],
            raw_z_lower_threshold=z_score_outlier_dict[ann][0],
            raw_z_upper_threshold=z_score_outlier_dict[ann][1],
            flag_postfix=ann,
        )
        for ann in annotation_dict
    }

    # Add a 'no_variants' flag indicating that there are zero observed variants summed
    # across pLoF, missense, and synonymous variants.
    constraint_flags_expr = {
        "no_variants": hl.sum(
            [hl.or_else(ht[ann].obs[0], 0) for ann in annotation_dict]
        )
        == 0,
    }
    for ann in annotation_dict:
        constraint_flags_expr.update(constraint_flags[ann])

    ht = ht.annotate(constraint_flags=add_filters_expr(filters=constraint_flags_expr))
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.constraint_flags.pop.ht",
        #_read_if_exists=True,
        overwrite=True,
    )

    # Add z-score 'sd' annotation to globals.
    ht = ht.annotate_globals(
        sd_raw_z=ht.aggregate(
            hl.struct(
                **{
                    ann: calculate_raw_z_score_sd(
                        raw_z_expr=ht[ann].z_raw[0],
                        flag_expr=ht.constraint_flags.intersection(
                            constraint_flags[ann].keys() | {"no_variants"}
                        ),
                        mirror_neg_raw_z=(ann != "syn"),
                    )
                    for ann in annotation_dict
                }
            )
        )
    )

    # Compute z-score from raw z-score and standard deviations.
    ht = ht.annotate(
        **{
            ann: ht[ann].annotate(z_score=ht[ann].z_raw / ht.sd_raw_z[ann])
            for ann in annotation_dict
        }
    )

    # Compute the rank and decile of the lof oe upper confidence
    # interval for MANE Select or canonical ensembl transcripts.
    ht = add_oe_lof_upper_rank_and_bin(
        ht, use_mane_select_over_canonical=use_mane_select_over_canonical
    )
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/constraint/constraint_metrics.oe_lof_upper_rank_and_bin.pop.ht",
        #_read_if_exists=True,
        overwrite=True,
    )
    """
    # Add transcript annotations from GENCODE.
    ht = add_gencode_transcript_annotations(ht, gencode_ht)

    return ht.naive_coalesce(1000)


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
