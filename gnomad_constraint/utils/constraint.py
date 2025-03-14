"""Script containing utility functions used in the constraint pipeline."""

import functools
import logging
import operator
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
import numpy as np
from gnomad.assessment.summary_stats import generate_filter_combinations
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS
from gnomad.utils.constraint import (
    add_gencode_transcript_annotations,
    aggregate_expected_variants_expr,
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    apply_models,
    apply_plateau_models,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    calibration_model_group_expr,
    compute_pli,
    count_observed_and_possible_by_group,
    coverage_correction_expr,
    get_constraint_flags,
    oe_confidence_interval,
    single_variant_count_expr,
    single_variant_observed_and_possible_expr,
    weighted_agg_sum_expr,
)
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.vep import filter_vep_transcript_csqs_expr
from hail.utils.misc import divide_null, new_temp_file

from gnomad_constraint.resources.resource_utils import (
    AGGREGATE_SUM_FIELDS,
    CALIBRATION_GROUPING,
    COVERAGE_CUTOFF,
    MU_GROUPING,
    MUTATION_TYPE_FIELDS,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


# TODO: For now I am leaving this here instead of moving to gnomad_methods because
#  there is another PR in gnomad_methods that might change the way this function is
#  implemented.
def filter_freq_for_constraint(
    freq_expr: hl.ArrayExpression,
    freq_meta_expr: List[Dict[str, str]],
    gen_ancs: Optional[List[str]] = None,
    downsamplings: Optional[List[int]] = None,
    downsampling_gen_ancs: Optional[List[str]] = None,
    gen_anc_label: str = "gen_anc",
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
    :param gen_ancs: Optional list of genetic ancestries to include in the frequency
        array. Default is None.
    :param downsamplings: Optional list of downsamplings to include in the frequency
        array. Default is None.
    :param downsampling_gen_ancs: Optional list of genetic ancestries to include
        downsamplings frequencies for. Default is None.
    :param gen_anc_label: Label for the genetic ancestry field in the frequency
        metadata. Default is "gen_anc".
    :return: Filtered frequency array and metadata.
    """
    freq_meta = hl.eval(freq_meta_expr)
    meta_keep = [{"group": "adj"}]

    if gen_ancs is not None:
        meta_keep += [{"group": "adj", gen_anc_label: pop} for pop in gen_ancs]

    if downsamplings is not None:
        downsampling_pops = ["global"] + (downsampling_gen_ancs or [])
        meta_keep += [
            {"group": "adj", gen_anc_label: pop, "downsampling": str(ds)}
            for pop in downsampling_pops
            for ds in downsamplings
        ]

    meta_keep = [m for m in meta_keep if m in freq_meta]
    freq_expr = hl.array([freq_expr[freq_meta.index(m)] for m in meta_keep])

    return freq_expr, meta_keep


def get_annotations_for_computing_mu(
    locus_expr: hl.expr.LocusExpression,
    genomes_filter_expr: hl.expr.StructExpression,
    genomes_freq_expr: hl.expr.ArrayExpression,
    genomes_freq_meta: List[Dict[str, str]],
    genomes_coverage_expr: hl.expr.Int32Expression,
    gerp_expr: hl.expr.Float64Expression,
    most_severe_consequence_expr: hl.expr.StringExpression,
    gen_ancs: Optional[List[str]] = None,
    downsampling_level: int = 1000,
    min_cov: int = 15,
    max_cov: int = 60,
    gerp_lower_cutoff: float = -3.9885,
    gerp_upper_cutoff: float = 2.6607,
    ac_cutoff: int = 5,
) -> Tuple[hl.expr.StructExpression, hl.expr.StructExpression]:
    """
    Get the annotations that are needed to compute the mutation rate.

    The function will return the following annotations:

        - genomes_freq: Frequency array for the genomes dataset, filtered to the
          requested genetic ancestries and downsampling level.
        - observed_variants: This annotation is an array, where each element
          corresponds to whether the variant is observed in the genomes dataset for the
          frequency group at the corresponding index in the `genomes_freq` array.
          Must PASS genome filters, have AC <= 'ac_cutoff' at the specified
          `downsampling_level`, and have a genome mean coverage >= `min_cov`
          and <= `max_cov`. The boolean value is stored as an integer (0 or 1).
        - possible_variants: Whether the variant is considered a possible variant in
          the genomes dataset. This includes variants not in the genome dataset (genome
          AF undefined), or also considered in the observed variant set. The boolean
          value is stored as an integer (0 or 1).

    The observed and possible variant annotations are set to missing if the variant
    does not meet the following criteria:

        - Is autosomal.
        - Has a most severe transcript consequence of: "intron_variant" or
          "intergenic_variant".
        - Is at a site with GERP > `gerp_lower_cutoff` and < `gerp_upper_cutoff`.

    The function also returns a struct of the mutation rate globals:

        - freq_meta: Frequency metadata for the genomes dataset, filtered to the
          requested genetic ancestries and downsampling level.
        - ac_cutoff: Allele count cutoff used for the mutation rate calculation.
        - min_cov: Minimum genome coverage used for the mutation rate calculation.
        - max_cov: Maximum genome coverage used for the mutation rate calculation.
        - gerp_lower_cutoff: Minimum GERP score used for the mutation rate calculation.
        - gerp_upper_cutoff: Maximum GERP score used for the mutation rate calculation.
        - downsampling_level: Downsampling level used for the mutation rate calculation.
        - downsampling_idx: Index of the downsampling level in the frequency metadata.
        - most_severe_consequence: List of most severe transcript consequences used for
          the mutation rate calculation.

    .. note::

        Values for `gerp_lower_cutoff` and `gerp_upper_cutoff` default to -3.9885 and
        2.6607, respectively. These values were precalculated on the GRCh37 context
        table and define the 5th and 95th percentiles.

    :param locus_expr: Locus expression.
    :param genomes_filter_expr: Filter expression for the genomes dataset.
    :param genomes_freq_expr: Frequency array for the genomes dataset.
    :param genomes_freq_meta: Frequency metadata for the genomes dataset.
    :param genomes_coverage_expr: Mean genome coverage expression.
    :param gerp_expr: GERP score expression.
    :param most_severe_consequence_expr: Most severe consequence expression.
    :param gen_ancs: List of genetic ancestries to filter the genome frequency array to.
        Default is None, which includes only the full genome dataset.
    :param downsampling_level: Downsampling level to use for the mutation rate
        calculation. Default is 1000.
    :param min_cov: Minimum genome coverage for variant to be included. Default is 15.
    :param max_cov: Maximum genome coverage for variant to be included. Default is 60.
    :param gerp_lower_cutoff: Minimum GERP score for variant to be included. Default
        is -3.9885.
    :param gerp_upper_cutoff: Maximum GERP score for variant to be included. Default
        is 2.6607.
    :param ac_cutoff: Allele count cutoff for variant to be included. Default is 5.
    :return: Tuple containing the observed and possible variant annotations and the
        globals.
    """
    genomes_freq_expr, genomes_freq_meta = filter_freq_for_constraint(
        genomes_freq_expr,
        genomes_freq_meta,
        gen_ancs=None,
        downsamplings=[downsampling_level],
        downsampling_gen_ancs=gen_ancs,
        gen_anc_label="pop",
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
        genomes_freq=genomes_freq_expr,
        **hl.or_missing(
            keep_expr, single_variant_observed_and_possible_expr(genomes_freq_expr)
        ),
    )
    obs_pos_globals = hl.struct(
        freq_meta=genomes_freq_meta,
        ac_cutoff=ac_cutoff,
        min_cov=min_cov,
        max_cov=max_cov,
        gerp_lower_cutoff=gerp_lower_cutoff,
        gerp_upper_cutoff=gerp_upper_cutoff,
        genetic_ancestry_groups=gen_ancs or hl.missing(hl.tarray(hl.tstr)),
        downsampling_level=downsampling_level,
        downsampling_idx=downsampling_idx,
        most_severe_consequence=["intron_variant", "intergenic_variant"],
    )

    return obs_pos_expr, obs_pos_globals


def get_exome_coverage_expr(
    ht: hl.Table,
    exome_coverage_metric: str = "AN_percent",
) -> hl.expr.Int32Expression:
    """
    Get the exome coverage expression based on the specified metric.

    The requested `exome_coverage_metric` is extracted from the exome coverage
    annotations in the input `ht`:

        - "median": the expression returned is "median_approx" if it exists in
          `ht.coverage.exomes`, otherwise "median".
        - "AN": the expression returned is the exomes allele number (`ht.AN.exomes`).
        - "AN_percent": the expression returned is the percent of samples with a
          non-missing genotype, which is the exomes allele number (`ht.AN.exomes`)
          divided by the total number of alleles in the exomes dataset (pulled from
          `ht.an_globals.exomes.strata_sample_count` * 2) multiplied by 100.

    :param ht: Input Table with exome coverage information.
    :param exome_coverage_metric: Metric to use for exome coverage. One of ["median",
        "AN", "AN_percent"]. Default is "AN_percent".
    :return: Exome coverage expression.
    """
    if exome_coverage_metric == "median":
        # Obtain field name for median exome coverage.
        exome_coverage_metric = (
            "median_approx" if "median_approx" in ht.coverage.exomes else "median"
        )
        cov_expr = ht.coverage.exomes[exome_coverage_metric]
    elif exome_coverage_metric == "AN":
        cov_expr = ht.AN.exomes
    elif exome_coverage_metric == "AN_percent":
        # Calculate total allele number from strata_sample_count and annotate
        # exomes_AN_percent (percent samples with AN).
        an_sample_count = ht.an_globals.exomes.strata_sample_count
        an_meta = ht.an_globals.exomes.strata_meta

        # Get total AN count taking into account XX and XY samples for X and Y non-PAR.
        xx_index = an_meta.index({"group": "adj", "sex": "XX"})
        xy_index = an_meta.index({"group": "adj", "sex": "XY"})
        xx_an_sample_count = an_sample_count[xx_index]
        xy_an_sample_count = an_sample_count[xy_index]
        an_count = (
            hl.case()
            .when(ht.locus.in_x_nonpar(), (xx_an_sample_count * 2) + xy_an_sample_count)
            .when(ht.locus.in_y_nonpar(), xy_an_sample_count)
            .default(an_sample_count[0] * 2)
        )

        cov_expr = hl.int((ht.AN.exomes / an_count) * 100)
    else:
        raise ValueError(
            f"Exome coverage metric must be one of ['median', 'AN', 'AN_percent'], not {exome_coverage_metric}"
        )

    logger.info("Setting 'exome_coverage' to %s", exome_coverage_metric)

    return cov_expr


def get_exomes_observed_and_possible(
    exomes_filter_expr: hl.expr.SetExpression,
    exomes_freq_expr: hl.expr.ArrayExpression,
    exomes_freq_meta: List[Dict[str, str]],
    exomes_coverage_expr: hl.expr.Int32Expression,
    gen_ancs: Optional[List[str]] = None,
    include_downsamplings: bool = False,
    max_af: float = 0.001,
) -> Tuple[hl.expr.StructExpression, hl.expr.StructExpression]:
    """
    Get the observed and possible variants for the exomes dataset.

    The function returns a struct indicating whether the variant should be included in
    the observed and possible variant counts for the exomes dataset. The struct includes
    the following fields:

        - observed_variants: This annotation is an array, where each element corresponds
          to whether the variant is observed in the exomes dataset for the frequency
          group at the corresponding index in the `exomes_freq` array and has an
          AF <= 0.001. The boolean value is stored as an integer (0 or 1).
        - possible_variants: Whether the variant is considered a possible variant in the
          exomes dataset. This includes variants not in the exome dataset (exome AF
          undefined), or also considered in the observed variant set. The boolean value
          is stored as an integer (0 or 1).

    The observed and possible variant annotations are set to missing if the exome
    coverage is undefined or the variant does not pass the exome filters.

    The function also returns a struct with the global parameters for the observed and
    possible variant annotations:

        - exomes_freq_meta: Frequency metadata for the exomes dataset.
        - genetic_ancestry_groups: List of genetic ancestry groups used for the
          observed and possible variant annotations.
        - downsamplings: List of downsamplings used for the observed and possible
          variant annotations.

    :param exomes_filter_expr: Filter expression for the exomes dataset.
    :param exomes_freq_expr: Frequency array for the exomes dataset.
    :param exomes_freq_meta: Frequency metadata for the exomes dataset.
    :param exomes_coverage_expr: Exome coverage expression.
    :param gen_ancs: List of genetic ancestries to filter the exome frequency array to.
        Default is None, which includes only the full exomes dataset.
    :param include_downsamplings: Whether to include downsamplings in the observed and
        possible variant annotations. Default is False.
    :param max_af: Maximum allele frequency to consider a variant as observed. Default
        is 0.001.
    :return: Tuple containing the observed and possible variant annotations and the
        globals.
    """
    # If downsamplings are requested and 'genetic_ancestry_groups' is not specified,
    # use the pared-down downsamplings list.
    downsamplings = [m["downsampling"] for m in exomes_freq_meta if "downsampling" in m]
    downsamplings = DOWNSAMPLINGS["v4"] if gen_ancs is None else downsamplings
    downsamplings = downsamplings if include_downsamplings else None
    logger.info("The following downsamplings will be used: %s", downsamplings)

    # Filter frequency array for computing the observed expression on all requested
    # populations and downsamplings.
    exomes_freq_expr, exomes_freq_meta = filter_freq_for_constraint(
        exomes_freq_expr,
        exomes_freq_meta,
        gen_ancs=gen_ancs,
        downsamplings=downsamplings,
        downsampling_gen_ancs=gen_ancs if downsamplings is not None else None,
    )

    # If the exome coverage is undefined or the variant does not pass the exome filters,
    # set the observed and possible variant annotations to missing. Otherwise, set the
    # observed and possible variant annotations based on the frequency array.
    obs_pos_expr = hl.struct(
        **hl.or_missing(
            hl.is_defined(exomes_coverage_expr)
            & hl.or_else(hl.len(exomes_filter_expr) == 0, True),
            single_variant_observed_and_possible_expr(exomes_freq_expr, max_af=max_af),
        )
    )
    obs_pos_globals = hl.struct(
        exomes_freq_meta=exomes_freq_meta,
        genetic_ancestry_groups=gen_ancs or hl.missing(hl.tarray(hl.tstr)),
        downsamplings=downsamplings or hl.missing(hl.tarray(hl.tstr)),
        max_af=max_af,
    )

    return obs_pos_expr, obs_pos_globals


def get_build_calibration_model_annotation(
    exomes_coverage_expr: hl.expr.Int32Expression,
    transcript_csq_expr: hl.expr.ArrayExpression,
    cpg_expr: hl.expr.BooleanExpression,
    genomic_region_expr: hl.expr.StringExpression,
    synonymous_transcript_filter_field: str = "mane_select",
    low_cov_cutoff: Optional[int] = None,
    high_cov_cutoff: int = COVERAGE_CUTOFF,
    upper_cov_cutoff: Optional[int] = None,
    skip_coverage_model: bool = False,
    additional_grouping_exprs: Optional[Dict[str, hl.expr.Expression]] = None,
) -> Tuple[hl.expr.StructExpression, hl.expr.StructExpression]:
    """
    Get the annotation and globals for building the calibration models.

    The build model grouping is set to missing if the variant is not a
    "synonymous_variant" in a canonical or MANE Select transcript (depending on
    `synonymous_transcript_filter_field`). Otherwise, it is a struct with the following
    fields detailed in `calibration_model_group_expr`.

    :param exomes_coverage_expr: Exome coverage expression.
    :param transcript_csq_expr: Transcript consequences expression.
    :param cpg_expr: CpG expression.
    :param genomic_region_expr: Genomic region expression.
    :param synonymous_transcript_filter_field: Field used to filter to variants with a
        transcript consequence of "synonymous_variant". Default is "mane_select".
    :param low_cov_cutoff: Low coverage cutoff for the build models step. Default is
        None.
    :param high_cov_cutoff: High coverage cutoff for the build models step. Default is
        COVERAGE_CUTOFF.
    :param upper_cov_cutoff: Upper coverage cutoff for the build models step. Default is
        None.
    :param skip_coverage_model: Whether the coverage model should be skipped during the
        build models step. Default is False.
    :param additional_grouping_exprs: Additional grouping expressions to include in the
        build model. Default is None.
    :return: Tuple containing the build model expression and the global parameters.
    """
    # Determine the canonical and mane_select parameters for
    # 'filter_vep_transcript_csqs_expr' based on 'synonymous_transcript_filter_field'.
    if synonymous_transcript_filter_field == "canonical":
        canonical, mane_select = True, False
    elif synonymous_transcript_filter_field == "mane_select":
        canonical, mane_select = False, True
    else:
        raise ValueError(
            "synonymous_transcript_filter_field must be either 'canonical' or "
            "'mane_select'"
        )

    # Filter the VEP transcript consequences to include only synonymous transcripts.
    syn_csq_expr = filter_vep_transcript_csqs_expr(
        transcript_csq_expr,
        synonymous=True,
        ensembl_only=True,
        canonical=canonical,
        mane_select=mane_select,
    )

    # Define whether the variant should be included in the high or low coverage model.
    build_expr = calibration_model_group_expr(
        exomes_coverage_expr,
        cpg_expr,
        low_cov_cutoff=0 if low_cov_cutoff is None else low_cov_cutoff,
        high_cov_cutoff=high_cov_cutoff,
        upper_cov_cutoff=upper_cov_cutoff,
        skip_coverage_model=skip_coverage_model,
        additional_grouping_exprs=additional_grouping_exprs,
        cpg_in_high_only=True,
    )

    return hl.or_missing(syn_csq_expr.length() > 0, build_expr)


# TODO: We don't really need this, I just found int helpful to look over the
#  chosen parameters.
def print_global_struct(t: Union[hl.Table, hl.Struct, hl.StructExpression]) -> None:
    """
    Print the global struct.

    :param t: Table with globals or globals struct to print.
    :return: None
    """
    if isinstance(t, hl.Table):
        t = t.globals
    if isinstance(t, hl.StructExpression):
        t = hl.eval(t)

    def _get_pretty_print_globals(global_struct: hl.Struct, level: int = 1) -> str:
        output = ""
        level_tab = "".join(["    "] * level)
        for k, v in global_struct.items():
            if isinstance(v, hl.Struct):
                v = f"\n{_get_pretty_print_globals(v, level + 1)}"

            output += f"{level_tab}{k}: {v}\n"

        return output

    logger.info(
        "\nThe following parameters were used: \n%s", _get_pretty_print_globals(t)
    )


def prepare_ht_for_constraint_calculations(
    ht: hl.Table,
    exome_coverage_metric: str = "median",
    gen_ancs: Optional[List[str]] = None,
    include_downsamplings: bool = False,
    mu_downsampling_level: int = 1000,
    calculate_mutation_rate_min_cov: int = 15,
    calculate_mutation_rate_max_cov: int = 60,
    calculate_mutation_rate_gerp_lower_cutoff: float = -3.9885,
    calculate_mutation_rate_gerp_upper_cutoff: float = 2.6607,
    calculate_mutation_rate_ac_cutoff: int = 5,
    max_af: float = 0.001,
    build_model_low_cov_cutoff: Optional[int] = None,
    build_model_high_cov_cutoff: int = COVERAGE_CUTOFF,
    build_model_upper_cov_cutoff: Optional[int] = None,
    apply_model_low_cov_cutoff: Optional[int] = None,
    apply_model_high_cov_cutoff: int = COVERAGE_CUTOFF,
    skip_coverage_model: bool = False,
    synonymous_transcript_filter_field: str = "mane_select",
    additional_grouping_exprs: Optional[Dict[str, hl.expr.Expression]] = None,
) -> hl.Table:
    """
    Prepare Table for constraint calculations.

    This function is a wrapper around the functions that generate the annotations
    required for the constraint calculations. Please see the following functions for
    more information on the annotations generated:

        - `get_annotations_for_computing_mu`
        - `get_exomes_observed_and_possible`
        - `get_build_calibration_model_annotation`
        - `get_apply_calibration_model_annotation`

    :param ht: Annotated context Table.
    :param exome_coverage_metric: Metric to use for exome coverage. One of ["median",
        "AN", "AN_percent"]. Default is "median".
    :param gen_ancs: List of genetic ancestries to filter the frequency arrays to.
        Default is None, which includes only the full dataset.
    :param include_downsamplings: Whether to include downsamplings in the observed and
        possible variant annotations. Default is False.
    :param mu_downsampling_level: Downsampling level to use for the mutation rate
        calculation. Default is 1000.
    :param calculate_mutation_rate_min_cov: Minimum genome coverage for variant to be
        included in the mutation rate calculation. Default is 15.
    :param calculate_mutation_rate_max_cov: Maximum genome coverage for variant to be
        included in the mutation rate calculation. Default is 60.
    :param calculate_mutation_rate_gerp_lower_cutoff: Minimum GERP score for variant to
        be included in the mutation rate calculation. Default is -3.9885.
    :param calculate_mutation_rate_gerp_upper_cutoff: Maximum GERP score for variant to
        be included in the mutation rate calculation. Default is 2.6607.
    :param calculate_mutation_rate_ac_cutoff: Allele count cutoff for variant to be
        included in the mutation rate calculation. Default is 5.
    :param max_af: Maximum allele frequency to consider a variant as observed. Default
        is 0.001.
    :param build_model_low_cov_cutoff: Low coverage cutoff for the build models step.
        Default is None.
    :param build_model_high_cov_cutoff: High coverage cutoff for the build models step.
        Default is COVERAGE_CUTOFF.
    :param build_model_upper_cov_cutoff: Upper coverage cutoff for the build models
        step. Default is None.
    :param apply_model_low_cov_cutoff: Low coverage cutoff for the apply models step.
        Default is COVERAGE_CUTOFF.
    :param apply_model_high_cov_cutoff: High coverage cutoff for the apply models step.
        Default is COVERAGE_CUTOFF.
    :param skip_coverage_model: Whether the coverage model should be skipped during the
        build and apply models steps. Default is False.
    :param synonymous_transcript_filter_field: Field used to filter to variants with a
        transcript consequence of "synonymous_variant". Default is "canonical".
    :param additional_grouping_exprs: Additional grouping expressions to include in the
        build model and apply model annotations. Default is None.
    :return: Table with the computed annotations.
    """
    apply_model_grouping_exprs = {"genomic_region": ht.genomic_region}
    build_model_grouping_exprs = {
        **apply_model_grouping_exprs,
        **(additional_grouping_exprs or {}),
    }

    # Get the annotations relevant for computing the mutation rate.
    compute_mu_expr, compute_mu_globals = get_annotations_for_computing_mu(
        ht.locus,
        ht.filters.genomes,
        ht.freq.genomes,
        ht.freq_globals.genomes.freq_meta,
        ht.coverage.genomes.mean,
        ht.gerp,
        ht.vep.most_severe_consequence,
        gen_ancs=gen_ancs,
        downsampling_level=mu_downsampling_level,
        min_cov=calculate_mutation_rate_min_cov,
        max_cov=calculate_mutation_rate_max_cov,
        gerp_lower_cutoff=calculate_mutation_rate_gerp_lower_cutoff,
        gerp_upper_cutoff=calculate_mutation_rate_gerp_upper_cutoff,
        ac_cutoff=calculate_mutation_rate_ac_cutoff,
    )

    # Get an observed and possible variant annotation for the exomes dataset.
    exomes_coverage_expr = get_exome_coverage_expr(ht, exome_coverage_metric)
    exomes_obs_pos_expr, exomes_obs_pos_globals = get_exomes_observed_and_possible(
        ht.filters.exomes,
        ht.freq.exomes,
        hl.eval(ht.freq_globals.exomes.freq_meta),
        exomes_coverage_expr,
        gen_ancs=gen_ancs,
        include_downsamplings=include_downsamplings,
        max_af=max_af,
    )

    # Get the annotations relevant for building the calibration models.
    build_expr = get_build_calibration_model_annotation(
        exomes_coverage_expr,
        ht.vep.transcript_consequences,
        ht.cpg,
        ht.genomic_region,
        synonymous_transcript_filter_field=synonymous_transcript_filter_field,
        low_cov_cutoff=build_model_low_cov_cutoff,
        high_cov_cutoff=build_model_high_cov_cutoff,
        upper_cov_cutoff=build_model_upper_cov_cutoff,
        skip_coverage_model=skip_coverage_model,
        additional_grouping_exprs=build_model_grouping_exprs,
    )

    # Get the annotations relevant for applying the calibration models.
    apply_expr = calibration_model_group_expr(
        exomes_coverage_expr,
        ht.cpg,
        low_cov_cutoff=apply_model_low_cov_cutoff,
        high_cov_cutoff=apply_model_high_cov_cutoff,
        skip_coverage_model=skip_coverage_model,
        additional_grouping_exprs=apply_model_grouping_exprs,
    )

    # Annotate the Table with the computed annotations, and select only the relevant
    # fields.
    ht = ht.annotate(
        exomes_coverage=exomes_coverage_expr,
        compute_mu=compute_mu_expr,
        calibrate_mu=hl.struct(
            **exomes_obs_pos_expr,
            build_model=hl.or_missing(hl.is_defined(ht.sfs_bin), build_expr),
            apply_model=hl.or_missing(hl.is_defined(ht.sfs_bin), apply_expr),
        ),
    )
    ht = ht.drop("freq")

    # Build a struct with the global parameters for building the calibration models.
    mis_int = hl.missing(hl.tint)
    handle_none = lambda x: x if x is not None else mis_int
    ht = ht.select_globals(
        calculate_mu_globals=compute_mu_globals,
        build_models_globals=hl.struct(
            synonymous_transcript_filter_field=synonymous_transcript_filter_field,
            low_cov_cutoff=handle_none(build_model_low_cov_cutoff),
            high_cov_cutoff=build_model_high_cov_cutoff,
            upper_cov_cutoff=handle_none(build_model_upper_cov_cutoff),
            skip_coverage_model=skip_coverage_model,
            additional_model_grouping=list(build_model_grouping_exprs.keys()),
        ),
        apply_models_globals=hl.struct(
            low_cov_cutoff=handle_none(apply_model_low_cov_cutoff),
            high_cov_cutoff=apply_model_high_cov_cutoff,
            skip_coverage_model=skip_coverage_model,
            additional_model_grouping=list(apply_model_grouping_exprs.keys()),
        ),
        **exomes_obs_pos_globals,
    )

    print_global_struct(ht)

    return ht


def create_training_set(
    ht: hl.Table,
    mutation_ht: hl.Table,
    partition_hint=100,
) -> hl.Table:
    """
    Create the training set for the constraint model.

    The input `ht` should be prepared using `prepare_ht_for_constraint_calculations`.
    The `ht` is filtered to include only the rows that have a build model annotation.
    The observed and possible variants are counted by group and annotated with the
    mutation rate. The Table is then checkpointed to avoid memory and shuffle issues.

    :param ht: Table prepared using `prepare_ht_for_constraint_calculations`.
    :param mutation_ht: Mutation rate Table.
    :param partition_hint: Partition hint for the Table. Default is 100.
    :return: Training set Table.
    """
    # Selecting the only fields that are needed for the training set and filtering out
    # the rows that are not needed, then checkpointing the Table. This is added to
    # help avoid memory and shuffle issues.
    ht = ht.transmute(**ht.calibrate_mu)
    # TODO: From Konrad's script parser.add_argument('--skip_af_filter_upfront',
    #  help='Skip AF filter up front (to be applied later to ensure that it is not
    #  affecting population-specific constraint): not generally recommended',
    #  action='store_true')
    ht = ht.filter(hl.is_defined(ht.build_model) & (ht.possible_variants > 0))
    select_fields = {*MU_GROUPING, *MUTATION_TYPE_FIELDS, *CALIBRATION_GROUPING}
    ht = ht.select(
        *select_fields,
        "observed_variants",
        "possible_variants",
    )
    ht = ht.checkpoint(new_temp_file("create_training_set", "ht"))

    # Aggregate and count the observed and possible variants by group.
    ht = count_observed_and_possible_by_group(
        ht,
        ht.possible_variants,
        ht.observed_variants,
        additional_grouping=("methylation_level",)
        + MUTATION_TYPE_FIELDS
        + CALIBRATION_GROUPING,
        partition_hint=partition_hint,
    )

    # Annotate with mutation rate.
    ht = annotate_with_mu(ht, mutation_ht)

    return ht


def create_per_variant_expected_ht(
    ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    log10_coverage: bool = True,
    filter_to_apply_variants: bool = True,
    custom_vep_annotation: str = "transcript_consequences",
    use_mane_select: bool = False,
) -> hl.Table:
    """
    Create the per-variant expected Table.

    The input `ht` should be prepared using `prepare_ht_for_constraint_calculations`.
    The `ht` is filtered to include only the rows that have an apply model annotation (
    if `filter_to_apply_variants` is True). The Table is then annotated with the
    expected number of variants using `apply_models`. See the function `apply_models`
    for more information on the expected annotations.

    :param ht: Table prepared using `prepare_ht_for_constraint_calculations`.
    :param mutation_ht: Mutation rate Table.
    :param plateau_models: Plateau models for the constraint calculations.
    :param coverage_model: Coverage model for the constraint calculations.
    :param log10_coverage: Whether to use log10 coverage. Default is True.
    :param filter_to_apply_variants: Whether to filter to only the rows with an apply
        model annotation. Default is True.
    :param custom_vep_annotation: Custom VEP annotation to use. Default is
    :param use_mane_select: Whether to include MANE Select as a group. Default is False.
    :return: Per-variant expected Table.
    """
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

    calibrate_mu_fields = set(ht.calibrate_mu.keys())
    ht = ht.annotate(**ht.calibrate_mu)

    if filter_to_apply_variants:
        # TODO: From Konrad's script parser.add_argument('--skip_af_filter_upfront',
        #  help='Skip AF filter up front (to be applied later to ensure that it is not
        #  affecting population-specific constraint): not generally recommended',
        #  action='store_true')
        ht = ht.filter(hl.is_defined(ht.apply_model) & (ht.possible_variants > 0))

    ht = annotate_with_mu(ht, mutation_ht)

    sfs_bins = range(8)
    ht = ht.annotate(
        expected_variants_by_sfs_bin=[
            apply_models(
                ht.mu_snp,
                plateau_models[ht.apply_model.model_group.annotate(sfs_bin=b)],
                ht.possible_variants,
                coverage_model=coverage_model,
                coverage_expr=ht.exomes_coverage,
                model_group_expr=ht.apply_model,
                log10_coverage=log10_coverage,
            ).annotate(sfs_bin=b)
            for b in sfs_bins
        ]
    )

    ht, groupings = annotate_exploded_vep_for_constraint_groupings(
        ht=ht,
        vep_annotation=vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )
    ht = ht.checkpoint(new_temp_file(prefix="constraint", extension="ht"))

    am_ht = hl.read_table(
        "gs://gnomad/v4.1/constraint/resources/alpha_missense.esm.filters.ht"
    )
    new_loftee = hl.read_table(
        "gs://gnomad/v4.1/constraint/resources/split10_gnomAD_LoF_Ppost_misannot.filters.ht"
    )
    am_keyed = am_ht[ht.locus, ht.alleles, ht.transcript]
    new_loftee_keyed = new_loftee[ht.locus, ht.alleles, ht.transcript]
    ann_expr = {
        "am_0_999": hl.or_else(am_keyed.alpha_missense.am_0_999, False),
        **{
            f"am_per_{p}": hl.or_else(am_keyed.alpha_missense[f"am_per_{p}"], False)
            for p in [90, 95, 98, 99]
        },
        **{
            f"am_tx_per_{p}": hl.or_else(
                am_keyed.alpha_missense[f"am_tx_per_{p}"], False
            )
            for p in [90, 95, 98, 99]
        },
        **{
            f"esm_per_{p}": hl.or_else(am_keyed.esm[f"esm_per_{p}"], False)
            for p in [90, 95, 98, 99]
        },
        **{
            f"esm_tx_per_{p}": hl.or_else(am_keyed.esm[f"esm_tx_per_{p}"], False)
            for p in [90, 95, 98, 99]
        },
        **{
            f"new_loftee_{int(p * 100)}": hl.or_else(
                new_loftee_keyed[f"misannot_Pposterior_{p}"], False
            )
            for p in [0.2, 0.5, 0.8]
        },
    }
    ht = ht.annotate(**ann_expr)
    ht = ht.explode("expected_variants_by_sfs_bin")
    zero_array = hl.zeros(ht.exomes_freq_meta.length())
    ht = ht.annotate(
        **hl.if_else(
            ht.expected_variants_by_sfs_bin.sfs_bin == ht.sfs_bin,
            {
                "possible_variants": ht.possible_variants,
                "observed_variants": ht.observed_variants,
            },
            {
                "possible_variants": 0,
                "observed_variants": zero_array,
            },
        ),
        **ht.expected_variants_by_sfs_bin,
    )
    ht = ht.annotate_globals(
        apply_models_globals=ht.apply_models_globals.annotate(
            plateau_models=plateau_models,
            coverage_model=coverage_model,
            log10_coverage=log10_coverage,
            groupings=groupings + tuple(ann_expr.keys()) + ("sfs_bin",),
        )
    )

    # tmp_path = new_temp_file(prefix="constraint", extension="ht")
    # ht.drop(*calibrate_mu_fields).write(tmp_path)

    # return hl.read_table(tmp_path, _n_partitions=2000)

    return ht


def aggregate_per_variant_expected_ht(
    ht,
    include_mu_annotations_in_grouping: bool = False,
):
    """
    Aggregate the per-variant expected Table.

    The input `ht` should be the Table returned by `create_per_variant_expected_ht`.
    The Table is exploded by the VEP annotation and aggregated by "genomic_region",
    "context", "ref", "alt", "methylation_level", groupings returned by
    `annotate_exploded_vep_for_constraint_groupings` and fields in
    `additional_grouping` to get the observed and expected counts.

    :param ht: Table returned by `create_per_variant_expected_ht`.
    :param include_mu_annotations_in_grouping: Whether to include the mutation rate
        key annotations in the grouping. Default is False.
    :return: Table with the observed and expected counts.
    """
    groupings = [
        *(MU_GROUPING if include_mu_annotations_in_grouping else []),
        *[
            g
            for g in hl.eval(ht.apply_models_globals.groupings)
            if g not in MU_GROUPING
        ],
    ]
    ht = ht.group_by(*groupings).aggregate(**aggregate_expected_variants_expr(ht))

    return ht.naive_coalesce(1000)


# TODO: Move this up after review in this location.
def calculate_mu_by_downsampling(
    ht: hl.Table,
    additional_grouping: Tuple[str] = ("methylation_level",),
    total_mu: float = 1.2e-08,
) -> hl.Table:
    """
    Calculate mutation rate.

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context.
        - ref - the reference allele.
        - alt - the alternate base.
        - methylation_level - methylation_level.
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`.
        - mu_snp - SNP mutation rate.
        - annotations added by `annotate_mutation_type`.

    :param ht: Table returned by `prepare_ht_for_constraint_calculations`.
    :param additional_grouping: Annotations other than 'context', 'ref', and 'alt'.
        Default is ('methylation_level',).
    :param total_mu: The per-generation mutation rate. Default is 1.2e-08.
    :return: Mutation rate Table.
    """
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by context, ref, alt, and 'additional_grouping'.
    ht = count_observed_and_possible_by_group(
        ht,
        ht.compute_mu.possible_variants,
        ht.compute_mu.observed_variants,
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


# TODO: I think we decided this isn't needed right? We can just use canonical.
def filter_to_mane_select_over_canonical(ht: hl.Table) -> hl.Table:
    """
    Filter to MANE Select over canonical transcripts.

    Filter to only ensembl transcripts of the specified transcript filter. If MANE
    select is specified, and a gene does not have a MANE select transcript, use
    canonical instead.

    :param ht: Table with the MANE Select and canonical annotations.
    :return: Table filtered to MANE Select over canonical transcripts.
    """
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
    ms_ht = ms_ht.filter(
        (ms_ht.transcript.startswith("ENST"))
        & (
            (ms_ht._mane_present & ms_ht.mane_select)
            | (ms_ht._only_canonical & ms_ht.canonical)
        )
    )

    return ms_ht


# TODO: Move to gnomad_methods?
def add_oe_upper_rank_and_decile(
    ht: hl.Table,
    len_meta: int,
    use_mane_select_over_canonical: bool = True,
) -> hl.Table:
    """
    Compute the rank and decile of the oe upper confidence interval.

    :param ht: Table with the oe upper confidence interval.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical
        transcripts for filtering the Table when determining ranks for the lof oe upper
        confidence interval. If a gene
        does not have a MANE Select transcript, the canonical transcript (if available)
        will be used instead. Default is True.
    :return: Struct containing the rank and decile of the oe upper confidence interval.
    """
    total_count = ht.count()

    if use_mane_select_over_canonical:
        ms_ht = filter_to_mane_select_over_canonical(ht)
    else:
        ms_ht = ht.filter((ht.canonical) & (ht.transcript.startswith("ENST")))

    ms_ht = ms_ht.checkpoint(new_temp_file("constraint_metrics.canonical"))

    n_transcripts = ms_ht.count()
    logger.info(
        "Retaining %d out of %d transcripts to use for rank annotations.",
        n_transcripts,
        total_count,
    )

    ms_ht = ms_ht.annotate(upper_rank=hl.empty_array(hl.tint64))
    for i in range(len_meta):
        # Rank in ascending order.
        ms_ht = ms_ht.order_by(ms_ht.constraint_groups[i].oe_info[0].oe_ci.upper)
        ms_ht = ms_ht.add_index(name="rank")
        ms_ht = ms_ht.annotate(upper_rank=ms_ht.upper_rank.append(ms_ht.rank))

    ms_ht = ms_ht.annotate(
        upper_bin_sextile=ms_ht.upper_rank.map(lambda x: hl.int(x * 6 / n_transcripts)),
        upper_bin_decile=ms_ht.upper_rank.map(lambda x: hl.int(x * 10 / n_transcripts)),
    )

    # Map rank and bin annotations back to original Table.
    ms_ht = ms_ht.key_by(*list(ht.key))
    ms_keyed = ms_ht[ht.key]

    return ht.annotate(
        constraint_groups=hl.enumerate(ht.constraint_groups).map(
            lambda x: x[1].annotate(
                **{
                    k: ms_keyed[k][x[0]]
                    for k in ["upper_rank", "upper_bin_sextile", "upper_bin_decile"]
                }
            )
        )
    )


# TODO: Move to gnomad_methods?
def build_constraint_consequence_groups(
    csq_expr: hl.expr.ArrayExpression,
    lof_modifier_expr: hl.expr.StringExpression,
    classic_lof_annotations: Tuple = (
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ),
    additional_groupings: Dict[str, Dict[str, hl.expr.BooleanExpression]] = None,
    additional_grouping_combinations: List[List[str]] = None,
) -> Tuple[List[hl.expr.BooleanExpression], List[Dict[str, str]]]:
    """
    Build constraint consequence groups.

    The function builds constraint groups based on the consequence expression and LoF
    modifier expression. By default, the following groups are built:

        - csq_set: synonymous_variant, missense_variant
        - lof: classic, hc_lc, classic_hc_lc, hc

    The resulting meta and cooresponding constraint group filters are:

        - {"csq_set": "syn"}: synonymous_variant
        - {"csq_set": "mis"}: missense_variant
        - {"lof": "classic"}: classic LoF annotations
        - {"lof": "hc_lc"}: LoFTEE HC or LC
        - {"lof": "classic_hc_lc"}: classic LoF annotations with LoFTEE HC or LC
        - {"lof": "hc"}: LoF annotations with LoFTEE HC

    Additional groupings can be added to the constraint groups by specifying the
    `additional_groupings` parameter, and grouping combinations can also be added
    by specifying the `additional_grouping_combinations` parameter.

    :param csq_expr: Consequence expression.
    :param lof_modifier_expr: LoF modifier expression.
    :param classic_lof_annotations: Classic LoF Annotations used to filter the input
        Table. Default is {"stop_gained", "splice_donor_variant",
        "splice_acceptor_variant"}.
    :param additional_groupings: Additional groupings to add to the constraint groups.
        Default is None.
    :param additional_grouping_combinations: Additional grouping combinations to add to
        the constraint groups. Default is None.
    :return: Tuple containing the constraint group filters and the meta.
    """
    lof_classic_expr = hl.literal(set(classic_lof_annotations)).contains(csq_expr)
    lof_hc_expr = lof_modifier_expr == "HC"
    lof_hc_lc_expr = lof_hc_expr | (lof_modifier_expr == "LC")
    mis_expr = csq_expr == "missense_variant"
    annotation_dict = {
        "csq_set": {"syn": csq_expr == "synonymous_variant", "mis": mis_expr},
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

    annotation_dict.update(additional_groupings or {})
    additional_grouping_combinations = additional_grouping_combinations or []

    grouping_combinations = [["csq_set"], ["lof"]]
    grouping_combinations.extend(additional_grouping_combinations)

    meta = generate_filter_combinations(
        grouping_combinations,
        {k: list(v.keys()) for k, v in annotation_dict.items()},
    )
    constraint_group_filters = [
        functools.reduce(operator.ior, [annotation_dict[k][v] for k, v in m.items()])
        for m in meta
    ]

    return constraint_group_filters, meta


# TODO: Move to gnomad_methods?
def convert_multi_array_to_array_of_structs(
    t: Union[hl.Table, hl.expr.StructExpression],
    array_fields_to_combine: List[str],
    new_array_field: str,
) -> hl.Table:
    """
    Convert multiple arrays to an array of structs.

    :param t: Table or Struct to convert.
    :param array_fields_to_combine: Array fields to combine.
    :param new_array_field: Name of the new array field.
    :return: Table with the array fields combined into an array of structs named
        `new_array_field`.
    """
    logger.warning("This function assumes that all arrays have the same length!")
    return t.annotate(
        **{
            new_array_field: hl.range(t[array_fields_to_combine[0]].length()).map(
                lambda i: hl.struct(**{f: t[f][i] for f in array_fields_to_combine})
            )
        }
    ).drop(*array_fields_to_combine)


def aggregate_by_constraint_groups(
    ht: hl.Table,
    keys: Tuple = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Tuple = (
        "stop_gained",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ),
    additional_groupings: Dict[str, Dict[str, hl.expr.BooleanExpression]] = None,
    additional_grouping_combinations: List[List[str]] = None,
) -> hl.Table:
    """
    Aggregate the observed and expected variant info for synonymous variants, missense variants, and predicted loss-of-function (pLoF) variants.

    .. note::

        The following annotations should be present in `ht`:

            - modifier
            - annotation
            - observed_variants
            - mu
            - possible_variants
            - expected_variants

    :param ht: Input Table with the number of expected variants (output of
        `get_proportion_observed()`).
    :param keys: The keys of the output Table, defaults to ('gene', 'transcript',
        'canonical').
    :param classic_lof_annotations: Classic LoF Annotations used to filter the input
        Table. Default is {"stop_gained", "splice_donor_variant",
        "splice_acceptor_variant"}.
    :param additional_groupings: Additional groupings to add to the constraint groups.
        Default is None.
    :param additional_grouping_combinations: Additional grouping combinations to add to
        the constraint groups. Default is None.
    :return: Table with the aggregated observed and expected variant info for synonymous
        variants, missense variants, and pLoF variants.
    """
    # Build constraint groups.
    constraint_group_filters_expr, meta = build_constraint_consequence_groups(
        ht.annotation,
        ht.modifier,
        classic_lof_annotations=classic_lof_annotations,
        additional_groupings=additional_groupings,
        additional_grouping_combinations=additional_grouping_combinations,
    )
    ht = ht.annotate(constraint_groups=constraint_group_filters_expr)
    ht = ht.annotate_globals(constraint_group_meta=meta)
    ht = ht.checkpoint(
        new_temp_file("constraint_metrics.constraint_group_filters", "ht")
    )

    # Group by keys and get an aggregate sum of mu_snp, observed_variants,
    # possible_variants, predicted_proportion_observed, coverage_correction, and
    # expected_variants for each constraint group.
    ht = ht.group_by(*keys).aggregate(
        constraint_groups=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, aggregate_expected_variants_expr(ht)),
            ht.constraint_groups,
        )
    )

    # Add a 'no_variants' annotation indicating that there are zero observed variants
    # summed across pLoF, missense, and synonymous variants.
    ht = ht.annotate(
        no_variants=hl.sum(
            ht.constraint_groups.map(lambda x: hl.or_else(x.observed_variants[0], 0))
        )
        == 0
    )

    # Filter to only rows with at least 1 obs or exp across all keys in annotation_dict.
    ht = ht.filter(
        ~ht.no_variants
        & hl.any(
            ht.constraint_groups.map(
                lambda x: (hl.or_else(x.expected_variants[0], 0) > 0)
            )
        )
    )

    # Change format of arrays in constraint_groups to an array of structs.
    array_fields_to_combine = [
        k
        for k, v in ht.constraint_groups.dtype._element_type.items()
        if isinstance(v, hl.tarray)
    ]
    ht = ht.annotate(
        constraint_groups=ht.constraint_groups.map(
            lambda x: convert_multi_array_to_array_of_structs(
                x, array_fields_to_combine, "oe_info"
            )
        )
    )

    ht = ht.annotate_globals(constraint_group_meta=meta)

    return ht


def compute_constraint_metrics(
    ht: hl.Table,
    gencode_ht: hl.Table,
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
    raw_z_outlier_threshold_lower_lof: float = -8.0,
    raw_z_outlier_threshold_lower_missense: float = -8.0,
    raw_z_outlier_threshold_lower_syn: float = -8.0,
    raw_z_outlier_threshold_upper_syn: float = 8.0,
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
    :param expected_values: Dictionary containing the expected values for 'Null',
        'Rec', and 'LI' to use as starting values.
    :param min_diff_convergence: Minimum iteration change in LI to consider the EM
        model convergence criteria as met. Default is 0.001.
    :param raw_z_outlier_threshold_lower_lof: Value at which the raw z-score is considered an outlier for lof variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_lower_missense: Value at which the raw z-score is considered an outlier for missense variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_lower_syn: Lower value at which the raw z-score is considered an outlier for synonymous variants. Values below this threshold will be considered outliers. Default is -8.0.
    :param raw_z_outlier_threshold_upper_syn: Upper value at which the raw z-score is considered an outlier for synonymous variants. Values above this threshold will be considered outliers. Default is  8.0.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical transcripts for filtering the Table when determining ranks for the lof oe upper confidence interval.
        If a gene does not have a MANE Select transcript, the canonical transcript (if available) will be used instead. Default is True.
    :param gencode_ht: Table containing GENCODE annotations.
    :return: Table with pLI scores, observed:expected ratio, confidence interval of the
        observed:expected ratio, and z scores.
    """

    def _add_oe_ci_z(
        oe_info: hl.expr.StructExpression,
        m: Dict[str, str],
    ) -> hl.expr.StructExpression:
        """
        Add oe, oe_ci, and z_raw to the oe_info struct.

        :param oe_info: Struct containing the observed and expected variants.
        :return: Struct containing oe, oe_ci, and z_raw.
        """
        obs = oe_info.observed_variants
        exp = oe_info.expected_variants
        z_raw = calculate_raw_z_score(obs, exp)
        z_threshold = dict(
            {
                "lof": (raw_z_outlier_threshold_lower_lof, None),
                "mis": (raw_z_outlier_threshold_lower_missense, None),
                "syn": (
                    raw_z_outlier_threshold_lower_syn,
                    raw_z_outlier_threshold_upper_syn,
                ),
            }
        )
        z_threshold = z_threshold.get(
            hl.coalesce(m.get("lof"), m.get("csq_set", "None")),
            (None, None),
        )
        flags = get_constraint_flags(exp, z_raw, z_threshold[0], z_threshold[1])
        return oe_info.annotate(
            oe=divide_null(obs, exp),
            oe_ci=oe_confidence_interval(obs, exp),
            z_raw=z_raw,
            flags=add_filters_expr(filters=flags),
        )

    # Annotate with the observed:expected ratio, 95% confidence interval around the
    # observed:expected ratio, and z scores for each constraint group.
    meta = hl.eval(ht.constraint_group_meta)
    ht = ht.annotate(
        constraint_groups=hl.map(
            lambda x, m: x.annotate(
                oe_info=x.oe_info.map(lambda oe: _add_oe_ci_z(oe, m))
            ),
            ht.constraint_groups,
            meta,
        )
    )
    # ht = ht.annotate(constraint_flags=...)
    ht = ht.checkpoint(new_temp_file("constraint_metrics.oe.oe_ci.z_raw", "ht"))

    # Add z-score 'sd' annotation to globals.
    # ht = ht.annotate_globals(
    #    sd_raw_z=ht.aggregate(
    #        hl.agg.filter(
    #            ~ht.no_variants,
    #            [
    #                calculate_raw_z_score_sd(
    #                    ht.constraint_groups[i].oe_info[0].z_raw,
    #                    ht.constraint_groups[i].oe_info[0].flags,
    #                    mirror_neg_raw_z=m.get("csq_set") != "syn",
    #                )
    #                for i, m in enumerate(meta)
    #            ],
    #        )
    #    )
    # )

    # Compute z-score from raw z-score and standard deviations.
    # TODO: Need to fix z_score
    # ht = ht.annotate(
    #    constraint_groups=hl.map(
    #        lambda x, sd_raw_z: x.annotate(z_score=x.oe_info[0].z_raw / sd_raw_z),
    #        ht.constraint_groups,
    #        ht.sd_raw_z,
    #    )
    # )

    # Add a rank and decile of the upper confidence interval for MANE Select or
    # canonical ensembl transcripts.
    ht = add_oe_upper_rank_and_decile(ht, len(meta), use_mane_select_over_canonical)

    # TODO: Add back pLI computation
    # Compute the observed:expected ratio.
    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}

    # Add transcript annotations from GENCODE.
    ht = add_gencode_transcript_annotations(ht, gencode_ht)

    return ht


# TODO: Move to gnomad_methods?
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
