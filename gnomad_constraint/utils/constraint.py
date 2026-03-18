"""Script containing utility functions used in the constraint pipeline."""

import logging
from typing import Dict, List, Optional, Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS
from gnomad.utils.constraint import (
    add_gencode_transcript_annotations,
    aggregate_constraint_metrics_expr,
    annotate_bins_by_threshold,
    annotate_exploded_vep_for_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    apply_models,
    assemble_constraint_context_ht,
    build_constraint_consequence_groups,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    calibration_model_group_expr,
    compute_percentile_thresholds,
    compute_pli,
    count_observed_and_possible_by_group,
    get_constraint_flags,
    oe_confidence_interval,
    rank_array_element_metrics,
    variant_observed_and_possible_expr,
)
from gnomad.utils.file_utils import (
    convert_multi_array_to_array_of_structs,
    print_global_struct,
)
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.vep import (
    CSQ_CODING,
    filter_vep_transcript_csqs_expr,
    mane_select_over_canonical_filter_expr,
    update_loftee_end_trunc_filter,
)
from hail.utils.misc import divide_null, new_temp_file

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    AGGREGATE_SUM_FIELDS,
    CALIBRATION_GROUPING,
    CLASSIC_LOF_ANNOTATIONS,
    CONSTRAINT_GRANULARITIES,
    COVERAGE_CUTOFF,
    GENCODE_FIELD_RENAMES,
    MU_GROUPING,
    MUTATION_TYPE_FIELDS,
    PLI_EXPECTED_VALUES,
    RELEASE_CG_RENAME,
    RELEASE_CG_SELECT,
    RELEASE_CI_FIELDS,
    RELEASE_CI_FIELDS_WITH_RANK,
    RELEASE_GROUP_NAMES,
    RELEASE_GROUP_RENAMES,
    RELEASE_GROUPS_WITH_PLI,
    RELEASE_GROUPS_WITH_RANK,
    RELEASE_KEY_ORDER,
    RELEASE_LOF_FIELDS,
    RELEASE_PIPELINE_PARAM_GLOBALS,
    RELEASE_TOP_LEVEL_ANNOTATIONS,
    SFS_BIN_CUTOFFS,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def prepare_context_ht(
    ht: hl.Table,
    coverage_hts: Dict[str, hl.Table],
    an_hts: Dict[str, hl.Table],
    freq_hts: Dict[str, hl.Table],
    filter_hts: Dict[str, hl.Table],
    methylation_ht: hl.Table,
    gerp_ht: hl.Table,
    adj_r_ht: hl.Table,
    syn_adj_r_ht: hl.Table,
    sfs_bin_cutoffs: Tuple[float, ...] = SFS_BIN_CUTOFFS,
) -> hl.Table:
    """
    Annotate the context Table with coverage, AN, frequency, and constraint annotations.

    Applies the LOFTEE END_TRUNC filter fix, assembles the constraint context
    Table via :func:`assemble_constraint_context_ht`, then adds genomic region,
    SFS bin, adj_r, syn_adj_r, and coverage/AN reshaping annotations.

    :param ht: VEP context Table.
    :param coverage_hts: Dict mapping data type ("exomes", "genomes") to coverage
        Tables.
    :param an_hts: Dict mapping data type to allele number Tables.
    :param freq_hts: Dict mapping data type to frequency Tables (with ``freq``
        field).
    :param filter_hts: Dict mapping data type to filter Tables (with ``filters``
        field).
    :param methylation_ht: Methylation sites Table.
    :param gerp_ht: GERP scores Table.
    :param adj_r_ht: Table with adj_r annotation keyed by locus.
    :param syn_adj_r_ht: Table with synonymous DNM adj_r annotation keyed by locus.
    :param sfs_bin_cutoffs: Allele frequency upper bounds defining site frequency
        spectrum bins. Default is ``SFS_BIN_CUTOFFS``.
    :return: Annotated context Table.
    """
    # There was a bug in the GERP cutoffs used to filter transcripts with the
    # "END_TRUNC" filter in the LOFTEE VEP plugin resulting in some transcripts
    # being considered "HC" when they should have been "LC". We use the
    # `update_loftee_end_trunc_filter` function to correct this issue.
    ht = ht.annotate(
        vep=ht.vep.annotate(
            transcript_consequences=update_loftee_end_trunc_filter(
                ht.vep.transcript_consequences
            )
        )
    )
    ht = assemble_constraint_context_ht(
        ht,
        coverage_hts=coverage_hts,
        an_hts=an_hts,
        freq_hts=freq_hts,
        filter_hts=filter_hts,
        methylation_ht=methylation_ht,
        gerp_ht=gerp_ht,
        transformation_funcs=None,
    )

    # Add annotation for genomic region (autosome/PAR, X non-PAR, Y non-PAR).
    genomic_region_expr = (
        hl.case()
        .when(ht.locus.in_autosome_or_par(), "autosome_or_par")
        .when(ht.locus.in_x_nonpar(), "chrx_nonpar")
        .when(ht.locus.in_y_nonpar(), "chry_nonpar")
        .or_missing()
    )

    # Add annotation for SFS bin.
    af_expr = ht.freq.exomes[0].AF
    sfs_bin_expr = hl.case().when(hl.is_missing(af_expr), 0)
    for i, af in enumerate(sfs_bin_cutoffs):
        sfs_bin_expr = sfs_bin_expr.when(af_expr <= af, i)
    sfs_bin_expr = sfs_bin_expr.or_missing()

    return ht.annotate(
        coverage=hl.struct(
            exomes=ht.coverage.exomes.select("mean", "median_approx"),
            genomes=ht.coverage.genomes.select("mean", "median_approx"),
        ),
        AN=hl.struct(
            exomes=ht.AN.exomes[0],
            genomes=ht.AN.genomes[0],
        ),
        genomic_region=genomic_region_expr,
        adj_r=adj_r_ht[ht.locus].adj_r[ht.context],
        syn_adj_r=syn_adj_r_ht[ht.locus].adj_r[ht.context],
        sfs_bin=sfs_bin_expr,
    )


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
    the genetic ancestry groups in ``gen_ancs`` and the downsamplings in
    ``downsamplings``, for the genetic ancestry groups in
    ``downsampling_gen_ancs``.

    No matter the input, the frequency array is always filtered to include the
    "adj" frequency for the full dataset.

    If ``downsamplings`` is None, no downsamplings are included. If
    ``downsamplings`` is provided, and ``downsampling_gen_ancs`` is None, only
    the "global" downsampling is included. If ``downsampling_gen_ancs`` is
    provided, the downsamplings for the genetic ancestry groups in
    ``downsampling_gen_ancs`` are included as well as the "global"
    downsampling.

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
    meta_keep = [ADJ_FREQ_META]

    if gen_ancs is not None:
        meta_keep += [{**ADJ_FREQ_META, gen_anc_label: gen_anc} for gen_anc in gen_ancs]

    if downsamplings is not None:
        downsampling_gen_ancs = ["global"] + (downsampling_gen_ancs or [])
        meta_keep += [
            {**ADJ_FREQ_META, gen_anc_label: gen_anc, "downsampling": str(ds)}
            for gen_anc in downsampling_gen_ancs
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
          frequency group at the corresponding index in the ``genomes_freq`` array.
          Must PASS genome filters, have AC <= ``ac_cutoff`` at the specified
          ``downsampling_level``, and have a genome mean coverage >= ``min_cov``
          and <= ``max_cov``. The boolean value is stored as an integer (0 or 1).
        - possible_variants: Whether the variant is considered a possible variant in
          the genomes dataset. This includes variants not in the genome dataset (genome
          AF undefined), or also considered in the observed variant set. The boolean
          value is stored as an integer (0 or 1).

    The observed and possible variant annotations are set to missing if the variant
    does not meet the following criteria:

        - Is autosomal.
        - Has a most severe transcript consequence of: "intron_variant" or
          "intergenic_variant".
        - Is at a site with GERP > ``gerp_lower_cutoff`` and < ``gerp_upper_cutoff``.

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

        Values for ``gerp_lower_cutoff`` and ``gerp_upper_cutoff`` default to -3.9885 and
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
    # Always include the global downsampling; gen_ancs only controls which
    # per-ancestry downsamplings are included.
    genomes_freq_expr, genomes_freq_meta = filter_freq_for_constraint(
        genomes_freq_expr,
        genomes_freq_meta,
        gen_ancs=None,
        downsamplings=[downsampling_level],
        downsampling_gen_ancs=gen_ancs,
        gen_anc_label="pop",
    )
    downsampling_idx = genomes_freq_meta.index(
        {**ADJ_FREQ_META, "pop": "global", "downsampling": str(downsampling_level)}
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
    # 'intergenic_variant'.
    keep_expr &= hl.any(
        most_severe_consequence_expr == c
        for c in ["intron_variant", "intergenic_variant"]
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
            keep_expr, variant_observed_and_possible_expr(genomes_freq_expr)
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

    The requested ``exome_coverage_metric`` is extracted from the exome coverage
    annotations in the input ``ht``:

        - "median": the expression returned is "median_approx" if it exists in
          ``ht.coverage.exomes``, otherwise "median".
        - "AN": the expression returned is the exomes allele number (``ht.AN.exomes``).
        - "AN_percent": the expression returned is the percent of samples with a
          non-missing genotype, which is the exomes allele number (``ht.AN.exomes``)
          divided by the total number of alleles in the exomes dataset (pulled from
          ``ht.an_globals.exomes.strata_sample_count`` * 2) multiplied by 100.

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
        xx_index = an_meta.index({**ADJ_FREQ_META, "sex": "XX"})
        xy_index = an_meta.index({**ADJ_FREQ_META, "sex": "XY"})
        xx_an_sample_count = an_sample_count[xx_index]
        xy_an_sample_count = an_sample_count[xy_index]
        an_count = (
            hl.case()
            .when(ht.locus.in_x_nonpar(), (xx_an_sample_count * 2) + xy_an_sample_count)
            .when(ht.locus.in_y_nonpar(), xy_an_sample_count)
            .default(an_sample_count[0] * 2)  # Index 0 is adj (all samples).
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
          group at the corresponding index in the ``exomes_freq`` array and has an
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
    downsamplings = sorted(map(int, list(set(downsamplings))))
    downsamplings = downsamplings if include_downsamplings else None
    logger.info("The following downsamplings will be used: %s", downsamplings)

    # Filter frequency array for computing the observed expression on all requested
    # genetic ancestry groups and downsamplings.
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
    exomes_freq_expr = hl.or_missing(hl.len(exomes_filter_expr) == 0, exomes_freq_expr)
    obs_pos_expr = hl.struct(
        exomes_freq=exomes_freq_expr,
        **hl.or_missing(
            hl.is_defined(exomes_coverage_expr),
            variant_observed_and_possible_expr(exomes_freq_expr, max_af=max_af),
        ),
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
) -> hl.expr.StructExpression:
    """
    Get the annotation for building the calibration models.

    The build model grouping is set to missing if the variant is not a
    "synonymous_variant" in a canonical or MANE Select transcript (depending on
    ``synonymous_transcript_filter_field``). Otherwise, it is a struct with the
    following fields detailed in ``calibration_model_group_expr``.

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
    :return: Build model struct expression, or missing if no synonymous transcripts.
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
        additional_grouping_exprs={"genomic_region": genomic_region_expr},
        cpg_in_high_only=True,
    )

    return hl.or_missing(syn_csq_expr.length() > 0, build_expr)


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
) -> hl.Table:
    """
    Prepare Table for constraint calculations.

    This function is a wrapper around the functions that generate the annotations
    required for the constraint calculations. Please see the following functions for
    more information on the annotations generated:

        - ``get_annotations_for_computing_mu``
        - ``get_exomes_observed_and_possible``
        - ``get_build_calibration_model_annotation``
        - ``calibration_model_group_expr`` (for apply model annotations)

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
        Default is None.
    :param apply_model_high_cov_cutoff: High coverage cutoff for the apply models step.
        Default is COVERAGE_CUTOFF.
    :param skip_coverage_model: Whether the coverage model should be skipped during the
        build and apply models steps. Default is False.
    :param synonymous_transcript_filter_field: Field used to filter to variants with a
        transcript consequence of "synonymous_variant". Default is "mane_select".
    :return: Table with the computed annotations.
    """
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
    )

    # Get the annotations relevant for applying the calibration models.
    apply_expr = calibration_model_group_expr(
        exomes_coverage_expr,
        ht.cpg,
        low_cov_cutoff=apply_model_low_cov_cutoff,
        high_cov_cutoff=apply_model_high_cov_cutoff,
        skip_coverage_model=skip_coverage_model,
        additional_grouping_exprs={"genomic_region": ht.genomic_region},
    )

    # Annotate the Table with the computed annotations, and select only the relevant
    # fields.
    ht = ht.annotate(
        exomes_coverage=exomes_coverage_expr,
        compute_mu=compute_mu_expr,
        calibrate_mu=hl.struct(
            **exomes_obs_pos_expr, build_model=build_expr, apply_model=apply_expr
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
        ),
        apply_models_globals=hl.struct(
            low_cov_cutoff=handle_none(apply_model_low_cov_cutoff),
            high_cov_cutoff=apply_model_high_cov_cutoff,
            skip_coverage_model=skip_coverage_model,
        ),
        **exomes_obs_pos_globals,
    )

    print_global_struct(ht)

    return ht


def create_training_set(
    ht: hl.Table,
    mutation_ht: hl.Table,
    partition_hint: int = 100,
) -> hl.Table:
    """
    Create the training set for the constraint model.

    The input ``ht`` should be prepared using
    ``prepare_ht_for_constraint_calculations``. The ``ht`` is filtered to include only
    the rows that have a build model annotation. The observed and possible variants are
    counted by group and annotated with the mutation rate. The Table is then
    checkpointed to avoid memory and shuffle issues.

    :param ht: Table prepared using ``prepare_ht_for_constraint_calculations``.
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


def _prepare_ht_for_apply_models(
    ht: hl.Table,
    custom_vep_annotation: str = "transcript_consequences",
    use_mane_select: bool = False,
) -> Tuple[hl.Table, List[str]]:
    """
    Prepare a preprocessed Table for model application.

    Promotes ``calibrate_mu`` fields, filters to rows with a defined apply model
    annotation and positive possible variant count, and explodes VEP annotations
    to per-transcript rows.

    :param ht: Table prepared using ``prepare_ht_for_constraint_calculations``.
    :param custom_vep_annotation: Custom VEP annotation to use. Default is
        ``"transcript_consequences"``.
    :param use_mane_select: Whether to include MANE Select as a group. Default is
        False.
    :return: Tuple of (prepared Table, list of VEP grouping field names).
    """
    if custom_vep_annotation == "worst_csq_by_gene" and use_mane_select:
        raise ValueError(
            "'mane_select' cannot be set to True when custom_vep_annotation is set"
            " to 'worst_csq_by_gene'."
        )
    include_canonical_group = custom_vep_annotation != "worst_csq_by_gene"
    include_mane_select_group = include_canonical_group and use_mane_select

    ht = ht.annotate(**ht.calibrate_mu)
    ht = ht.filter(hl.is_defined(ht.apply_model) & (ht.possible_variants > 0))

    ht, groupings = annotate_exploded_vep_for_constraint_groupings(
        ht=ht,
        vep_annotation=custom_vep_annotation,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )

    return ht, groupings


def _apply_constraint_models(
    ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    log10_coverage: bool = True,
) -> hl.Table:
    """
    Apply plateau and coverage models to a Table with ``mu_snp`` and ``apply_model``.

    :param ht: Table with ``mu_snp``, ``possible_variants``, ``exomes_coverage``,
        and ``apply_model`` fields.
    :param plateau_models: Plateau models for the constraint calculations.
    :param coverage_model: Coverage model for the constraint calculations.
    :param log10_coverage: Whether to use log10 coverage. Default is True.
    :return: Table annotated with model outputs (``mu``,
        ``predicted_proportion_observed``, ``expected_variants``,
        ``coverage_correction``).
    """
    return ht.annotate(
        **apply_models(
            ht.mu_snp,
            plateau_models.get(ht.apply_model.model_group),
            ht.possible_variants,
            coverage_model=coverage_model,
            coverage_expr=ht.exomes_coverage,
            model_group_expr=ht.apply_model,
            log10_coverage=log10_coverage,
        )
    )


def _annotate_apply_models_globals(
    ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    log10_coverage: bool,
    groupings: List[str],
) -> hl.Table:
    """
    Annotate the Table with model parameters in ``apply_models_globals``.

    If the Table already has ``apply_models_globals``, the new fields are added
    to the existing struct. Otherwise a new struct is created.

    :param ht: Input Table.
    :param plateau_models: Plateau models used.
    :param coverage_model: Coverage model used.
    :param log10_coverage: Whether log10 coverage was used.
    :param groupings: List of grouping field names.
    :return: Table with updated ``apply_models_globals`` global.
    """
    model_params = hl.struct(
        plateau_models=plateau_models,
        coverage_model=coverage_model,
        log10_coverage=log10_coverage,
        groupings=groupings,
    )
    if "apply_models_globals" in ht.globals:
        ht = ht.annotate_globals(
            apply_models_globals=ht.apply_models_globals.annotate(**model_params)
        )
    else:
        ht = ht.annotate_globals(apply_models_globals=model_params)
    return ht


def create_per_variant_expected_ht(
    ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    log10_coverage: bool = True,
    custom_vep_annotation: str = "transcript_consequences",
    use_mane_select: bool = False,
) -> hl.Table:
    """
    Create the per-variant expected Table.

    The input ``ht`` should be prepared using
    ``prepare_ht_for_constraint_calculations``. The ``ht`` is filtered to include only
    the rows that have an apply model annotation. The Table is then annotated with the
    expected number of variants using ``apply_models``. See the function
    ``apply_models`` for more information on the expected annotations.

    :param ht: Table prepared using ``prepare_ht_for_constraint_calculations``.
    :param mutation_ht: Mutation rate Table.
    :param plateau_models: Plateau models for the constraint calculations.
    :param coverage_model: Coverage model for the constraint calculations.
    :param log10_coverage: Whether to use log10 coverage. Default is True.
    :param custom_vep_annotation: Custom VEP annotation to use. Default is
        ``"transcript_consequences"``.
    :param use_mane_select: Whether to include MANE Select as a group. Default is False.
    :return: Per-variant expected Table.
    """
    calibrate_mu_fields = set(ht.calibrate_mu.keys())

    ht, groupings = _prepare_ht_for_apply_models(
        ht, custom_vep_annotation, use_mane_select
    )

    ht = annotate_with_mu(ht, mutation_ht)
    ht = _apply_constraint_models(ht, plateau_models, coverage_model, log10_coverage)
    ht = _annotate_apply_models_globals(
        ht, plateau_models, coverage_model, log10_coverage, groupings
    )

    return ht.drop(*calibrate_mu_fields)


def aggregate_per_variant_expected_ht(
    ht: hl.Table,
    include_mu_annotations_in_grouping: bool = False,
) -> hl.Table:
    """
    Aggregate the per-variant expected Table.

    The input ``ht`` should be the Table returned by ``create_per_variant_expected_ht``.
    The Table is aggregated by the groupings stored in
    ``apply_models_globals.groupings`` to get the observed and expected counts.

    :param ht: Table returned by ``create_per_variant_expected_ht``.
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
    if "calibrate_mu" in ht.row:
        ht = ht.annotate(**ht.calibrate_mu)

    ht = ht.filter(hl.set(CSQ_CODING).contains(ht.annotation))
    ht = ht.key_by().select(*groupings, *AGGREGATE_SUM_FIELDS)
    ht = ht.checkpoint(new_temp_file("pre_aggregation", "ht"))

    ht = ht.group_by(*groupings).aggregate(**aggregate_constraint_metrics_expr(ht))
    ht = ht.checkpoint(new_temp_file("post_aggregation", "ht"))

    return ht.naive_coalesce(1000)


def create_aggregated_expected_ht(
    ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    log10_coverage: bool = True,
    custom_vep_annotation: str = "transcript_consequences",
    use_mane_select: bool = False,
    partition_hint: int = 100,
) -> hl.Table:
    """
    Create aggregated expected variant counts by first aggregating, then applying models.

    Unlike :func:`create_per_variant_expected_ht`, which applies models per-variant and
    then aggregates, this function first aggregates observed and possible variant counts
    by VEP groupings and coverage, then applies plateau and coverage models on the
    aggregated counts. The output is compatible with
    :func:`aggregate_by_constraint_groups`.

    The steps are:

        1. Explode VEP annotations to get per-transcript rows.
        2. Aggregate observed and possible counts by VEP groupings, coverage, and
           mutation rate context (using ``count_observed_and_possible_by_group``).
        3. Annotate with mutation rate and apply plateau/coverage models on the
           aggregated counts.
        4. Aggregate by VEP groupings only (summing model outputs across coverage
           and context groups).

    :param ht: Table prepared using ``prepare_ht_for_constraint_calculations``.
    :param mutation_ht: Mutation rate Table.
    :param plateau_models: Plateau models for the constraint calculations.
    :param coverage_model: Coverage model for the constraint calculations.
    :param log10_coverage: Whether to use log10 coverage. Default is True.
    :param custom_vep_annotation: Custom VEP annotation to use. Default is
        ``"transcript_consequences"``.
    :param use_mane_select: Whether to include MANE Select as a group. Default is
        False.
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :return: Table with aggregated expected variant counts, compatible with
        ``aggregate_by_constraint_groups``.
    """
    ht, groupings = _prepare_ht_for_apply_models(
        ht, custom_vep_annotation, use_mane_select
    )

    # Filter to coding consequences.
    ht = ht.filter(hl.set(CSQ_CODING).contains(ht.annotation))

    # Aggregate observed and possible counts by VEP groupings, coverage, and mutation
    # rate context. The additional_grouping includes VEP groupings (gene, transcript,
    # annotation, etc.) plus the apply_model fields needed for model application.
    vep_groupings = tuple(g for g in groupings if g not in MU_GROUPING)
    ht = count_observed_and_possible_by_group(
        ht,
        ht.possible_variants,
        ht.observed_variants,
        additional_grouping=vep_groupings
        + ("exomes_coverage", "apply_model")
        + tuple(f for f in MUTATION_TYPE_FIELDS if f not in MU_GROUPING),
        partition_hint=partition_hint,
    )

    # Annotate with mutation rate and apply models on the aggregated counts.
    ht = annotate_with_mu(ht, mutation_ht)
    ht = _apply_constraint_models(ht, plateau_models, coverage_model, log10_coverage)
    ht = ht.checkpoint(new_temp_file("aggregated_apply_models", "ht"))

    # Aggregate by VEP groupings only, summing model outputs across coverage and
    # context groups to produce the same schema as aggregate_per_variant_expected_ht.
    ht = ht.key_by().select(*vep_groupings, *AGGREGATE_SUM_FIELDS)
    ht = ht.group_by(*vep_groupings).aggregate(**aggregate_constraint_metrics_expr(ht))

    ht = _annotate_apply_models_globals(
        ht, plateau_models, coverage_model, log10_coverage, list(vep_groupings)
    )

    return ht


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
        - downsampling_counts_{gen_anc} - variant counts in downsamplings for genetic
          ancestry groups in ``gen_ancs``.
        - mu_snp - SNP mutation rate.
        - annotations added by ``annotate_mutation_type``.

    :param ht: Table returned by ``prepare_ht_for_constraint_calculations``.
    :param additional_grouping: Annotations other than "context", "ref", and "alt".
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


def get_transcript_filter_expr(
    ht: hl.Table,
    use_mane_select_over_canonical: bool = True,
    mane_select_only: bool = False,
) -> hl.expr.BooleanExpression:
    """
    Return a filter expression for selecting one representative transcript per gene.

    Operates on an exploded, transcript-keyed table (one row per gene/transcript
    pair) — not on VEP ``transcript_consequences`` arrays.

    :param ht: Table with ``transcript``, ``mane_select``, ``canonical``, and
        ``gene_id`` annotations.
    :param use_mane_select_over_canonical: When ``True`` (default), prefer MANE
        Select transcripts, falling back to canonical for genes without a MANE
        Select entry. When ``False``, use canonical transcripts only. Ignored when
        ``mane_select_only`` is ``True``.
    :param mane_select_only: When ``True``, restrict to ENST MANE Select transcripts
        only, with no canonical fallback. Default is ``False``.
    :return: Boolean expression that is ``True`` for the selected transcripts.
    """
    if mane_select_only:
        return ht.transcript.startswith("ENST") & ht.mane_select
    elif use_mane_select_over_canonical:
        return mane_select_over_canonical_filter_expr(
            ht.transcript, ht.mane_select, ht.canonical, ht.gene_id
        )
    else:
        return ht.transcript.startswith("ENST") & ht.canonical


def add_oe_upper_rank_and_bins(
    ht: hl.Table,
    use_mane_select_over_canonical: bool = True,
    mane_select_only: bool = False,
    bin_granularities: Optional[Dict[str, int]] = None,
) -> hl.Table:
    """
    Compute the rank and bins of the oe upper confidence interval.

    Thin wrapper around :func:`rank_array_element_metrics` that extracts the
    discretized Poisson and gamma upper CI values from each constraint group's
    first oe_info element.

    :param ht: Table with the oe upper confidence interval.
    :param use_mane_select_over_canonical: Use MANE Select over canonical transcripts
        for ranking, falling back to canonical when MANE Select is absent for a gene.
        Default is True. Ignored when ``mane_select_only`` is True.
    :param mane_select_only: Restrict ranking to ENST MANE Select transcripts only,
        with no canonical fallback. Default is False.
    :param bin_granularities: Mapping of bin name to multiplier used to assign each
        transcript to a bin (``hl.int(rank * multiplier / n_transcripts)``). Each entry
        produces a ``bin_{name}`` field. Default is
        ``{"percentile": 100, "decile": 10, "sextile": 6}``.
    :return: Input table with ``oe_ci_{ci}_rank`` fields added at the constraint-group
        level (e.g., ``oe_ci_discretized_poisson_rank``, ``oe_ci_gamma_rank``), each
        a struct with ``rank`` and ``bin_{name}`` fields for every entry in
        ``bin_granularities``. Transcripts excluded from ranking have these fields set
        to missing.
    """
    ci_fields = ["discretized_poisson", "gamma"]

    return rank_array_element_metrics(
        ht,
        array_field="constraint_groups",
        element_value_fn=lambda x: {
            f"oe_ci_{ci}": x.oe_info[0][f"oe_ci_{ci}"].upper for ci in ci_fields
        },
        filter_fn=lambda t: get_transcript_filter_expr(
            t, use_mane_select_over_canonical, mane_select_only
        ),
        bin_granularities=bin_granularities,
    )


def aggregate_by_constraint_groups(
    ht: hl.Table,
    keys: Tuple = ("gene", "transcript", "canonical"),
    classic_lof_annotations: Tuple = CLASSIC_LOF_ANNOTATIONS,
    additional_groupings: Optional[
        Dict[str, Dict[str, hl.expr.BooleanExpression]]
    ] = None,
    additional_grouping_combinations: Optional[List[List[str]]] = None,
) -> hl.Table:
    """
    Aggregate observed and expected variant info for synonymous, missense, and pLoF variants.

    .. note::

        The following annotations should be present in ``ht``:

            - modifier
            - annotation
            - observed_variants
            - mu
            - possible_variants
            - expected_variants

    :param ht: Input Table with observed and expected variant counts (output of the
        apply models step).
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
            lambda f: hl.agg.filter(f, aggregate_constraint_metrics_expr(ht)),
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
        | hl.any(
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

    return ht


def _compute_coverage_metrics(
    ht: hl.Table,
    gencode_cds_ht: hl.Table,
    an_coverage_threshold: int = 90,
) -> hl.Table:
    """Compute per-transcript proportion of CDS bases with adequate coverage.

    Uses the ``exomes_coverage`` field from the preprocessed context table
    (AN as a percentage of total alleles) to determine what fraction of CDS
    bases per transcript meet the coverage threshold.

    :param ht: Preprocessed context Hail Table with ``exomes_coverage`` (AN percent,
        0-100) per position.
    :param gencode_cds_ht: GENCODE CDS positions Table keyed by locus with
        ``transcript_id`` array.
    :param an_coverage_threshold: Minimum ``exomes_coverage`` value (0-100) for a
        position to be considered adequately covered. Default is 90.
    :return: Table keyed by ``transcript`` with ``prop_bp_AN90``.
    """
    # Deduplicate context table by locus (3 SNV alts per position share coverage).
    ht = ht.key_by("locus").select("exomes_coverage").distinct()

    # Join CDS positions with context coverage.
    ht = gencode_cds_ht.annotate(
        exomes_coverage=ht[gencode_cds_ht.locus].exomes_coverage
    )
    ht = ht.filter(hl.is_defined(ht.exomes_coverage))

    # Explode by transcript and aggregate.
    ht = ht.explode("transcript_id").cache()

    return ht.group_by(transcript=ht.transcript_id).aggregate(
        prop_bp_AN90=hl.agg.fraction(ht.exomes_coverage >= an_coverage_threshold),
    )


def _compute_site_quality_metrics(
    ht: hl.Table,
    gencode_cds_ht: hl.Table,
) -> hl.Table:
    """Compute per-transcript mapping quality and region flag metrics.

    Computes mean AS_MQ, proportion of sites in segmental duplications, and
    proportion of sites in low-complexity regions from variant sites within
    CDS regions.

    :param ht: gnomAD exomes sites Hail Table.
    :param gencode_cds_ht: GENCODE CDS positions Table keyed by locus with
        ``transcript_id`` array.
    :return: Table keyed by ``transcript`` with ``mean_AS_MQ``,
        ``prop_segdup``, and ``prop_LCR``.
    """
    # Deduplicate sites by locus.
    ht = ht.select("region_flags", AS_MQ=ht.info.AS_MQ)
    ht = ht.key_by("locus").select("AS_MQ", "region_flags").distinct()

    # Join with GENCODE CDS to get per-locus transcript IDs.
    ht = ht.annotate(transcript_id=gencode_cds_ht[ht.locus].transcript_id)
    ht = ht.filter(hl.is_defined(ht.transcript_id)).explode("transcript_id").cache()

    return ht.group_by(transcript=ht.transcript_id).aggregate(
        mean_AS_MQ=hl.agg.mean(ht.AS_MQ),
        prop_segdup=hl.agg.fraction(ht.region_flags.segdup),
        prop_LCR=hl.agg.fraction(ht.region_flags.lcr),
    )


def compute_gene_quality_metrics(
    context_ht: hl.Table,
    exomes_ht: hl.Table,
    gencode_cds_ht: hl.Table,
    an_coverage_threshold: int = 90,
) -> hl.Table:
    """Compute per-transcript gene quality metrics.

    Combines coverage metrics from :func:`_compute_coverage_metrics` and
    site quality metrics from :func:`_compute_site_quality_metrics` into a
    single Table with release-ready fields:

    - ``gene_quality_metrics``: struct with ``exome_prop_bp_AN90``,
      ``exome_mean_AS_MQ``, ``exome_prop_segdup``, ``exome_prop_LCR``.
    - ``gene_flags``: set of flag strings (``low_exome_mapping_quality``
      when mean AS_MQ < 50, ``low_exome_coverage`` when
      prop_bp_AN90 < 0.1).

    :param context_ht: Preprocessed context Hail Table with
        ``exomes_coverage`` (AN percent, 0-100) per position.
    :param exomes_ht: gnomAD exomes sites Hail Table.
    :param gencode_cds_ht: GENCODE CDS positions Table keyed by locus with
        ``transcript_id`` array (output of
        :func:`~gnomad_constraint.resources.resource_utils.get_gencode_cds_ht`).
    :param an_coverage_threshold: Minimum ``exomes_coverage`` value (0-100)
        for a position to be considered adequately covered. Default is 90.
    :return: Table keyed by ``transcript`` with ``gene_quality_metrics``
        and ``gene_flags``.
    """
    an90_ht = _compute_coverage_metrics(
        context_ht, gencode_cds_ht, an_coverage_threshold
    ).cache()
    sites_ht = _compute_site_quality_metrics(exomes_ht, gencode_cds_ht).cache()

    ht = an90_ht.annotate(**sites_ht[an90_ht.transcript])
    ht = ht.select(
        gene_quality_metrics=hl.struct(**{f"exome_{f}": ht[f] for f in ht.row_value}),
        gene_flags=add_filters_expr(
            {
                "low_exome_mapping_quality": ht.mean_AS_MQ < 50,
                "low_exome_coverage": ht.prop_bp_AN90 < 0.1,
            }
        ),
    )

    return ht.key_by("transcript")


def _annotate_oe_ci_z(
    ht: hl.Table,
    z_thresholds: Dict[str, Tuple[Optional[float], Optional[float]]],
) -> hl.Table:
    """
    Annotate constraint groups with OE ratio, confidence intervals, z-scores, and flags.

    For each constraint group's ``oe_info`` entries, adds:

        - ``oe`` — observed / expected ratio.
        - ``oe_ci_discretized_poisson`` — discretized Poisson CI.
        - ``oe_ci_gamma`` — gamma-distribution CI.
        - ``z_raw`` — raw z-score.

    Then adds per-group ``flags`` based on z-score outlier thresholds.

    :param ht: Table with ``constraint_groups`` array.
    :param z_thresholds: Mapping from constraint category (``"lof"``, ``"mis"``,
        ``"syn"``) to ``(lower, upper)`` raw z-score outlier thresholds.
    :return: Table with OE, CI, z-score, and flag annotations.
    """
    ht = ht.annotate(
        constraint_groups=ht.constraint_groups.map(
            lambda x: x.annotate(
                oe_info=x.oe_info.map(
                    lambda oe_info: oe_info.annotate(
                        oe=divide_null(
                            oe_info.observed_variants, oe_info.expected_variants
                        ),
                        oe_ci_discretized_poisson=oe_confidence_interval(
                            oe_info.observed_variants,
                            oe_info.expected_variants,
                            method="poisson",
                        ),
                        oe_ci_gamma=oe_confidence_interval(
                            oe_info.observed_variants,
                            oe_info.expected_variants,
                            method="gamma",
                        ),
                        z_raw=calculate_raw_z_score(
                            oe_info.observed_variants, oe_info.expected_variants
                        ),
                    )
                )
            )
        )
    )

    meta = hl.eval(ht.constraint_group_meta)
    freq_meta = hl.eval(ht.exomes_freq_meta)
    all_freq_idx = freq_meta.index(ADJ_FREQ_META)
    ht = ht.annotate(
        constraint_groups=[
            ht.constraint_groups[i].annotate(
                flags=add_filters_expr(
                    get_constraint_flags(
                        ht.constraint_groups[i].oe_info[all_freq_idx].expected_variants,
                        ht.constraint_groups[i].oe_info[all_freq_idx].z_raw,
                        z_thresholds.get(
                            "lof" if m.get("lof") else m.get("csq_set", "None"),
                            (None, None),
                        )[0],
                        z_thresholds.get(
                            "lof" if m.get("lof") else m.get("csq_set", "None"),
                            (None, None),
                        )[1],
                        flag_postfix="lof" if m.get("lof") else m.get("csq_set", None),
                    )
                )
            )
            for i, m in enumerate(meta)
        ]
    )

    return ht


def _compute_z_scores(ht: hl.Table) -> hl.Table:
    """
    Compute normalized z-scores and union per-group constraint flags.

    Computes the standard deviation of raw z-scores (stored as a global), normalizes
    each group's raw z-score by its standard deviation, and unions the syn, mis, and
    lof flags into a single ``constraint_flags`` set.

    :param ht: Table output by :func:`_annotate_oe_ci_z`.
    :return: Table with ``z_score`` and ``constraint_flags`` annotations.
    """
    meta = hl.eval(ht.constraint_group_meta)
    freq_meta = hl.eval(ht.exomes_freq_meta)
    syn_idx = meta.index({"csq_set": "syn"})
    mis_idx = meta.index({"csq_set": "mis"})
    lof_idx = meta.index({"lof": "hc"})
    all_freq_idx = freq_meta.index(ADJ_FREQ_META)

    ht = ht.annotate_globals(
        sd_raw_z=ht.aggregate(
            hl.agg.filter(
                ~ht.no_variants,
                [
                    calculate_raw_z_score_sd(
                        ht.constraint_groups[i].oe_info[all_freq_idx].z_raw,
                        ht.constraint_groups[i].flags,
                        mirror_neg_raw_z=m.get("csq_set") != "syn",
                    )
                    for i, m in enumerate(meta)
                ],
            )
        )
    )

    ht = ht.annotate(
        constraint_groups=hl.map(
            lambda x, sd_raw_z: x.annotate(
                z_score=x.oe_info[all_freq_idx].z_raw / sd_raw_z
            ),
            ht.constraint_groups,
            ht.sd_raw_z,
        ),
        constraint_flags=(
            ht.constraint_groups[syn_idx].flags
            | ht.constraint_groups[mis_idx].flags
            | ht.constraint_groups[lof_idx].flags
        ),
    )

    return ht


def _compute_percentile_bins(
    ht: hl.Table,
    use_mane_select_over_canonical: bool = True,
) -> hl.Table:
    """
    Add OE upper CI rank and percentile bin annotations.

    Adds rank and bin annotations via :func:`add_oe_upper_rank_and_bins`,
    then computes percentile thresholds across all granularities defined in
    ``CONSTRAINT_GRANULARITIES`` and annotates bins via
    :func:`annotate_constraint_percentile_bins`.

    :param ht: Table output by :func:`_compute_z_scores`.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical
        transcripts for filtering when determining ranks. Default is True.
    :return: Table with rank, decile, and percentile bin annotations.
    """
    ht = add_oe_upper_rank_and_bins(ht, use_mane_select_over_canonical)

    meta = hl.eval(ht.constraint_group_meta)
    metric_group_idx = {
        "syn": next(i for i, m in enumerate(meta) if m == {"csq_set": "syn"}),
        "mis": next(i for i, m in enumerate(meta) if m == {"csq_set": "mis"}),
        "lof": next(i for i, m in enumerate(meta) if m == {"lof": "hc"}),
    }
    outlier_expr = ht.constraint_flags.length() > 0

    gran_percentiles: Dict[str, List[float]] = {}
    all_qs = []
    for gran_name, bins in CONSTRAINT_GRANULARITIES.items():
        n_bins = len(bins) + 1
        pcts = [b / n_bins * 100 for b in bins]
        gran_percentiles[gran_name] = pcts
        all_qs.extend(pcts)

    thresholds = {}
    for metric, idx in metric_group_idx.items():
        vals = compute_percentile_thresholds(
            ht,
            percentiles=all_qs,
            metric_expr=ht.constraint_groups[idx].oe_info[0].oe_ci_gamma.upper,
            outlier_expr=outlier_expr,
            transcript_filter_expr=get_transcript_filter_expr(
                ht, mane_select_only=True
            ),
        )
        for gran_name, pcts in gran_percentiles.items():
            thresholds[(gran_name, metric)] = [vals[p] for p in pcts]

    return annotate_constraint_percentile_bins(ht, thresholds, metric_group_idx)


def _compute_pli_scores(
    ht: hl.Table,
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
) -> hl.Table:
    """
    Compute pLI, pNull, and pRec scores for the HC LoF constraint group.

    :param ht: Table output by :func:`_compute_percentile_bins`.
    :param expected_values: Dictionary containing the expected OE values for 'Null',
        'Rec', and 'LI' to use as starting values. Default is ``PLI_EXPECTED_VALUES``.
    :param min_diff_convergence: Minimum iteration change in LI to consider the EM
        model convergence criteria as met. Default is 0.001.
    :return: Table with pLI, pNull, and pRec annotations.
    """
    if expected_values is None:
        expected_values = PLI_EXPECTED_VALUES

    meta = hl.eval(ht.constraint_group_meta)
    freq_meta = hl.eval(ht.exomes_freq_meta)
    lof_idx = meta.index({"lof": "hc"})
    all_freq_idx = freq_meta.index(ADJ_FREQ_META)

    hc_lof_expr = ht.constraint_groups[lof_idx].oe_info[all_freq_idx]
    return ht.annotate(
        **compute_pli(
            ht,
            obs_expr=hc_lof_expr.observed_variants,
            exp_expr=hc_lof_expr.expected_variants,
            expected_values=expected_values,
            min_diff_convergence=min_diff_convergence,
        )
    )


def compute_constraint_metrics(
    ht: hl.Table,
    gencode_ht: hl.Table,
    gene_quality_metrics_ht: hl.Table,
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
    raw_z_outlier_threshold_lower_lof: float = -8.0,
    raw_z_outlier_threshold_lower_missense: float = -8.0,
    raw_z_outlier_threshold_lower_syn: float = -8.0,
    raw_z_outlier_threshold_upper_syn: float = 8.0,
    use_mane_select_over_canonical: bool = True,
) -> hl.Table:
    """
    Compute constraint metrics for synonymous, missense, and pLoF variants.

    Orchestrates the following steps:

    1. Annotate OE ratios, confidence intervals, raw z-scores, and per-group flags
       (:func:`_annotate_oe_ci_z`).
    2. Normalize z-scores and union constraint flags (:func:`_compute_z_scores`).
    3. Add OE upper CI rank, decile, and percentile bins
       (:func:`_compute_percentile_bins`).
    4. Compute pLI / pNull / pRec scores (:func:`_compute_pli_scores`).
    5. Annotate with gene quality metrics and GENCODE transcript annotations.

    .. note::

        The input ``ht`` should be the output of
        :func:`aggregate_by_constraint_groups`, which has a
        ``constraint_groups`` array, ``constraint_group_meta``,
        ``exomes_freq_meta``, and ``no_variants`` annotations.

    :param ht: Table output by :func:`aggregate_by_constraint_groups`.
    :param gencode_ht: Table containing GENCODE annotations.
    :param gene_quality_metrics_ht: Table keyed by transcript with
        ``gene_quality_metrics`` and ``gene_flags`` fields (output of
        :func:`compute_gene_quality_metrics`).
    :param expected_values: Dictionary containing the expected OE values for 'Null',
        'Rec', and 'LI' to use as starting values.
    :param min_diff_convergence: Minimum iteration change in LI to consider the EM
        model convergence criteria as met. Default is 0.001.
    :param raw_z_outlier_threshold_lower_lof: Lower raw z-score outlier threshold for
        LoF variants. Default is -8.0.
    :param raw_z_outlier_threshold_lower_missense: Lower raw z-score outlier threshold
        for missense variants. Default is -8.0.
    :param raw_z_outlier_threshold_lower_syn: Lower raw z-score outlier threshold for
        synonymous variants. Default is -8.0.
    :param raw_z_outlier_threshold_upper_syn: Upper raw z-score outlier threshold for
        synonymous variants. Default is 8.0.
    :param use_mane_select_over_canonical: Use MANE Select rather than canonical
        transcripts for filtering when determining ranks. Default is True.
    :return: Table with pLI scores, OE ratios, confidence intervals, z-scores,
        percentile bins, gene quality metrics, and GENCODE annotations.
    """
    z_thresholds = {
        "lof": (raw_z_outlier_threshold_lower_lof, None),
        "mis": (raw_z_outlier_threshold_lower_missense, None),
        "syn": (
            raw_z_outlier_threshold_lower_syn,
            raw_z_outlier_threshold_upper_syn,
        ),
    }

    ht = _annotate_oe_ci_z(ht, z_thresholds)
    ht = ht.checkpoint(new_temp_file("constraint_metrics.oe_ci_z", "ht"))

    ht = _compute_z_scores(ht)
    ht = ht.checkpoint(new_temp_file("constraint_metrics.z_scores", "ht"))

    ht = _compute_percentile_bins(ht, use_mane_select_over_canonical)
    ht = ht.checkpoint(new_temp_file("constraint_metrics.percentile_bins", "ht"))

    ht = _compute_pli_scores(ht, expected_values, min_diff_convergence)
    ht = ht.checkpoint(new_temp_file("constraint_metrics.pli", "ht"))

    # Add per-transcript gene quality metrics and flags.
    ht = ht.annotate(**gene_quality_metrics_ht[ht.transcript])

    # Add transcript annotations from GENCODE.
    ht = add_gencode_transcript_annotations(ht, gencode_ht)

    return ht


def _restructure_release_rows(
    ht: hl.Table,
    field_names: List[str],
    all_freq_idx: int,
    gen_anc_ds_indices: Dict[str, List[int]],
) -> hl.Table:
    """
    Restructure ``constraint_groups`` into named top-level release fields.

    For each constraint group, builds a flat release struct by:

        - Flattening the adjusted-frequency ``oe_info`` entry onto the group
        struct, keeping ``oe`` and ``z_raw`` under their original names and
        overriding ``oe_ci`` with ``oe_ci_gamma``.
        - Applying ``RELEASE_CG_RENAME`` to rename group-level and ``oe_info``
        fields (e.g. ``mu_snp`` -> ``mu``, ``observed_variants`` -> ``obs``).
        - When downsampling data is present, adding ``gen_anc_obs`` /
        ``gen_anc_exp`` structs keyed by genetic ancestry with arrays of
        values ordered by downsampling level.

    Annotates the Table with one top-level field per group (applying
    ``RELEASE_GROUP_RENAMES``, e.g. ``lof_hc`` -> ``lof``), trims the
    ``oe_ci`` struct to ranked or unranked CI fields depending on the group,
    and adds ``pLI`` / ``pNull`` / ``pRec`` for the LoF group. Finally
    selects release row fields, re-keys, and filters out transcripts with
    no possible variants in any group.

    :param ht: Table with ``constraint_groups`` and associated annotations.
    :param field_names: Internal name for each constraint group, derived
        from ``constraint_group_meta``.
    :param all_freq_idx: Index into ``oe_info`` for the adjusted allele
        frequency group.
    :param gen_anc_ds_indices: Mapping from genetic ancestry label to list
        of ``oe_info`` indices for its downsampling entries. Empty dict when
        no downsampling data is present.
    :return: Table with named top-level constraint group structs, release
        row fields selected, re-keyed, and filtered.
    """
    add_ds_fields = (
        ["observed_variants", "expected_variants"] if gen_anc_ds_indices else []
    )

    cg_expr = ht.constraint_groups.map(
        lambda cg: cg.annotate(
            **cg.oe_info[all_freq_idx],
            oe_ci=cg.oe_info[all_freq_idx].oe_ci_gamma,
            # Per-genetic-ancestry downsampling obs/exp arrays, one value per
            # downsampling level, keyed by genetic ancestry.
            **{
                f"gen_anc_{RELEASE_CG_RENAME[f]}": hl.struct(
                    **{
                        gen_anc: hl.array([cg.oe_info[j][f] for j in indices])
                        for gen_anc, indices in gen_anc_ds_indices.items()
                    }
                )
                for f in add_ds_fields
            },
        )
    )
    cg_select = (
        RELEASE_CG_SELECT
        if gen_anc_ds_indices
        else [f for f in RELEASE_CG_SELECT if not f.startswith("gen_anc_")]
    )
    cg_expr = cg_expr.map(
        lambda cg: cg.annotate(
            **{RELEASE_CG_RENAME[k]: cg[k] for k in RELEASE_CG_RENAME}
        ).select(*cg_select)
    )

    # Build a top-level release field for each constraint group, applying
    # group renames (e.g. lof_hc -> lof), trimming oe_ci sub-fields, and
    # adding pLI/pNull/pRec for the LoF group.
    cg_fields = {
        RELEASE_GROUP_RENAMES.get(name, name): cg_expr[i]
        .annotate(
            oe_ci=cg_expr[i].oe_ci.select(
                *(
                    RELEASE_CI_FIELDS_WITH_RANK
                    if name in RELEASE_GROUPS_WITH_RANK
                    else RELEASE_CI_FIELDS
                )
            ),
            **{
                k: ht[k]
                for k in (RELEASE_LOF_FIELDS if name in RELEASE_GROUPS_WITH_PLI else [])
            },
        )
        .select()  # Drop intermediate fields, keep only annotated release fields.
        for i, name in enumerate(field_names)
    }
    ht = ht.annotate(**cg_fields)
    ht = ht.select(*RELEASE_TOP_LEVEL_ANNOTATIONS, *RELEASE_GROUP_NAMES)

    available_keys = [k for k in RELEASE_KEY_ORDER if k in ht.key]
    if list(ht.key) != available_keys:
        ht = ht.key_by(*available_keys)

    ht = ht.filter(hl.any([ht[k].possible != 0 for k in RELEASE_GROUP_NAMES]))

    return ht


def _restructure_release_globals(
    ht: hl.Table,
    field_names: List[str],
    freq_meta: List[Dict],
    gen_anc_ds_indices: Dict[str, List[int]],
    sd_raw_z_arr: List[float],
    release_version: Optional[str],
) -> hl.Table:
    """
    Restructure globals for public release.

    Replaces internal globals with a clean release set:

        - Pipeline parameter globals are renamed and stripped of internal-only
        fields via ``RELEASE_PIPELINE_PARAM_GLOBALS``.
        - ``sd_raw_z`` is converted from an ordered array (one entry per
        constraint group) to a named struct keyed by release group name,
        retaining only groups in ``RELEASE_GROUP_NAMES``.
        - When downsampling data is present, a ``downsamplings`` struct is added
        keyed by genetic ancestry, with arrays of integer downsampling levels
        matching the order of ``gen_anc_obs`` / ``gen_anc_exp`` in the rows.
        - ``max_af`` is preserved unchanged if present.
        - ``version`` is set to ``release_version`` if provided, otherwise
        carried over from the existing global.

    :param ht: Table whose globals are being restructured.
    :param field_names: Internal name for each constraint group (parallel to
        ``sd_raw_z_arr``), used to map array positions to release group names.
    :param freq_meta: Evaluated ``exomes_freq_meta`` global, used to extract
        downsampling levels for each genetic ancestry.
    :param gen_anc_ds_indices: Mapping from genetic ancestry label to list
        of ``oe_info`` indices for its downsampling entries. Empty dict when
        no downsampling data is present.
    :param sd_raw_z_arr: Evaluated ``sd_raw_z`` global array, parallel to
        ``field_names``.
    :param release_version: Version string for the ``version`` global. When
        *None*, the existing ``version`` global is retained if present.
    :return: Table with release-formatted globals.
    """
    sd_raw_z_name_map = {n: RELEASE_GROUP_RENAMES.get(n, n) for n in field_names}
    sd_raw_z_struct = hl.struct(
        **{
            sd_raw_z_name_map[field_names[i]]: sd_raw_z_arr[i]
            for i in range(len(field_names))
            if sd_raw_z_name_map[field_names[i]] in RELEASE_GROUP_NAMES
        }
    )

    global_kwargs = {}
    if release_version is not None:
        global_kwargs["version"] = release_version
    elif "version" in ht.globals:
        global_kwargs["version"] = ht.globals.version

    for src, dest, drop_fields in RELEASE_PIPELINE_PARAM_GLOBALS:
        if src in ht.globals:
            global_kwargs[dest] = ht.globals[src].drop(*drop_fields)

    if gen_anc_ds_indices:
        global_kwargs["downsamplings"] = hl.struct(
            **{
                gen_anc: [int(freq_meta[j]["downsampling"]) for j in indices]
                for gen_anc, indices in gen_anc_ds_indices.items()
            }
        )

    if "max_af" in ht.globals:
        global_kwargs["max_af"] = ht.globals.max_af

    global_kwargs["sd_raw_z"] = sd_raw_z_struct
    return ht.select_globals(**global_kwargs)


def prepare_release_ht(
    ht: hl.Table,
    release_version: Optional[str] = None,
) -> hl.Table:
    """
    Prepare the constraint metrics Table for public release.

    Computes shared metadata needed by both restructuring steps, then
    delegates row and global restructuring to
    :func:`_restructure_release_rows` and
    :func:`_restructure_release_globals`.

    The internal ``constraint_groups`` schema has:

        - Group-level fields: ``mu_snp``, ``possible_variants``, ``z_score``.
        - Per-frequency ``oe_info`` array (one entry per ``exomes_freq_meta``
          element): ``observed_variants``, ``expected_variants``, ``oe``,
          ``oe_ci_gamma``, ``z_raw``.

    The release schema exposes one top-level struct per group
    (``syn``, ``mis``, ``lof_hc_lc``, ``lof``; ``lof_hc`` is renamed to
    ``lof``), with fields ``mu``, ``possible``, ``obs``, ``exp``, ``oe``,
    ``oe_ci``, ``z_raw``, ``z_score``, and optionally ``gen_anc_obs`` /
    ``gen_anc_exp`` when downsampling data is present.

    :param ht: Internal constraint metrics Table (output of
        ``compute_constraint_metrics``). Expected to already contain GENCODE
        transcript annotations (``transcript_id_version``, ``level``, etc.)
        and gene quality metric annotations (``gene_quality_metrics``,
        ``gene_flags``).
    :param release_version: Version string for the ``version`` global. When
        *None*, the existing ``version`` global is retained if present.
    :return: Release-formatted Table.
    """
    ht = ht.rename(GENCODE_FIELD_RENAMES)

    constraint_meta = hl.eval(ht.constraint_group_meta)
    freq_meta = hl.eval(ht.exomes_freq_meta)
    all_freq_idx = freq_meta.index(ADJ_FREQ_META)

    field_names = [
        "_".join(f"{k}_{v}" for k, v in m.items()).replace("csq_set_", "")
        for m in constraint_meta
    ]
    logger.info("Release constraint group field names: %s", field_names)

    gen_anc_ds_indices: Dict[str, List[int]] = {}
    if "downsamplings" in ht.globals:
        for j, m in enumerate(freq_meta):
            gen_anc = m.get("gen_anc")
            if gen_anc is not None and "downsampling" in m:
                gen_anc_ds_indices.setdefault(gen_anc, []).append(j)

    # Evaluate sd_raw_z before the row select (globals persist through it).
    sd_raw_z_arr = hl.eval(ht.sd_raw_z)

    ht = _restructure_release_rows(ht, field_names, all_freq_idx, gen_anc_ds_indices)
    ht = _restructure_release_globals(
        ht, field_names, freq_meta, gen_anc_ds_indices, sd_raw_z_arr, release_version
    )
    return ht


def flatten_release_ht(ht: hl.Table) -> hl.Table:
    """
    Flatten the release constraint metrics Table for TSV export.

    Drops per-genetic-ancestry downsampling fields (``gen_anc_obs``,
    ``gen_anc_exp``) when present and calls :meth:`~hail.Table.flatten`
    to expand nested struct fields using ``.`` as the separator
    (e.g. ``lof.obs``, ``lof.oe_ci.upper``).

    :param ht: Release-format constraint metrics Table (output of
        :func:`prepare_release_ht`).
    :return: Flat Table suitable for :meth:`~hail.Table.export`.
    """
    # Drop struct/array fields not suitable for flat TSV export.
    drop_fields = [f for f in ["gen_anc_obs", "gen_anc_exp"] if f in ht.row]
    if drop_fields:
        ht = ht.drop(*drop_fields)

    return ht.flatten()


def annotate_constraint_percentile_bins(
    ht: hl.Table,
    thresholds: Dict[Tuple[str, str], List[float]],
    metric_group_idx: Dict[str, int],
) -> hl.Table:
    """
    Annotate each transcript with its percentile bin for all metric/granularity combinations.

    Thin wrapper around :func:`annotate_bins_by_threshold` that extracts the
    gamma upper CI value from each constraint group's first oe_info element.

    Annotates ``constraint_bins.{granularity}.{metric}`` for each combination.
    Bin 0 is the most constrained (value below all thresholds); bin N equals
    the number of boundaries the value exceeds.

    :param ht: Constraint metrics Table with a ``constraint_groups`` array field.
    :param thresholds: Mapping of ``(granularity, metric)`` to an ordered list
        of threshold values, as produced by
        :func:`compute_percentile_thresholds`.
    :param metric_group_idx: Mapping of metric name to its index in
        ``constraint_groups`` (e.g. ``{"lof": 5, "mis": 1, "syn": 0}``).
    :return: Annotated Table with an added ``constraint_bins`` struct field.
    """
    logger.info(
        "Annotating bins for %d (granularity, metric) combinations.",
        len(thresholds),
    )

    metric_exprs = {
        metric: ht.constraint_groups[idx].oe_info[0].oe_ci_gamma.upper
        for metric, idx in metric_group_idx.items()
    }

    return annotate_bins_by_threshold(
        ht,
        metric_exprs=metric_exprs,
        thresholds=thresholds,
        granularities=list(CONSTRAINT_GRANULARITIES),
    )
