"""
Functions to build the regional mutation rate (rr) site features, select features, and compute rr.

The models are fit, and rr is computed and normalized, with
:mod:`gnomad.utils.constraint.regional_mutation_rate`.
"""

import logging
from typing import Dict, List, Optional

import hail as hl
from gnomad.utils.constraint import (
    normalize_regional_mutation_rate_expr,
    regional_mutation_rate_expr,
    regional_mutation_rate_normalization_factors_agg_expr,
)
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_SCALES,
    RR_BONFERRONI_CUTOFF,
    RR_CONTEXT_FEATURES,
    RR_CONTEXTS,
    RR_FEATURE_SELECTION_ALPHA,
)
from gnomad_constraint.utils.constraint import get_transcript_filter_expr

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_utils")
logger.setLevel(logging.INFO)


def create_rr_site_features_ht(
    preprocessed_ht: hl.Table,
    dnm_ht: hl.Table,
    features_ht: hl.Table,
    background_fraction: float,
    seed: int,
    gene_flags_ht: Optional[hl.Table] = None,
    gene_flags_version: Optional[str] = None,
) -> hl.Table:
    """
    Create the rr site features: synonymous DNM and non-mutated background sites.

    Only synonymous sites in the plateau models' high coverage range are used:
    ``exomes_coverage`` of at least ``build_models_globals.high_cov_cutoff`` (and at
    most ``upper_cov_cutoff``, if set), as for the sites the plateau models are
    trained on. Elsewhere a DNM can't be called in the exome cohorts and a site is
    absent from gnomAD because it isn't covered, so a site there could only be a
    background site, and the models would learn the features of poorly covered
    regions as a low mutation rate.

    A site is synonymous if that is its most severe consequence in any transcript,
    which keeps more DNMs for training than the Ensembl MANE Select synonymous sites
    rr is normalized on (:func:`get_rr_normalization_sites_ht`), the transcripts the
    constraint metrics rank.

    The DNM sites are those sites that are in ``dnm_ht``. The background is a
    random sample (``background_fraction``) of those sites absent from gnomAD
    exomes (``filters.exomes`` is missing) at the other loci, with one site kept per
    locus, so no locus is both a DNM and a background site. DNM sites get ``y`` = 1
    and background sites ``y`` = 0, and every site gets the features of the 1kb
    window containing it.
    With ``gene_flags_ht``, sites in genes with a QC gene flag are left out
    (:func:`remove_qc_flagged_genes`). The settings that chose the sites are
    recorded in an ``rr_training_sites`` global, which the later rr steps carry
    to the rr Table and the release.

    Each DNM variant is one row, however many DNMs were reported there, since the DNM
    Table can't tell a recurrent DNM from one reported by two cohorts.

    A preprocessed Table that doesn't record its exome coverage metric, written
    before ``apply_models_globals`` had ``exome_coverage_metric``, raises an error,
    since the coverage range of the sites couldn't be recorded.

    :param preprocessed_ht: Preprocessed Table with ``vep``, ``context``,
        ``filters``, ``exomes_coverage``, a ``build_models_globals`` global with
        ``high_cov_cutoff`` and ``upper_cov_cutoff``, and an
        ``apply_models_globals`` global with ``exome_coverage_metric``.
    :param dnm_ht: DNM Table keyed by ``locus`` and ``alleles``.
    :param features_ht: Genomic feature Table keyed by ``interval``.
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
    :param gene_flags_ht: Constraint metrics Table with gene flags, whose flagged
        genes' sites are left out (see :func:`remove_qc_flagged_genes`), or None to
        keep every gene. Default is None.
    :param gene_flags_version: Version of the constraint metrics release
        ``gene_flags_ht`` comes from, recorded in the globals. Default is None.
    :return: Site features Table keyed by ``locus`` with ``alleles``, ``context``,
        ``y``, and ``GENOMIC_FEATURE_COLUMNS``, and an ``rr_training_sites`` global
        with ``background_fraction``, ``background_seed``,
        ``qc_flagged_gene_filter``, ``gene_flags_version``, and the coverage range
        (``exome_coverage_metric``, ``high_cov_cutoff``, ``upper_cov_cutoff``).
    """
    if "exome_coverage_metric" not in preprocessed_ht.apply_models_globals:
        raise ValueError(
            "The preprocessed Table doesn't record its exome coverage metric in"
            " apply_models_globals. Rerun prepare_data.py --preprocess-data."
        )

    coverage = preprocessed_ht.exomes_coverage
    build_globals = preprocessed_ht.build_models_globals
    ht = preprocessed_ht.filter(
        (preprocessed_ht.vep.most_severe_consequence == "synonymous_variant")
        & (coverage >= build_globals.high_cov_cutoff)
        & (
            hl.is_missing(build_globals.upper_cov_cutoff)
            | (coverage <= build_globals.upper_cov_cutoff)
        )
    )
    ht = ht.select("context", in_exomes=hl.is_defined(ht.filters.exomes))
    ht = ht.checkpoint(new_temp_file("rr_syn_sites", "ht"))

    dnm_sites_ht = ht.semi_join(dnm_ht).key_by("locus")
    # TODO: Whether a site is absent from gnomAD depends on its own mutation rate
    #  (about half of CpG synonymous sites are seen in gnomAD exomes) and on coverage,
    #  so this background makes the models fit the DNM rate ratio times the ratio of
    #  absence probabilities, which steepens the regional effect, most at CpG. It
    #  partly offsets the compression of the rr form (see regional_mutation_rate_expr
    #  in gnomad.utils.constraint). Consider sampling from all non-DNM synonymous
    #  sites; the ridge penalty would need retuning.
    bg_ht = ht.filter(~ht.in_exomes).key_by("locus").distinct()
    bg_ht = bg_ht.anti_join(dnm_sites_ht).sample(background_fraction, seed=seed)

    ht = dnm_sites_ht.select("alleles", "context", y=1).union(
        bg_ht.select("alleles", "context", y=0)
    )
    if gene_flags_ht is not None:
        ht = remove_qc_flagged_genes(ht, gene_flags_ht)

    ht = ht.annotate(**features_ht.select(*GENOMIC_FEATURE_COLUMNS)[ht.locus])
    preprocessed_globals = preprocessed_ht.index_globals()

    return ht.select_globals(
        rr_training_sites=hl.struct(
            background_fraction=hl.float64(background_fraction),
            background_seed=hl.int32(seed),
            qc_flagged_gene_filter=gene_flags_ht is not None,
            gene_flags_version=hl.or_missing(
                gene_flags_ht is not None, hl.str(gene_flags_version or "")
            ),
            exome_coverage_metric=(
                preprocessed_globals.apply_models_globals.exome_coverage_metric
            ),
            high_cov_cutoff=preprocessed_globals.build_models_globals.high_cov_cutoff,
            upper_cov_cutoff=(
                preprocessed_globals.build_models_globals.upper_cov_cutoff
            ),
        )
    )


def remove_qc_flagged_genes(ht: hl.Table, metrics_ht: hl.Table) -> hl.Table:
    """
    Remove sites in genes with a release QC gene flag from a locus-keyed Table.

    Genes flagged ``low_exome_coverage`` or ``low_exome_mapping_quality`` lose
    variant calls in both the DNM cohorts and gnomAD, so leaving their sites in the
    rr training set can teach the models that those regions mutate less. The span
    (``start_position`` to ``end_position``) of every Ensembl MANE Select transcript
    with a non-empty ``gene_flags`` is removed.

    :param ht: Table keyed by ``locus``.
    :param metrics_ht: Constraint metrics Table with ``transcript``, ``mane_select``,
        ``gene_flags``, ``chromosome``, ``start_position``, and ``end_position``.
    :return: ``ht`` without the sites in flagged genes.
    """
    metrics_ht = metrics_ht.filter(
        get_transcript_filter_expr(metrics_ht, mane_select_only=True)
        & (hl.len(metrics_ht.gene_flags) > 0)
    )
    span_expr = hl.locus_interval(
        metrics_ht.chromosome,
        metrics_ht.start_position,
        metrics_ht.end_position,
        includes_end=True,
        reference_genome="GRCh38",
    )
    intervals = metrics_ht.aggregate(
        hl.agg.filter(hl.is_defined(span_expr), hl.agg.collect(span_expr))
    )
    logger.info("Removing sites in %d QC-flagged genes.", len(intervals))
    if not intervals:
        return ht

    return hl.filter_intervals(ht, intervals, keep=False)


def get_rr_fold_expr(locus_expr: hl.expr.LocusExpression) -> hl.expr.StringExpression:
    """
    Get the chromosome-split fold used to tune the rr ridge penalty.

    The fold splits the chromosomes into ``odd`` (chr1, chr3, ...) and ``even``. The
    genomic features cover chr1-22 only, so chrX and chrY sites, which have no
    features, are left out of every fit and score whatever their fold.

    :param locus_expr: Locus of each training site.
    :return: ``odd`` or ``even``.
    """
    return hl.if_else(locus_expr.contig_idx % 2 == 0, "odd", "even")


def select_rr_features(screen_ht: hl.Table, cutoff: str) -> Dict[str, List[str]]:
    """
    Select each context's rr model features from the univariate feature screen.

    This is Gnocchi's feature selection. With ``RR_BONFERRONI_CUTOFF``, a context's
    feature is selected when its p-value is at most
    ``RR_FEATURE_SELECTION_ALPHA / len(GENOMIC_FEATURE_SCALES) / n``, where ``n`` is
    the number of features in the context's ``RR_CONTEXT_FEATURES`` at one scale (13,
    or 8 in CpG contexts). Otherwise ``cutoff`` is the p-value cutoff itself. For each
    feature only the smallest selected scale is kept. A screen fit that didn't converge
    has no p-value, so it is never selected, and a context with no selected feature is
    left out, so it gets no model. The screen must be of ``RR_CONTEXT_FEATURES``, or a
    feature it didn't screen would go unselected without an error, so a screen of other
    features raises one.

    :param screen_ht: Feature screen Table keyed by ``context`` and ``feature``, with
        ``p_value`` and a ``features`` global of the features screened
        (``screen_regional_mutation_rate_features`` in gnomad_methods).
    :param cutoff: ``RR_BONFERRONI_CUTOFF``, or a p-value cutoff such as ``"1e-2"``.
    :return: Dictionary of context to its selected features, in the order of its
        ``RR_CONTEXT_FEATURES``.
    """
    if (
        "features" not in screen_ht.globals
        or hl.eval(screen_ht.features) != RR_CONTEXT_FEATURES
    ):
        raise ValueError(
            "The feature screen isn't a screen of RR_CONTEXT_FEATURES. Rerun"
            " expected_models.py --screen-rr-features."
        )
    p_values = {
        (r.context, r.feature): r.p_value for r in screen_ht.select("p_value").collect()
    }
    selected = {}
    for context, features in RR_CONTEXT_FEATURES.items():
        if cutoff == RR_BONFERRONI_CUTOFF:
            n_features = len(features) // len(GENOMIC_FEATURE_SCALES)
            threshold = RR_FEATURE_SELECTION_ALPHA / len(GENOMIC_FEATURE_SCALES)
            threshold /= n_features
        else:
            threshold = float(cutoff)
        # RR_CONTEXT_FEATURES lists every feature at one scale before the next, from
        # the smallest, so the first selected column of a feature is its smallest.
        kept, seen = [], set()
        for column in features:
            feature = column.rsplit("_", 1)[0]
            p_value = p_values.get((context, column))
            if feature in seen or p_value is None or p_value > threshold:
                continue
            seen.add(feature)
            kept.append(column)
        if kept:
            selected[context] = kept

    return selected


def get_rr_normalization_sites_ht(per_variant_ht: hl.Table) -> hl.Table:
    """
    Get the sites rr is normalized on: the exome synonymous Ensembl MANE Select sites.

    The RefSeq MANE Select transcripts are the same transcripts, so keeping them
    would count every site twice.

    :param per_variant_ht: Per-variant expected Table with ``annotation``,
        ``transcript``, ``mane_select``, ``context``, and ``expected_variants``, and
        an ``exomes_freq_meta`` global.
    :return: Table keyed as ``per_variant_ht`` with ``context`` and ``exp``, the adj
        expected count.
    """
    adj_idx = hl.eval(per_variant_ht.exomes_freq_meta).index(ADJ_FREQ_META)
    ht = per_variant_ht.filter(
        (per_variant_ht.annotation == "synonymous_variant")
        & get_transcript_filter_expr(per_variant_ht, mane_select_only=True)
    )

    return ht.select("context", exp=ht.expected_variants[adj_idx])


def compute_rr_ht(
    features_ht: hl.Table,
    models_ht: hl.Table,
    per_variant_ht: hl.Table,
    normalize: bool = True,
) -> hl.Table:
    """
    Compute rr per 1kb window and context, normalized on the per-variant Table's sites.

    Each window's rr in each of ``RR_CONTEXTS`` comes from that context's model, and
    windows no model covers are dropped. A window missing a feature of a context's
    model has no rr in that context, so its expected counts there aren't adjusted.
    A context whose model didn't converge gets rr 1.0, as if it had no model,
    rather than a fit that may be far from the best one. rr is then divided, per
    context, by its mean over the sites of
    :func:`get_rr_normalization_sites_ht` that have one, weighted by their adj
    expected count, so it moves those sites' expected counts between windows without
    changing their total. Without ``normalize``, rr is left as the models give it, as
    Gnocchi's was, and no normalization factors are computed.

    :param features_ht: Genomic feature Table keyed by ``interval``.
    :param models_ht: rr models Table with ``ridge_penalty``, ``rr_training_sites``,
        and ``rr_feature_selection`` globals.
    :param per_variant_ht: Per-variant expected Table to normalize on, with an
        ``apply_models_globals`` global (see :func:`get_rr_normalization_sites_ht`),
        as read from disk, so counting its rows reads only metadata.
    :param normalize: Whether to normalize rr per context. Default is True.
    :return: Table keyed by ``interval`` with ``adj_r``, a dict of context to rr, and
        an ``rr_globals`` global with ``ridge_penalty``, ``feature_selection`` (the
        models' ``rr_feature_selection``), ``normalized``, ``normalization_factors``
        (empty without ``normalize``), ``training_sites`` (the models'
        ``rr_training_sites``), ``non_converged_contexts``,
        ``per_variant_apply_models_globals``, and ``per_variant_row_count``.
    """
    # The models Table has one small row per context.
    non_converged = models_ht.aggregate(
        hl.agg.filter(~models_ht.converged, hl.agg.collect(models_ht.context))
    )
    if non_converged:
        logger.warning(
            "The rr models of %s didn't converge, so those contexts get rr 1.",
            sorted(non_converged),
        )
    ht = features_ht.select(
        adj_r=regional_mutation_rate_expr(
            features_ht.row, models_ht, RR_CONTEXTS, skip_non_converged=True
        )
    )
    ht = ht.filter(hl.is_defined(ht.adj_r))
    ht = ht.checkpoint(new_temp_file("rr_unnormalized", "ht"))

    factors = {}
    if normalize:
        sites_ht = get_rr_normalization_sites_ht(per_variant_ht)
        sites_ht = sites_ht.checkpoint(new_temp_file("rr_normalization_sites", "ht"))
        sites_ht = sites_ht.naive_coalesce(200)
        factors = sites_ht.aggregate(
            regional_mutation_rate_normalization_factors_agg_expr(
                sites_ht.context,
                ht[sites_ht.locus].adj_r.get(sites_ht.context),
                sites_ht.exp,
            )
        )
        logger.info("rr normalization factors: %s", factors)
        ht = ht.annotate(adj_r=normalize_regional_mutation_rate_expr(ht.adj_r, factors))

    models_globals = models_ht.index_globals()
    return ht.annotate_globals(
        rr_globals=hl.struct(
            ridge_penalty=models_globals.ridge_penalty,
            feature_selection=models_globals.rr_feature_selection,
            normalized=normalize,
            normalization_factors=hl.literal(factors, hl.tdict(hl.tstr, hl.tfloat64)),
            training_sites=models_globals.rr_training_sites,
            non_converged_contexts=hl.literal(
                sorted(non_converged), hl.tarray(hl.tstr)
            ),
            # --apply-regional-rr checks that the per-variant Table it aggregates
            # still has these models and this many rows.
            per_variant_apply_models_globals=(
                per_variant_ht.index_globals().apply_models_globals
            ),
            per_variant_row_count=hl.int64(per_variant_ht.count()),
        )
    )
