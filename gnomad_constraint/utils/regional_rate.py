"""
Functions to build the regional mutation rate (rr) training set and compute rr.

The models are fit, and rr is computed and normalized, with
:mod:`gnomad.utils.constraint.regional_mutation_rate`.
"""

import logging
from typing import Optional

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
    RR_CONTEXTS,
)
from gnomad_constraint.utils.constraint import get_transcript_filter_expr

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_utils")
logger.setLevel(logging.INFO)


def create_rr_training_ht(
    context_ht: hl.Table,
    dnm_ht: hl.Table,
    features_ht: hl.Table,
    background_fraction: float,
    seed: int,
    gene_flags_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Create the rr training Table of synonymous DNM and non-mutated background sites.

    The DNM sites are the synonymous sites of ``context_ht`` that are in ``dnm_ht``.
    The background is a random sample (``background_fraction``) of the synonymous
    sites absent from gnomAD exomes (``filters.exomes`` is missing) at the other
    loci, with one site kept per locus, so no locus is both a DNM and a background
    site. DNM sites get ``y`` = 1 and background
    sites ``y`` = 0, and every site gets the features of the 1kb window containing it.
    With ``gene_flags_ht``, sites in genes with a QC gene flag are left out
    (:func:`remove_qc_flagged_genes`). The Table has no globals.

    Each DNM variant is one row, however many DNMs were reported there, since the DNM
    Table can't tell a recurrent DNM from one reported by two cohorts.

    :param context_ht: Annotated context Table with ``vep``, ``context``, and
        ``filters``.
    :param dnm_ht: DNM Table keyed by ``locus`` and ``alleles``.
    :param features_ht: Genomic feature Table keyed by ``interval``.
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
    :param gene_flags_ht: Constraint metrics Table with gene flags, whose flagged
        genes' sites are left out (see :func:`remove_qc_flagged_genes`), or None to
        keep every gene. Default is None.
    :return: Training Table keyed by ``locus`` with ``alleles``, ``context``, ``y``,
        and ``GENOMIC_FEATURE_COLUMNS``.
    """
    ht = context_ht.filter(
        context_ht.vep.most_severe_consequence == "synonymous_variant"
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

    return ht.select_globals()


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
    features_ht: hl.Table, models_ht: hl.Table, per_variant_ht: hl.Table
) -> hl.Table:
    """
    Compute rr per 1kb window and context, normalized on the per-variant Table's sites.

    Each window's rr in each of ``RR_CONTEXTS`` comes from that context's model, and
    windows no model covers are dropped. rr is then divided, per context, by its mean
    over the sites of :func:`get_rr_normalization_sites_ht`, weighted by their adj
    expected count, so it moves those sites' expected counts between windows without
    changing their total.

    :param features_ht: Genomic feature Table keyed by ``interval``.
    :param models_ht: rr models Table with a ``ridge_penalty`` global.
    :param per_variant_ht: Per-variant expected Table to normalize on, with an
        ``apply_models_globals`` global (see :func:`get_rr_normalization_sites_ht`).
    :return: Table keyed by ``interval`` with ``adj_r``, a dict of context to rr, and
        an ``rr_globals`` global with ``ridge_penalty``, ``normalization_factors``,
        and ``per_variant_apply_models_globals``.
    """
    # TODO: Windows missing a context's features get 1.0 here, the model's rate for
    #  the average window. Normalization then divides it by the context's factor, so
    #  these windows end up at 1 / factor instead of 1, unlike sites without a
    #  window. Consider missing_feature_rate=None, which leaves them missing, and
    #  filling 1.0 after normalization.
    ht = features_ht.select(
        adj_r=regional_mutation_rate_expr(features_ht.row, models_ht, RR_CONTEXTS)
    )
    ht = ht.filter(hl.is_defined(ht.adj_r))
    ht = ht.checkpoint(new_temp_file("rr_unnormalized", "ht"))

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

    return ht.annotate_globals(
        rr_globals=hl.struct(
            ridge_penalty=models_ht.index_globals().ridge_penalty,
            normalization_factors=hl.literal(factors, hl.tdict(hl.tstr, hl.tfloat64)),
            # --apply-regional-rr checks that the per-variant Table it aggregates
            # still has these models.
            per_variant_apply_models_globals=(
                per_variant_ht.index_globals().apply_models_globals
            ),
        )
    )
