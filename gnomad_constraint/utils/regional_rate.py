"""
Functions to build the training set of the regional mutation rate (rr) models.

The models are fit, and rr is computed and normalized, with
:mod:`gnomad.utils.constraint.regional_mutation_rate`.
"""

import logging

import hail as hl
import pandas as pd
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.constants import GENOMIC_FEATURE_COLUMNS
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
) -> hl.Table:
    """
    Create the rr training Table of synonymous DNM and non-mutated background sites.

    The DNM sites are the synonymous sites of ``context_ht`` that are in ``dnm_ht``.
    The background is a random sample (``background_fraction``) of the synonymous
    sites absent from gnomAD exomes (``filters.exomes`` is missing) at the other
    loci, with one site kept per locus, so no locus is both a DNM and a background
    site. DNM sites get ``y`` = 1 and background
    sites ``y`` = 0, and every site gets the features of the 1kb window containing it.
    The Table has no globals.

    Each DNM variant is one row, however many DNMs were reported there, since the DNM
    Table can't tell a recurrent DNM from one reported by two cohorts.

    :param context_ht: Annotated context Table with ``vep``, ``context``, and
        ``filters``.
    :param dnm_ht: DNM Table keyed by ``locus`` and ``alleles``.
    :param features_ht: Genomic feature Table keyed by ``interval``.
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
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


def rr_training_ht_to_pandas(ht: hl.Table) -> pd.DataFrame:
    """
    Collect the rr training Table into a pandas DataFrame with a tuning fold.

    The fold splits the chromosomes into ``odd`` (chr1, chr3, ...) and ``even``. The
    genomic features cover chr1-22 only, so chrX and chrY sites, which have no
    features, are left out of every fit and score whatever their fold.

    :param ht: rr training Table (:func:`create_rr_training_ht`).
    :return: DataFrame with ``context``, ``y``, ``GENOMIC_FEATURE_COLUMNS``, and
        ``fold``.
    """
    ht = ht.key_by()
    ht = ht.select(
        "context",
        "y",
        *GENOMIC_FEATURE_COLUMNS,
        fold=hl.if_else(ht.locus.contig_idx % 2 == 0, "odd", "even"),
    )
    # to_pandas returns nullable extension dtypes (e.g. Float64); use NumPy dtypes
    # (missing values as NaN).
    return ht.to_pandas().astype(
        {"y": "int64", **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, "float64")}
    )
