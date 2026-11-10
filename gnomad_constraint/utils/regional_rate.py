"""
Utility functions for the regional mutation rate (rr) adjustment of expected counts.

The rr adjustment is a port of the Gnocchi regional adjustment
(https://github.com/atgu/gnomad_nc_constraint) retrained on synonymous de novo
mutations (DNMs):

    - For each strand-collapsed trinucleotide context, a ridge-penalized logistic
      model separates synonymous DNM sites from a sample of non-mutated synonymous
      sites, using 13 regional genomic features of the 1kb window containing the
      site, each at 4 flanking scales.
    - The rr of a window is the model's predicted probability for the window divided
      by its prediction for the average window, then divided by its mean over exome
      synonymous sites in each context, so rr moves expected counts between regions
      without changing each context's total.

rr is applied to the per-variant expected counts after model calibration; see
:func:`gnomad_constraint.utils.constraint.aggregate_per_variant_expected_ht`.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
import numpy as np
import pandas as pd
import scipy.stats
from gnomad.utils.liftover import default_lift_data
from hail.utils.misc import new_temp_file
from sklearn.decomposition import IncrementalPCA
from sklearn.linear_model import LogisticRegression

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    CPG_CONTEXTS,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_SCALES,
    METHYLATION_CORRELATED_FEATURES,
    RR_CONTEXTS,
    RR_FEATURE_ORDER,
    RR_RIDGE_PENALTY,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_utils")
logger.setLevel(logging.INFO)

SYNONYMOUS_CSQ = "synonymous_variant"
"""VEP most severe consequence used to define synonymous sites and DNMs."""

ENSEMBL_TRANSCRIPT_PREFIX = "ENST"
"""Prefix of Ensembl transcript IDs (RefSeq MANE Select transcripts are excluded)."""


# ---------------------------------------------------------------------------
# Synonymous de novo mutations
# ---------------------------------------------------------------------------


def format_dnm_cohort(cohort: str, ht: hl.Table) -> hl.Table:
    """
    Format one imported DNM cohort as a GRCh38 Table keyed by ``locus`` and ``alleles``.

    Supported cohorts (see
    :func:`gnomad_constraint.resources.resource_utils.get_dnm_cohort_resources`):

        - ``decode_psychencode``: deCODE + PsychENCODE DNMs used by Gnocchi
          (GRCh38). Each DNM is listed in both allele orientations; the wrong
          orientation is removed when joining to the context Table.
        - ``satterstrom_2026_asd``, ``satterstrom_2026_ddd``: Satterstrom 2026 ASD
          and DD DNMs (GRCh38, with a "chrom:pos:ref:alt" ``Variant`` column).
          Indels are removed when joining to the context Table.
        - ``gnomad_v4.1``: gnomAD v4.1 high-quality coding de novos (GRCh38 Hail
          Table).
        - ``chd_2017``, ``ee_2014``, ``cdh_2021``: GRCh37, lifted over to GRCh38.
        - ``ofc_2020``: GRCh38 positions.

    :param cohort: Cohort name.
    :param ht: Imported cohort Table.
    :return: Table keyed by ``locus`` and ``alleles`` with a ``cohort`` field.
    """
    if cohort == "gnomad_v4.1":
        ht = ht.key_by("locus", "alleles").select()
    elif cohort in {"satterstrom_2026_asd", "satterstrom_2026_ddd"}:
        variant = ht.Variant.split(":")
        ht = ht.select(
            locus=hl.locus(variant[0], hl.int(variant[1]), reference_genome="GRCh38"),
            alleles=variant[2:],
        ).key_by("locus", "alleles")
    elif cohort == "decode_psychencode":
        ht = ht.select(
            locus=hl.parse_locus(ht.locus, reference_genome="GRCh38"),
            alleles=hl.array([ht.ref, ht.alt]),
        ).key_by("locus", "alleles")
    elif cohort == "chd_2017":
        ht = ht.select(
            locus=hl.parse_locus(
                ht.CHROM + ":" + hl.str(ht.POS), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.REF, ht.ALT]),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    elif cohort == "ofc_2020":
        ht = ht.select(
            locus=hl.parse_locus(
                "chr" + hl.str(ht.Chromosome) + ":" + hl.str(ht["Position (build 38)"]),
                reference_genome="GRCh38",
            ),
            alleles=hl.array([ht.Ref, ht.Alt]),
        ).key_by("locus", "alleles")
    elif cohort == "ee_2014":
        ht = ht.select(
            locus=hl.parse_locus(
                ht["hg19 coordinates (chr:position)"].replace("chr", ""),
                reference_genome="GRCh37",
            ),
            alleles=ht["Ref/Alt alleles"].split("/"),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    elif cohort == "cdh_2021":
        ht = ht.select(
            locus=hl.parse_locus(
                hl.str(ht.chr) + ":" + hl.str(ht.position), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.ref, ht.alt]),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    else:
        raise ValueError(f"Unknown DNM cohort: {cohort}!")

    return ht.annotate(cohort=cohort)


def combine_syn_dnms(
    cohort_hts: Dict[str, hl.Table],
    context_ht: hl.Table,
) -> Tuple[hl.Table, Dict[str, int]]:
    """
    Combine DNM cohorts, filter to synonymous DNMs, and deduplicate by locus/alleles.

    DNMs are annotated with the strand-collapsed ``context``, ``ref``, ``alt``, and
    ``methylation_level`` from the annotated context Table (inner join on ``locus``
    and ``alleles``) and kept if the VEP most severe consequence is synonymous.
    Duplicate DNMs (same locus and alleles) across or within cohorts are collapsed
    into one row, with the reporting cohorts in ``cohorts``.

    :param cohort_hts: Dictionary of cohort name to Table returned by
        :func:`format_dnm_cohort`.
    :param context_ht: Annotated context Table from the constraint pipeline
        (``--prepare-context-ht`` output).
    :return: Tuple of the deduplicated synonymous DNM Table (keyed by ``locus`` and
        ``alleles``) and a dictionary with the number of synonymous DNMs per cohort
        before deduplication.
    """
    hts = [ht.select("cohort") for ht in cohort_hts.values()]
    ht = hts[0].union(*hts[1:])
    ht = ht.checkpoint(new_temp_file("syn_dnm.all_cohorts", "ht"))

    # The cohort Tables have few partitions (text files are imported as one), so
    # joining the context Table into them would scan most of the context Table in
    # each DNM partition. Instead, filter the context Table to the DNM sites in a
    # scan over its own partitions, then join the small filtered Table.
    context_ht = context_ht.semi_join(ht.select().distinct())
    context_ht = context_ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        csq=context_ht.vep.most_severe_consequence,
    )
    context_ht = context_ht.naive_coalesce(50)
    context_ht = context_ht.checkpoint(new_temp_file("syn_dnm.context", "ht"))
    ht = ht.annotate(**context_ht[ht.key])
    ht = ht.filter(ht.csq == SYNONYMOUS_CSQ).drop("csq")
    ht = ht.checkpoint(new_temp_file("syn_dnm.by_cohort", "ht"))
    cohort_counts = dict(ht.aggregate(hl.agg.counter(ht.cohort)))

    ht = ht.group_by("locus", "alleles").aggregate(
        cohorts=hl.agg.collect_as_set(ht.cohort),
        context=hl.agg.take(ht.context, 1)[0],
        ref=hl.agg.take(ht.ref, 1)[0],
        alt=hl.agg.take(ht.alt, 1)[0],
        methylation_level=hl.agg.take(ht.methylation_level, 1)[0],
    )

    return ht, cohort_counts


# ---------------------------------------------------------------------------
# Regional genomic features
# ---------------------------------------------------------------------------


def format_genomic_features(ht: hl.Table) -> hl.Table:
    """
    Key an imported 1kb regional genomic feature matrix by window interval.

    The matrix has an ``element_id`` column ("chr-start-end", BED coordinates)
    followed by ``GENOMIC_FEATURE_COLUMNS``. NaN feature values become missing.

    As in the Gnocchi rr Tables, each window's interval is ``[start, end)`` in
    1-based positions, which is 1bp left of the BED window, and windows starting at 0
    are removed.

    :param ht: Imported feature matrix with ``element_id`` and the feature columns
        as floats.
    :return: Feature Table keyed by ``interval`` with ``element_id`` and the feature
        columns.
    """
    ht = ht.select("element_id", *GENOMIC_FEATURE_COLUMNS)
    ht = ht.annotate(
        **{
            c: hl.if_else(hl.is_nan(ht[c]), hl.missing(hl.tfloat64), ht[c])
            for c in GENOMIC_FEATURE_COLUMNS
        }
    )
    coords = ht.element_id.split("-")
    ht = ht.annotate(
        _contig=coords[0], _start=hl.int(coords[1]), _end=hl.int(coords[2])
    )
    ht = ht.filter(ht._start > 0)
    ht = ht.key_by(
        interval=hl.locus_interval(
            ht._contig, ht._start, ht._end, reference_genome="GRCh38"
        )
    )

    return ht.drop("_contig", "_start", "_end")


def summarize_genomic_features(ht: hl.Table) -> Dict[str, Any]:
    """
    Summarize a genomic feature Table per contig.

    Returns the number of windows and the features that are entirely missing per
    contig, and the number of duplicated window intervals.

    :param ht: Feature Table returned by :func:`format_genomic_features`.
    :return: Dictionary with ``n_duplicated_windows`` and ``contigs`` (contig ->
        dict with ``n_windows`` and ``missing_features``).
    """
    n_dup = ht.count() - ht.select().distinct().count()

    contig_summary = ht.aggregate(
        hl.agg.group_by(
            ht.interval.start.contig,
            hl.struct(
                n_windows=hl.agg.count(),
                n_missing=hl.struct(
                    **{
                        c: hl.agg.count_where(hl.is_missing(ht[c]))
                        for c in GENOMIC_FEATURE_COLUMNS
                    }
                ),
            ),
        )
    )
    contigs = {}
    for contig, s in contig_summary.items():
        missing = sorted(
            {
                c.rsplit("_", 1)[0]
                for c in GENOMIC_FEATURE_COLUMNS
                if s.n_missing[c] == s.n_windows
            }
        )
        contigs[contig] = {"n_windows": s.n_windows, "missing_features": missing}

    return {"n_duplicated_windows": n_dup, "contigs": contigs}


# ---------------------------------------------------------------------------
# rr training set
# ---------------------------------------------------------------------------


def get_syn_sites_ht(context_ht: hl.Table) -> hl.Table:
    """
    Get all possible synonymous SNVs from the annotated context Table.

    :param context_ht: Annotated context Table from the constraint pipeline.
    :return: Table keyed by ``locus`` and ``alleles`` with ``context``, ``ref``,
        ``alt``, and ``methylation_level``.
    """
    ht = context_ht.filter(context_ht.vep.most_severe_consequence == SYNONYMOUS_CSQ)
    return ht.select("context", "ref", "alt", "methylation_level")


def create_rr_training_ht(
    syn_sites_ht: hl.Table,
    syn_dnm_ht: hl.Table,
    exomes_sites_ht: hl.Table,
    features_ht: hl.Table,
    background_fraction: float,
    seed: int,
) -> hl.Table:
    """
    Create the rr training Table of synonymous DNM and non-mutated background sites.

    The background is a random sample (``background_fraction``) of the possible
    synonymous sites that are not a DNM and are absent from the gnomAD exomes sites
    Table, with one site kept per locus. DNM sites get ``y`` = 1 and background
    sites ``y`` = 0. Each site is annotated with the features of the 1kb window
    containing it. ``row_idx`` gives the locus order within each class, which is the
    row order used to train the models.

    :param syn_sites_ht: Possible synonymous SNVs (:func:`get_syn_sites_ht`).
    :param syn_dnm_ht: Synonymous DNM Table (:func:`combine_syn_dnms`).
    :param exomes_sites_ht: gnomAD exomes sites Table.
    :param features_ht: Feature Table (:func:`format_genomic_features`).
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
    :return: Training Table keyed by ``locus`` with ``alleles``, ``context``, ``y``,
        ``row_idx``, and the feature columns.
    """
    bg_ht = syn_sites_ht.anti_join(syn_dnm_ht)
    bg_ht = bg_ht.anti_join(exomes_sites_ht.select().select_globals())
    bg_ht = bg_ht.key_by("locus").distinct()
    bg_ht = bg_ht.sample(background_fraction, seed=seed)
    bg_ht = bg_ht.key_by("locus").select("alleles", "context", y=0)
    bg_ht = bg_ht.add_index("row_idx")

    dnm_ht = syn_dnm_ht.key_by("locus").select("alleles", "context", y=1)
    dnm_ht = dnm_ht.add_index("row_idx")

    ht = dnm_ht.union(bg_ht)
    features_ht = features_ht.select(*GENOMIC_FEATURE_COLUMNS)

    return ht.annotate(**features_ht[ht.locus])


def remove_qc_flagged_genes(ht: hl.Table, metrics_ht: hl.Table) -> hl.Table:
    """
    Remove sites in genes with a release QC gene flag from a locus-keyed Table.

    Genes flagged ``low_exome_coverage`` or ``low_exome_mapping_quality`` lose
    variant calls in both the DNM cohorts and gnomAD, so leaving their sites in the
    rr training set can teach the models that those regions mutate less. The span
    (``start_position`` to ``end_position``) of every Ensembl MANE Select transcript
    with a non-empty ``gene_flags`` is removed.

    :param ht: Table keyed by ``locus`` (e.g. the rr training Table).
    :param metrics_ht: Constraint metrics Table with ``transcript``, ``mane_select``,
        ``gene_flags``, ``chromosome``, ``start_position`` and ``end_position``.
    :return: ``ht`` without the sites in flagged genes.
    """
    flagged = metrics_ht.filter(
        metrics_ht.mane_select
        & metrics_ht.transcript.startswith(ENSEMBL_TRANSCRIPT_PREFIX)
        & (hl.len(metrics_ht.gene_flags) > 0)
        & hl.is_defined(metrics_ht.start_position)
        & hl.is_defined(metrics_ht.end_position)
    )
    spans = flagged.key_by().select("chromosome", "start_position", "end_position")
    intervals = [
        hl.Interval(
            hl.Locus(r.chromosome, r.start_position, reference_genome="GRCh38"),
            hl.Locus(r.chromosome, r.end_position, reference_genome="GRCh38"),
            includes_end=True,
        )
        for r in spans.collect()
    ]
    logger.info("Removing sites in %d QC-flagged genes.", len(intervals))

    return hl.filter_intervals(ht, intervals, keep=False)


def rr_training_ht_to_pandas(ht: hl.Table) -> pd.DataFrame:
    """
    Collect the rr training Table into a pandas DataFrame in training row order.

    Rows are ordered by DNM sites first, then background sites, each in ``row_idx``
    order.

    :param ht: rr training Table (:func:`create_rr_training_ht`).
    :return: DataFrame with ``locus`` (as "contig:position"), ``context``, ``y``,
        ``row_idx``, and the feature columns.
    """
    ht = ht.key_by()
    ht = ht.select(
        "context", "y", "row_idx", *GENOMIC_FEATURE_COLUMNS, locus=hl.str(ht.locus)
    )
    df = ht.to_pandas()
    # to_pandas returns nullable extension dtypes (e.g. Float64); use NumPy dtypes
    # (missing values as NaN).
    df = df.astype(
        {
            "y": "int64",
            "row_idx": "int64",
            **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, "float64"),
        }
    )
    df = df.sort_values(["y", "row_idx"], ascending=[False, True], kind="stable")

    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# rr models
# ---------------------------------------------------------------------------


def get_rr_model_features(context: str) -> List[str]:
    """
    Get the ordered model features of a context.

    Every feature goes in at every scale, ordered by scale and then by
    ``RR_FEATURE_ORDER``. Methylation-correlated features are dropped for CpG
    contexts, since the mutation rates already model methylation there.

    :param context: Trinucleotide context.
    :return: Feature column names (``{feature}_{scale}``).
    """
    return [
        f"{f}_{s}"
        for s in GENOMIC_FEATURE_SCALES
        for f in RR_FEATURE_ORDER
        if not (context in CPG_CONTEXTS and f in METHYLATION_CORRELATED_FEATURES)
    ]


def fit_rr_context_model(
    df: pd.DataFrame,
    features: List[str],
    ridge_penalty: float = RR_RIDGE_PENALTY,
) -> Optional[Dict[str, Any]]:
    """
    Fit the ridge-penalized rr logistic model of one context.

    Sites with any missing feature are dropped, DNM rows are followed by background
    rows (training row order), the features are z-scored (ddof=0) and rotated with a
    full-rank ``IncrementalPCA``, and a logistic regression is fit that minimizes
    mean log loss + ``ridge_penalty`` / 2 * sum(coef ** 2), with the intercept
    unpenalized. The PCA rotation is orthogonal, so this is the same ridge fit as on
    the z-scored features.

    The PCA is folded into the model coefficients so that the linear predictor of a
    window is ``intercept + sum(coefs * (x - mean) / sd)``. As in the Gnocchi models,
    ``sd`` is the ddof=1 standard deviation (applied at prediction) while the model
    is trained on ddof=0 z-scores. ``ave`` is the predicted probability of the
    average window (all principal components equal to 0).

    :param df: Training rows of one context (:func:`rr_training_ht_to_pandas`).
    :param features: Model features (:func:`get_rr_model_features`).
    :param ridge_penalty: Per-site ridge penalty. Default is ``RR_RIDGE_PENALTY``.
    :return: Dictionary of model parameters, or None if there is no training data.
    """
    cols = ["locus", *features]
    df_1 = df[df["y"] == 1][cols].dropna()
    df_0 = df[df["y"] == 0][cols].dropna()
    if df_1.empty or df_0.empty:
        return None
    df_01 = pd.concat([df_1.assign(group=1), df_0.assign(group=0)])
    x = df_01[features]
    y = df_01["group"]

    x_scaled = x.apply(scipy.stats.zscore)
    pca = IncrementalPCA(n_components=len(features))
    x_pca = pca.fit_transform(x_scaled)
    # sklearn minimizes C * sum(log loss) + sum(coef ** 2) / 2, so C = 1 / (n *
    # penalty) gives the per-site penalty.
    fit = LogisticRegression(C=1.0 / (len(y) * ridge_penalty), max_iter=10000)
    fit.fit(x_pca, y)
    params = np.concatenate([fit.intercept_, fit.coef_[0]])

    coefs = pca.components_.T @ params[1:]
    intercept = params[0] - pca.mean_ @ coefs

    return {
        "features": list(features),
        "mean": [float(v) for v in x.mean()],
        "sd": [float(v) for v in x.std()],
        "coefs": [float(v) for v in coefs],
        "intercept": float(intercept),
        "ave": float(1.0 / (1.0 + np.exp(-params[0]))),
        "n_dnm": int(len(df_1)),
        "n_background": int(len(df_0)),
        "converged": bool(fit.n_iter_[0] < fit.max_iter),
    }


def fit_rr_models(
    df: pd.DataFrame,
    ridge_penalty: float = RR_RIDGE_PENALTY,
    contexts: List[str] = RR_CONTEXTS,
) -> List[Dict[str, Any]]:
    """
    Fit the rr model of every context.

    :param df: Training DataFrame (:func:`rr_training_ht_to_pandas`).
    :param ridge_penalty: Per-site ridge penalty. Default is ``RR_RIDGE_PENALTY``.
    :param contexts: Contexts to model. Default is ``RR_CONTEXTS``.
    :return: List of model dictionaries (:func:`fit_rr_context_model`) with a
        ``context`` field.
    """
    models = []
    for context in contexts:
        model = fit_rr_context_model(
            df[df["context"] == context], get_rr_model_features(context), ridge_penalty
        )
        if model is not None:
            models.append({"context": context, **model})
    logger.info("%d of %d contexts modeled.", len(models), len(contexts))

    return models


def rr_models_to_ht(models: List[Dict[str, Any]]) -> hl.Table:
    """
    Convert rr model dictionaries to a Hail Table keyed by context.

    :param models: Model dictionaries (:func:`fit_rr_models`).
    :return: rr model Table.
    """
    schema = hl.tstruct(
        context=hl.tstr,
        features=hl.tarray(hl.tstr),
        mean=hl.tarray(hl.tfloat64),
        sd=hl.tarray(hl.tfloat64),
        coefs=hl.tarray(hl.tfloat64),
        intercept=hl.tfloat64,
        ave=hl.tfloat64,
        n_dnm=hl.tint32,
        n_background=hl.tint32,
        converged=hl.tbool,
    )

    return hl.Table.parallelize(models, schema=schema).key_by("context")


# ---------------------------------------------------------------------------
# rr per window and context
# ---------------------------------------------------------------------------


def _rr_model_expr(ht: hl.Table, model: hl.Struct) -> hl.expr.Float64Expression:
    """
    Get the rr expression of one context model for each window of a feature Table.

    rr = sigmoid(intercept + sum(coefs * (x - mean) / sd)) / ave, and is missing when
    any model feature is missing for the window.

    :param ht: Feature Table keyed by ``interval``.
    :param model: Model row of the rr model Table.
    :return: rr expression.
    """
    eta = hl.float64(model.intercept)
    for feature, mean, sd, coef in zip(
        model.features, model.mean, model.sd, model.coefs
    ):
        z = (ht[feature] - mean) / sd
        eta = eta + coef * z

    return (1.0 / (1.0 + hl.exp(-eta))) / model.ave


def compute_rr_ht(features_ht: hl.Table, models: List[hl.Struct]) -> hl.Table:
    """
    Compute rr per 1kb window and trinucleotide context.

    Contexts without a model, and windows missing a model feature for a context, get
    an rr of 1.0. Windows are kept only if at least one modeled context has an rr for
    the window.

    :param features_ht: Feature Table (:func:`format_genomic_features`).
    :param models: Rows of the rr model Table.
    :return: rr Table keyed by ``interval`` with an ``adj_r`` dict.
    """
    rr_exprs = {m.context: _rr_model_expr(features_ht, m) for m in models}
    if not rr_exprs:
        logger.warning("No context has an rr model, the rr Table will be empty.")

    ht = features_ht.select(_rr=hl.struct(**rr_exprs))
    ht = ht.filter(hl.any([hl.is_defined(ht._rr[c]) for c in rr_exprs]))

    return ht.select(
        adj_r=hl.dict(
            [
                (c, hl.or_else(ht._rr[c], 1.0) if c in rr_exprs else hl.float64(1.0))
                for c in RR_CONTEXTS
            ]
        )
    )


def get_syn_mane_select_ht(per_variant_ht: hl.Table) -> hl.Table:
    """
    Get the synonymous Ensembl MANE Select rows of the per-variant expected Table.

    Only Ensembl transcripts are kept: the RefSeq MANE Select transcripts are the
    same transcripts, so keeping both would count every site twice. Keeps the
    ``context`` and adj ``exp`` used to normalize rr.

    :param per_variant_ht: Per-variant expected Table from the constraint pipeline.
    :return: Table keyed by ``locus`` and ``alleles`` with ``context`` and ``exp``.
    """
    adj_idx = hl.eval(per_variant_ht.exomes_freq_meta).index(ADJ_FREQ_META)
    ht = per_variant_ht.filter(
        (per_variant_ht.annotation == SYNONYMOUS_CSQ)
        & per_variant_ht.mane_select
        & per_variant_ht.transcript.startswith(ENSEMBL_TRANSCRIPT_PREFIX)
    )
    ht = ht.select("context", exp=ht.expected_variants[adj_idx]).select_globals()

    return ht.naive_coalesce(200)


def compute_rr_normalization_factors(
    syn_ht: hl.Table, rr_ht: hl.Table
) -> Dict[str, float]:
    """
    Compute the per-context rr normalization factors on exome synonymous sites.

    The factor of a context is the expected-weighted mean rr over synonymous MANE
    Select sites of the per-variant expected Table that fall in an rr window:
    ``sum(exp * rr) / sum(exp)``.

    :param syn_ht: Synonymous MANE Select per-variant rows
        (:func:`get_syn_mane_select_ht`).
    :param rr_ht: rr Table (:func:`compute_rr_ht`).
    :return: Dictionary of context to normalization factor.
    """
    rr_expr = rr_ht[syn_ht.locus].adj_r.get(syn_ht.context)
    sums = syn_ht.aggregate(
        hl.agg.filter(
            hl.is_defined(rr_expr) & hl.is_defined(syn_ht.exp),
            hl.agg.group_by(
                syn_ht.context,
                hl.struct(
                    num=hl.agg.sum(syn_ht.exp * rr_expr), den=hl.agg.sum(syn_ht.exp)
                ),
            ),
        )
    )

    return {c: s.num / s.den for c, s in sums.items() if s.den > 0}


def normalize_rr_ht(rr_ht: hl.Table, factors: Dict[str, float]) -> hl.Table:
    """
    Divide rr by its per-context normalization factor in every window.

    Contexts without a factor are left unchanged.

    :param rr_ht: rr Table.
    :param factors: Dictionary of context to normalization factor
        (:func:`compute_rr_normalization_factors`).
    :return: Normalized rr Table.
    """
    factors_expr = hl.literal(factors, hl.tdict(hl.tstr, hl.tfloat64))
    return rr_ht.annotate(
        adj_r=hl.dict(
            rr_ht.adj_r.items().map(
                lambda kv: (kv[0], kv[1] / factors_expr.get(kv[0], 1.0))
            )
        )
    )
