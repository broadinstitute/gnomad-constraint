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

The DNMs and genomic features are imported by ``import_data.py`` (see
:mod:`gnomad_constraint.utils.import_data`). rr is applied to the per-variant
expected counts after model calibration; see
:func:`gnomad_constraint.utils.constraint.aggregate_per_variant_expected_ht`.
"""

import logging
from typing import Any, Dict, List, Optional

import hail as hl
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    CPG_CONTEXTS,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_SCALES,
    GENOMIC_FEATURES,
    METHYLATION_CORRELATED_FEATURES,
    RR_CONTEXTS,
    RR_FOLDS,
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
    dnm_ht: hl.Table,
    exomes_sites_ht: hl.Table,
    features_ht: hl.Table,
    background_fraction: float,
    seed: int,
) -> hl.Table:
    """
    Create the rr training Table of synonymous DNM and non-mutated background sites.

    The DNM sites are the synonymous sites of the DNM Table. The background is a
    random sample (``background_fraction``) of the possible synonymous sites that are
    not a DNM and are absent from the gnomAD exomes sites Table, with one site kept
    per locus. DNM sites get ``y`` = 1 and background sites ``y`` = 0. Each site is
    annotated with the features of the 1kb window containing it.

    :param syn_sites_ht: Possible synonymous SNVs (:func:`get_syn_sites_ht`).
    :param dnm_ht: DNM Table
        (:func:`~gnomad_constraint.utils.import_data.combine_dnms`).
    :param exomes_sites_ht: gnomAD exomes sites Table.
    :param features_ht: Feature Table
        (:func:`~gnomad_constraint.utils.import_data.format_genomic_features`).
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
    :return: Training Table keyed by ``locus`` with ``alleles``, ``context``, ``y``,
        and the feature columns.
    """
    syn_dnm_ht = dnm_ht.filter(dnm_ht.most_severe_consequence == SYNONYMOUS_CSQ)
    bg_ht = syn_sites_ht.anti_join(syn_dnm_ht)
    # TODO: Whether a site is absent from gnomAD depends on its own mutation rate
    #  (about half of CpG synonymous sites are seen in gnomAD exomes) and on coverage,
    #  so this background makes the models fit the DNM rate ratio times the ratio of
    #  absence probabilities, which steepens the regional effect, most at CpG. It
    #  partly offsets the compression of the sigmoid form of rr (see _rr_model_expr).
    #  Consider sampling from all non-DNM synonymous sites; the ridge penalty would
    #  need retuning.
    bg_ht = bg_ht.anti_join(exomes_sites_ht.select().select_globals())
    bg_ht = bg_ht.key_by("locus").distinct()
    bg_ht = bg_ht.sample(background_fraction, seed=seed)
    bg_ht = bg_ht.key_by("locus").select("alleles", "context", y=0)

    syn_dnm_ht = syn_dnm_ht.key_by("locus").select("alleles", "context", y=1)

    ht = syn_dnm_ht.union(bg_ht)
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
    Collect the rr training Table into a pandas DataFrame.

    :param ht: rr training Table (:func:`create_rr_training_ht`).
    :return: DataFrame with ``contig``, ``context``, ``y``, and the feature columns.
    """
    ht = ht.key_by()
    df = ht.select(
        "context", "y", *GENOMIC_FEATURE_COLUMNS, contig=ht.locus.contig
    ).to_pandas()
    # to_pandas returns nullable extension dtypes (e.g. Float64); use NumPy dtypes
    # (missing values as NaN).
    return df.astype(
        {"y": "int64", **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, "float64")}
    )


# ---------------------------------------------------------------------------
# rr models
# ---------------------------------------------------------------------------


def get_rr_model_features(context: str) -> List[str]:
    """
    Get the model features of a context.

    Every feature goes in at every scale. Methylation-correlated features are dropped
    for CpG contexts, since the mutation rates already model methylation there.

    :param context: Trinucleotide context.
    :return: Feature column names (``{feature}_{scale}``).
    """
    return [
        f"{f}_{s}"
        for s in GENOMIC_FEATURE_SCALES
        for f in GENOMIC_FEATURES
        if not (context in CPG_CONTEXTS and f in METHYLATION_CORRELATED_FEATURES)
    ]


def fit_rr_context_model(
    df: pd.DataFrame,
    features: List[str],
    ridge_penalty: float,
) -> Optional[Dict[str, Any]]:
    """
    Fit the ridge-penalized rr logistic model of one context.

    Sites with any missing feature are dropped, and features with a single value in
    the remaining sites are left out of the model, since they carry no information
    and can't be standardized. The features are z-scored and a logistic regression
    is fit that minimizes mean log loss + ``ridge_penalty`` / 2 * sum(coef ** 2),
    with the intercept unpenalized. The linear predictor of a window is
    ``intercept + sum(coefs * (x - mean) / sd)``, and ``ave`` is the predicted
    probability of the average window (every feature at its mean).

    :param df: Training rows of one context (:func:`rr_training_ht_to_pandas`).
    :param features: Model features (:func:`get_rr_model_features`).
    :param ridge_penalty: Per-site ridge penalty (:func:`choose_rr_penalty`).
    :return: Dictionary of model parameters, or None if the context has no DNM or
        no background sites, or no feature varies.
    """
    df = df[["y", *features]].dropna()
    n_dnm = int((df["y"] == 1).sum())
    n_background = int((df["y"] == 0).sum())
    mean = df[features].mean()
    sd = df[features].std(ddof=0)
    features = [f for f in features if sd[f] > 0]
    if not n_dnm or not n_background or not features:
        return None

    z = (df[features] - mean[features]) / sd[features]
    # sklearn minimizes C * sum(log loss) + sum(coef ** 2) / 2, so C = 1 / (n *
    # penalty) gives the per-site penalty.
    fit = LogisticRegression(C=1.0 / (len(df) * ridge_penalty), max_iter=10000)
    fit.fit(z.to_numpy(), df["y"].to_numpy())
    intercept = float(fit.intercept_[0])

    return {
        "features": features,
        "mean": [float(mean[f]) for f in features],
        "sd": [float(sd[f]) for f in features],
        "coefs": [float(v) for v in fit.coef_[0]],
        "intercept": intercept,
        "ave": float(1.0 / (1.0 + np.exp(-intercept))),
        "n_dnm": n_dnm,
        "n_background": n_background,
        "converged": bool(fit.n_iter_[0] < fit.max_iter),
    }


def fit_rr_models(
    df: pd.DataFrame,
    ridge_penalty: float,
    contexts: List[str] = RR_CONTEXTS,
) -> List[Dict[str, Any]]:
    """
    Fit the rr model of every context.

    :param df: Training DataFrame (:func:`rr_training_ht_to_pandas`).
    :param ridge_penalty: Per-site ridge penalty (:func:`choose_rr_penalty`).
    :param contexts: Contexts to model. Default is ``RR_CONTEXTS``.
    :return: List of model dictionaries (:func:`fit_rr_context_model`) with a
        ``context`` field.
    """
    models = []
    for context in contexts:
        features = get_rr_model_features(context)
        model = fit_rr_context_model(
            df[df["context"] == context], features, ridge_penalty
        )
        if model is None:
            logger.info(
                "No rr model for %s: no DNM or background sites, or no feature varies.",
                context,
            )
            continue
        constant = [f for f in features if f not in model["features"]]
        if constant:
            logger.info("%s: features left out (constant): %s", context, constant)
        models.append({"context": context, **model})
    logger.info("%d of %d contexts modeled.", len(models), len(contexts))

    return models


def get_chrom_fold(contig: str) -> Optional[str]:
    """
    Get the chromosome-split fold of a contig (see ``RR_FOLDS``).

    :param contig: Contig name (e.g. "chr12").
    :return: "odd" or "even" for an autosome, otherwise None (the genomic features,
        and so the training sites, are autosomal).
    """
    name = contig.replace("chr", "")
    if not name.isdigit():
        return None

    return RR_FOLDS[0] if int(name) % 2 == 1 else RR_FOLDS[1]


def _predict_rr_model(model: Dict[str, Any], df: pd.DataFrame) -> np.ndarray:
    """
    Get the predicted DNM probability of each site from an rr model.

    :param model: Model dictionary (:func:`fit_rr_context_model`).
    :param df: Sites with the model's feature columns.
    :return: Predicted probabilities.
    """
    z = (df[model["features"]].to_numpy() - model["mean"]) / model["sd"]
    eta = model["intercept"] + z @ np.array(model["coefs"])

    return 1.0 / (1.0 + np.exp(-eta))


def tune_rr_penalty(
    df: pd.DataFrame,
    penalties: List[float],
    contexts: List[str] = RR_CONTEXTS,
) -> pd.DataFrame:
    """
    Get the held-out log loss of chromosome-split rr models for each ridge penalty.

    For each penalty, context, and fold in ``RR_FOLDS``, the context's model is fit
    on the training sites outside the fold's chromosomes and scored on the sites in
    them. Only DNMs and the training background are used, no gnomAD data.

    :param df: Training DataFrame (:func:`rr_training_ht_to_pandas`).
    :param penalties: Per-site ridge penalties to try.
    :param contexts: Contexts to model. Default is ``RR_CONTEXTS``.
    :return: DataFrame with ``penalty``, ``context``, ``fold``, the number of held-out
        sites (``n_test``), and their summed log loss (``log_loss``).
    """
    df = df.assign(_fold=df["contig"].map(get_chrom_fold))
    rows = []
    for penalty in penalties:
        for context in contexts:
            features = get_rr_model_features(context)
            ctx_df = df[df["context"] == context]
            for fold in RR_FOLDS:
                model = fit_rr_context_model(
                    ctx_df[ctx_df["_fold"] != fold], features, penalty
                )
                test_df = ctx_df[ctx_df["_fold"] == fold][["y", *features]].dropna()
                if model is None or test_df.empty:
                    continue
                rows.append(
                    {
                        "penalty": float(penalty),
                        "context": context,
                        "fold": fold,
                        "n_test": len(test_df),
                        "log_loss": float(
                            log_loss(
                                test_df["y"],
                                _predict_rr_model(model, test_df),
                                labels=[0, 1],
                                normalize=False,
                            )
                        ),
                    }
                )

    return pd.DataFrame(
        rows, columns=["penalty", "context", "fold", "n_test", "log_loss"]
    )


def choose_rr_penalty(
    tuning_df: pd.DataFrame, default: float = RR_RIDGE_PENALTY
) -> float:
    """
    Choose the ridge penalty with the lowest held-out log loss.

    The log loss is summed over contexts and folds, so one penalty is used for every
    context. If no held-out sites could be scored (the training sites are all on one
    side of the chromosome split, as in a ``--test`` run), ``default`` is returned.

    :param tuning_df: Held-out log loss per penalty (:func:`tune_rr_penalty`).
    :param default: Penalty used when there is no held-out log loss. Default is
        ``RR_RIDGE_PENALTY``.
    :return: Chosen per-site ridge penalty.
    """
    if tuning_df.empty:
        logger.warning(
            "No held-out sites to score the ridge penalties on; using %s.", default
        )
        return default

    return float(tuning_df.groupby("penalty")["log_loss"].sum().idxmin())


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

    rr = sigmoid(intercept + sum(coefs * (x - mean) / sd)) / ave, as in Gnocchi, and
    is missing when any model feature is missing for the window.

    Because the background is a sample of non-mutated sites, the model is
    case-control, and the ratio of mutation rates between a window and the average
    window is closer to the odds ratio ``exp(eta - intercept)`` than to this ratio of
    probabilities. The sigmoid form compresses the regional effect toward 1 by about
    ``1 - ave`` on the log scale and caps rr at ``1 / ave``. That is negligible for
    non-CpG contexts, where DNMs are a few percent of the training sites, but
    damps CpG contexts, where they are 20-50%, and the damping changes with
    ``RR_BACKGROUND_FRACTION``. The sigmoid form is kept because it is the form the
    rr was evaluated with.

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

    :param features_ht: Feature Table
        (:func:`~gnomad_constraint.utils.import_data.format_genomic_features`).
    :param models: Rows of the rr model Table.
    :return: rr Table keyed by ``interval`` with an ``adj_r`` dict.
    """
    rr_exprs = {m.context: _rr_model_expr(features_ht, m) for m in models}
    if not rr_exprs:
        logger.warning("No context has an rr model, the rr Table will be empty.")

    ht = features_ht.select(_rr=hl.struct(**rr_exprs))
    ht = ht.filter(hl.any([hl.is_defined(ht._rr[c]) for c in rr_exprs]))

    # TODO: Windows missing a context's model features get an rr of 1.0 here, which is
    #  the model's rr for the average-feature window. normalize_rr_ht then divides it
    #  by the context's factor, so these windows end up at 1 / factor (on average about
    #  1.4% below 1), unlike every other site without an rr (chrX/Y, sites outside a
    #  window, and windows with no rr for any context). Consider leaving them missing
    #  here, computing the factors over predicted sites only, and filling 1.0 after
    #  normalization.
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
    return ht.select("context", exp=ht.expected_variants[adj_idx]).select_globals()


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
