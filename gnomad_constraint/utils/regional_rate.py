"""
Utility functions for the regional mutation rate (rr) adjustment of expected counts.

The rr adjustment is a port of the Gnocchi regional adjustment
(https://github.com/atgu/gnomad_nc_constraint) retrained on synonymous de novo
mutations (DNMs):

    - Univariate logistic screen (DNM vs non-mutated synonymous background sites) of
      13 regional genomic features at 4 flanking scales, per strand-collapsed
      trinucleotide context, to select features.
    - Per-context logistic model on the selected (z-scored, PCA-rotated) features.
    - Per 1kb window and context: rr = P(window features) / P(average window).

rr is applied to the per-variant expected counts after model calibration
(``expected_variants * rr``); see
:func:`gnomad_constraint.utils.constraint.annotate_regional_rr_expected`.
"""

import json
import logging
import os
import urllib.request
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.api as sm
from gnomad.resources.resource_utils import TableResource
from gnomad.utils.constraint import calculate_raw_z_score
from gnomad.utils.liftover import default_lift_data
from hail.utils.misc import new_temp_file
from sklearn.decomposition import IncrementalPCA
from sklearn.metrics import log_loss, roc_auc_score

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    CPG_CONTEXTS,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_COVERAGE_TRACKS,
    GENOMIC_FEATURE_DISTANCE_TRACKS,
    GENOMIC_FEATURE_MID_BEDS,
    GENOMIC_FEATURE_SCALES,
    GENOMIC_FEATURE_WINDOW_BEDS,
    GENOMIC_FEATURES,
    METHYLATION_CORRELATED_FEATURES,
    RR_BONFERRONI_CUTOFF,
    RR_CONTEXTS,
    RR_FEATURE_SELECTION_ALPHA,
    RR_FOLDS,
    RR_SCREEN_FEATURE_ORDER,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_utils")
logger.setLevel(logging.INFO)

SYNONYMOUS_CSQ = "synonymous_variant"
"""VEP most severe consequence used to define synonymous sites and DNMs."""

UCSC_SEQUENCE_API = "https://api.genome.ucsc.edu/getData/sequence"
"""UCSC REST API endpoint used to fetch hg38 reference sequence for GC content."""

MODEL_FOLD_ALL = "all"
"""Fold label of rr models trained on all training sites."""


# ---------------------------------------------------------------------------
# Synonymous de novo mutations
# ---------------------------------------------------------------------------


def _import_dnm_text(path: str) -> hl.Table:
    """
    Import a DNM cohort text file.

    :param path: Path to the tab-delimited DNM file.
    :return: Imported Table.
    """
    return hl.import_table(path, impute=True)


def _lift_to_grch38(ht: hl.Table) -> hl.Table:
    """
    Lift a GRCh37 Table keyed by ``locus`` and ``alleles`` over to GRCh38.

    Sites that fail liftover are removed.

    :param ht: GRCh37 Table with ``locus`` and ``alleles`` fields.
    :return: GRCh38 Table keyed by ``locus`` and ``alleles`` with no other fields.
    """
    return default_lift_data(ht.key_by("locus", "alleles")).select()


def load_dnm_cohort(cohort: str, resource: Union[str, TableResource]) -> hl.Table:
    """
    Load one DNM cohort as a GRCh38 Table keyed by ``locus`` and ``alleles``.

    Supported cohorts (see
    :func:`gnomad_constraint.resources.resource_utils.get_dnm_cohort_resources`):

        - ``decode_psychencode``: deCODE + PsychENCODE DNMs used by Gnocchi
          (GRCh38). Each DNM is listed in both allele orientations; the wrong
          orientation is removed when joining to the context Table.
        - ``ddd_asd``: DDD + ASD DNMs (GRCh38 Hail Table).
        - ``gnomad_v4.1``: gnomAD v4.1 high-quality coding de novos (GRCh38 Hail
          Table).
        - ``chd_2017``, ``ee_2014``, ``cdh_2021``: GRCh37 text files, lifted over.
        - ``ofc_2020``: text file with GRCh38 positions.

    :param cohort: Cohort name.
    :param resource: Path (text cohorts) or TableResource (Hail Table cohorts).
    :return: Table keyed by ``locus`` and ``alleles`` with a ``cohort`` field.
    """
    if cohort in {"ddd_asd", "gnomad_v4.1"}:
        ht = resource.ht().key_by("locus", "alleles").select()
    elif cohort == "decode_psychencode":
        ht = _import_dnm_text(resource)
        ht = ht.select(
            locus=hl.parse_locus(ht.locus, reference_genome="GRCh38"),
            alleles=hl.array([ht.ref, ht.alt]),
        ).key_by("locus", "alleles")
    elif cohort == "chd_2017":
        ht = _import_dnm_text(resource)
        ht = ht.select(
            locus=hl.parse_locus(
                ht.CHROM + ":" + hl.str(ht.POS), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.REF, ht.ALT]),
        )
        ht = _lift_to_grch38(ht)
    elif cohort == "ofc_2020":
        ht = _import_dnm_text(resource)
        ht = ht.select(
            locus=hl.parse_locus(
                "chr" + hl.str(ht.Chromosome) + ":" + hl.str(ht["Position (build 38)"]),
                reference_genome="GRCh38",
            ),
            alleles=hl.array([ht.Ref, ht.Alt]),
        ).key_by("locus", "alleles")
    elif cohort == "ee_2014":
        ht = _import_dnm_text(resource)
        ht = ht.select(
            locus=hl.parse_locus(
                ht["hg19 coordinates (chr:position)"].replace("chr", ""),
                reference_genome="GRCh37",
            ),
            alleles=ht["Ref/Alt alleles"].split("/"),
        )
        ht = _lift_to_grch38(ht)
    elif cohort == "cdh_2021":
        ht = _import_dnm_text(resource)
        ht = ht.select(
            locus=hl.parse_locus(
                hl.str(ht.chr) + ":" + hl.str(ht.position), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.ref, ht.alt]),
        )
        ht = _lift_to_grch38(ht)
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
    ``methylation_level`` from the annotated context Table (inner join on
    ``locus`` and ``alleles``) and kept if the VEP most severe consequence is
    synonymous. DNMs are kept on all contigs. Duplicate DNMs (same locus and
    alleles) across or within cohorts are always collapsed into one row, with the
    reporting cohorts in ``cohorts``.

    :param cohort_hts: Dictionary of cohort name to Table returned by
        :func:`load_dnm_cohort`.
    :param context_ht: Annotated context Table from the constraint pipeline
        (``--prepare-context-ht`` output).
    :return: Tuple of the deduplicated synonymous DNM Table (keyed by ``locus`` and
        ``alleles``) and a dictionary with the number of synonymous DNMs per cohort
        before deduplication.
    """
    hts = [ht.select("cohort") for ht in cohort_hts.values()]
    ht = hts[0].union(*hts[1:])

    context_ht = context_ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        csq=context_ht.vep.most_severe_consequence,
    )
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


def window_interval_expr(
    contig: hl.expr.StringExpression,
    start: hl.expr.Int32Expression,
    end: hl.expr.Int32Expression,
    corrected_window_coords: bool = False,
) -> hl.expr.IntervalExpression:
    """
    Get the locus interval of a 1kb window from its BED (0-based, half-open) coordinates.

    The legacy coordinates (used to build the existing ``syn_adj_r`` Tables and
    required to reproduce them) are ``[start, end)`` in 1-based positions, which is
    1bp left of the BED window. The corrected coordinates are ``[start + 1, end]``,
    exactly the BED window.

    :param contig: Window contig.
    :param start: Window BED start.
    :param end: Window BED end.
    :param corrected_window_coords: Whether to use the corrected coordinates.
        Default is False.
    :return: Locus interval expression.
    """
    if corrected_window_coords:
        return hl.locus_interval(
            contig,
            start + 1,
            end,
            includes_start=True,
            includes_end=True,
            reference_genome="GRCh38",
        )

    return hl.locus_interval(contig, start, end, reference_genome="GRCh38")


def import_genomic_features(
    paths: List[str],
    corrected_window_coords: bool = False,
    min_partitions: int = 100,
) -> hl.Table:
    """
    Import 1kb regional genomic feature matrices into a Table keyed by window interval.

    Each matrix has an ``element_id`` column ("chr-start-end", BED coordinates)
    followed by ``GENOMIC_FEATURE_COLUMNS``. "NA", "nan", "NaN", and empty values
    (and NaN floats) are imported as missing. With the legacy coordinates, windows
    starting at 0 are removed (as when the legacy Tables were built).

    :param paths: Paths of the feature matrices to import and union.
    :param corrected_window_coords: Whether to use the corrected window coordinates.
        See :func:`window_interval_expr`. Default is False.
    :param min_partitions: Minimum number of partitions when importing. Default is
        100.
    :return: Feature Table keyed by ``interval`` with ``element_id`` and the feature
        columns.
    """
    hts = []
    for path in paths:
        ht = hl.import_table(
            path,
            types={c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
            missing=["NA", "nan", "NaN", ""],
            min_partitions=min_partitions,
            force=path.endswith(".gz"),
        )
        ht = ht.select("element_id", *GENOMIC_FEATURE_COLUMNS)
        hts.append(ht)
    ht = hts[0].union(*hts[1:]) if len(hts) > 1 else hts[0]

    ht = ht.annotate(
        **{c: hl.or_missing(~hl.is_nan(ht[c]), ht[c]) for c in GENOMIC_FEATURE_COLUMNS}
    )
    coords = ht.element_id.split("-")
    ht = ht.annotate(
        _contig=coords[0], _start=hl.int(coords[1]), _end=hl.int(coords[2])
    )
    if not corrected_window_coords:
        ht = ht.filter(ht._start > 0)

    ht = ht.key_by(
        interval=window_interval_expr(
            ht._contig, ht._start, ht._end, corrected_window_coords
        )
    )

    return ht.drop("_contig", "_start", "_end")


def summarize_genomic_features(ht: hl.Table) -> Dict[str, Any]:
    """
    Summarize a genomic feature Table per contig.

    Returns the number of windows and the features that are entirely missing per
    contig, and the number of duplicated window intervals.

    :param ht: Feature Table returned by :func:`import_genomic_features`.
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
    synonymous sites that are not a DNM and are absent from the gnomAD exomes
    sites Table, with one site kept per locus. DNM sites get ``y`` = 1 and
    background sites ``y`` = 0. Each site is annotated with the features of the 1kb
    window containing it. ``row_idx`` gives the locus order within each class,
    which is the row order used to train the models.

    :param syn_sites_ht: Possible synonymous SNVs (:func:`get_syn_sites_ht`).
    :param syn_dnm_ht: Synonymous DNM Table (:func:`combine_syn_dnms`).
    :param exomes_sites_ht: gnomAD exomes sites Table.
    :param features_ht: Feature Table (:func:`import_genomic_features`).
    :param background_fraction: Fraction of non-mutated sites to sample.
    :param seed: Random seed used for sampling.
    :return: Training Table keyed by ``locus`` with ``alleles``, ``context``, ``y``,
        ``row_idx``, and the feature columns.
    """
    bg_ht = syn_sites_ht.anti_join(syn_dnm_ht)
    bg_ht = bg_ht.anti_join(exomes_sites_ht.select().select_globals())
    bg_ht = bg_ht.key_by("locus").distinct()
    bg_ht = bg_ht.sample(background_fraction, seed=seed)
    bg_ht = bg_ht.select("alleles", "context", y=0).add_index("row_idx")

    dnm_ht = syn_dnm_ht.key_by("locus").select("alleles", "context", y=1)
    dnm_ht = dnm_ht.add_index("row_idx")

    ht = dnm_ht.union(bg_ht)
    features_ht = features_ht.select(*GENOMIC_FEATURE_COLUMNS)

    return ht.annotate(**features_ht[ht.locus])


def rr_training_ht_to_pandas(ht: hl.Table) -> pd.DataFrame:
    """
    Collect the rr training Table into a pandas DataFrame in training row order.

    Rows are ordered by DNM sites first, then background sites, each in ``row_idx``
    order (the legacy training row order).

    :param ht: rr training Table (:func:`create_rr_training_ht`).
    :return: DataFrame with ``locus`` (as "contig:position"), ``context``, ``y``,
        ``row_idx``, and the feature columns.
    """
    ht = ht.key_by().select(
        "context", "y", "row_idx", *GENOMIC_FEATURE_COLUMNS, locus=hl.str(ht.locus)
    )
    df = ht.to_pandas()
    df = df.sort_values(["y", "row_idx"], ascending=[False, True], kind="stable")

    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# rr feature screen and models
# ---------------------------------------------------------------------------


def _chrom_fold(contig: str) -> str:
    """
    Get the chromosome-split fold of a contig.

    Autosomes are assigned by chromosome number parity, chrX to "odd", and chrY to
    "even" (see ``RR_FOLDS``).

    :param contig: Contig name (e.g. "chr12").
    :return: Fold name.
    """
    name = contig.replace("chr", "")
    if name == "X":
        return RR_FOLDS[0]
    if name == "Y":
        return RR_FOLDS[1]

    return RR_FOLDS[0] if int(name) % 2 == 1 else RR_FOLDS[1]


def _chrom_fold_expr(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    """
    Get the chromosome-split fold of a contig expression (see :func:`_chrom_fold`).

    :param contig: Contig expression.
    :return: Fold name expression.
    """
    name = contig.replace("chr", "")
    return (
        hl.case()
        .when(name == "X", RR_FOLDS[0])
        .when(name == "Y", RR_FOLDS[1])
        .when(hl.int(name) % 2 == 1, RR_FOLDS[0])
        .default(RR_FOLDS[1])
    )


def run_rr_feature_screen(
    df: pd.DataFrame, contexts: List[str] = RR_CONTEXTS
) -> pd.DataFrame:
    """
    Run the univariate logistic screen of each feature, scale, and context.

    For each context, scale, and feature (in ``RR_SCREEN_FEATURE_ORDER``), sites
    with a missing feature value are dropped, the feature is z-scored, and a
    logistic regression of ``y`` (DNM vs background) on the feature is fit with
    ``statsmodels`` ``Logit.fit_regularized`` (defaults, i.e. unpenalized). Fits
    that fail are recorded with NaN values (and are never selected).

    :param df: Training DataFrame (:func:`rr_training_ht_to_pandas`).
    :param contexts: Contexts to screen. Default is ``RR_CONTEXTS``.
    :return: DataFrame with ``context``, ``window``, ``feature``, ``coef``, ``se``,
        and ``pval`` in screening order.
    """
    rows = []
    for context in contexts:
        df_1_ = df[(df["y"] == 1) & (df["context"] == context)]
        df_0_ = df[(df["y"] == 0) & (df["context"] == context)]
        for scale in GENOMIC_FEATURE_SCALES:
            for feature in RR_SCREEN_FEATURE_ORDER:
                col = f"{feature}_{scale}"
                df_1 = df_1_[["locus", col]].dropna().assign(group=1)
                df_0 = df_0_[["locus", col]].dropna().assign(group=0)
                df_01 = pd.concat([df_1, df_0])
                df_x = df_01[[col]].apply(scipy.stats.zscore)
                try:
                    fit = sm.Logit(
                        df_01[["group"]], sm.add_constant(df_x, has_constant="add")
                    ).fit_regularized(disp=False)
                    coef = fit.params[col]
                    pval = fit.pvalues[col]
                    hci = fit.conf_int().loc[col].to_list()[1]
                    se = (hci - coef) / 1.96
                except Exception:  # pylint: disable=broad-except
                    coef, se, pval = [np.nan] * 3
                rows.append(
                    {
                        "context": context,
                        "window": scale,
                        "feature": feature,
                        "coef": coef,
                        "se": se,
                        "pval": pval,
                    }
                )

    return pd.DataFrame(rows)


def select_rr_features(screen_df: pd.DataFrame, cutoff: str) -> pd.DataFrame:
    """
    Select rr model features from the univariate feature screen.

    With the ``RR_BONFERRONI_CUTOFF`` token, a feature is selected when its p-value
    is at most ``RR_FEATURE_SELECTION_ALPHA / 4 / 13`` (``/ 8`` for CpG contexts).
    Otherwise ``cutoff`` is a raw p-value cutoff. Only the first passing scale (in
    screening order, i.e. the smallest) is kept for each context and feature.

    :param screen_df: Feature screen DataFrame (:func:`run_rr_feature_screen`).
    :param cutoff: ``RR_BONFERRONI_CUTOFF`` or a p-value cutoff (e.g. "1e-3").
    :return: Selected rows of ``screen_df``.
    """
    df = screen_df.copy()
    if cutoff == RR_BONFERRONI_CUTOFF:
        n_scales = len(GENOMIC_FEATURE_SCALES)
        n_cpg = len(GENOMIC_FEATURES) - len(METHYLATION_CORRELATED_FEATURES)
        df["threshold"] = np.where(
            df["context"].isin(CPG_CONTEXTS),
            RR_FEATURE_SELECTION_ALPHA / n_scales / n_cpg,
            RR_FEATURE_SELECTION_ALPHA / n_scales / len(GENOMIC_FEATURES),
        )
    else:
        df["threshold"] = float(cutoff)
    df = df[df["pval"] <= df["threshold"]]

    return df.drop_duplicates(subset=["context", "feature"]).drop(columns="threshold")


def get_rr_model_features(selected_df: pd.DataFrame, context: str) -> List[str]:
    """
    Get the ordered model features of a context from the selected features.

    Methylation-correlated features are dropped for CpG contexts.

    :param selected_df: Selected features (:func:`select_rr_features`).
    :param context: Trinucleotide context.
    :return: Feature column names (``{feature}_{scale}``) in selection order.
    """
    sel = selected_df[selected_df["context"] == context]
    if context in CPG_CONTEXTS:
        sel = sel[~sel["feature"].isin(METHYLATION_CORRELATED_FEATURES)]

    return list(sel["feature"] + "_" + sel["window"])


def fit_rr_context_model(
    df: pd.DataFrame, features: List[str]
) -> Optional[Dict[str, Any]]:
    """
    Fit the rr logistic model of one context.

    Port of the Gnocchi per-context model fit: sites with any missing feature are
    dropped, DNM rows are followed by background rows (training row order), the
    features are z-scored (ddof=0), rotated with a full-rank ``IncrementalPCA``, and
    a logistic regression is fit with ``fit_regularized(method="l1", alpha=0)``.

    The PCA is folded into the model coefficients so that the linear predictor of a
    window is ``intercept + sum(coefs * (x - mean) / sd)``. As in the legacy models,
    ``sd`` is the ddof=1 standard deviation (applied at prediction) while the model
    is trained on ddof=0 z-scores; this is kept to reproduce the legacy rr exactly.
    ``ave`` is the predicted probability of the average window (all principal
    components equal to 0).

    :param df: Training rows of one context (:func:`rr_training_ht_to_pandas`).
    :param features: Model features (:func:`get_rr_model_features`).
    :return: Dictionary of model parameters, or None if there is no training data.
    """
    df_1 = df[df["y"] == 1][["locus", *features]].dropna()
    df_0 = df[df["y"] == 0][["locus", *features]].dropna()
    if df_1.empty or df_0.empty:
        return None
    df_01 = pd.concat([df_1.assign(group=1), df_0.assign(group=0)])
    x = df_01[features]
    y = df_01["group"]

    x_scaled = x.apply(scipy.stats.zscore)
    pca = IncrementalPCA(n_components=len(features))
    x_pca = pca.fit_transform(x_scaled)
    fit = sm.Logit(y, sm.add_constant(x_pca, has_constant="add")).fit_regularized(
        method="l1", alpha=0.0, disp=False
    )

    params = np.asarray(fit.params, dtype=float)
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
        "converged": bool(fit.mle_retvals.get("converged", False)),
    }


def predict_rr_model(model: Dict[str, Any], df: pd.DataFrame) -> np.ndarray:
    """
    Predict the DNM probability of sites with a fitted rr model.

    :param model: Model parameters (:func:`fit_rr_context_model`).
    :param df: DataFrame with the model feature columns.
    :return: Predicted probabilities.
    """
    z = (df[model["features"]].to_numpy() - np.asarray(model["mean"])) / np.asarray(
        model["sd"]
    )
    eta = model["intercept"] + z @ np.asarray(model["coefs"])

    return 1.0 / (1.0 + np.exp(-eta))


def fit_rr_models(
    df: pd.DataFrame,
    screen_df: pd.DataFrame,
    cutoffs: List[str],
    cv_chrom_split: bool = False,
    contexts: List[str] = RR_CONTEXTS,
) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
    """
    Fit the rr models of all contexts for each feature selection cutoff.

    Models trained on all training sites have fold ``MODEL_FOLD_ALL``. With
    ``cv_chrom_split``, models are also fit with one chromosome-split fold held out
    (see ``RR_FOLDS``); these models have ``fold`` set to the held-out fold (the
    fold whose windows they are applied to) and their held-out log loss and AUC are
    returned.

    :param df: Training DataFrame (:func:`rr_training_ht_to_pandas`).
    :param screen_df: Feature screen DataFrame (:func:`run_rr_feature_screen`).
    :param cutoffs: Feature selection cutoff tokens.
    :param cv_chrom_split: Whether to also fit chromosome-split models. Default is
        False.
    :param contexts: Contexts to model. Default is ``RR_CONTEXTS``.
    :return: Tuple of the list of model dictionaries and a DataFrame of held-out
        metrics (empty if ``cv_chrom_split`` is False).
    """
    df = df.assign(_fold=df["locus"].str.split(":").str[0].map(_chrom_fold))
    folds = [MODEL_FOLD_ALL, *(RR_FOLDS if cv_chrom_split else [])]

    models = []
    metrics = []
    for cutoff in cutoffs:
        selected_df = select_rr_features(screen_df, cutoff)
        n_modeled = 0
        for context in contexts:
            features = get_rr_model_features(selected_df, context)
            if not features:
                continue
            ctx_df = df[df["context"] == context]
            for fold in folds:
                train_df = (
                    ctx_df
                    if fold == MODEL_FOLD_ALL
                    else ctx_df[ctx_df["_fold"] != fold]
                )
                model = fit_rr_context_model(train_df, features)
                if model is None:
                    continue
                models.append(
                    {"cutoff": cutoff, "fold": fold, "context": context, **model}
                )
                if fold == MODEL_FOLD_ALL:
                    n_modeled += 1
                    continue

                test_df = ctx_df[ctx_df["_fold"] == fold][["y", *features]].dropna()
                if test_df["y"].nunique() < 2:
                    continue
                pred = predict_rr_model(model, test_df)
                metrics.append(
                    {
                        "cutoff": cutoff,
                        "fold": fold,
                        "context": context,
                        "n_features": len(features),
                        "n_test": len(test_df),
                        "log_loss": log_loss(test_df["y"], pred, labels=[0, 1]),
                        "auc": roc_auc_score(test_df["y"], pred),
                    }
                )
        logger.info(
            "Feature selection cutoff %s: %d of %d contexts modeled.",
            cutoff,
            n_modeled,
            len(contexts),
        )

    return models, pd.DataFrame(metrics)


def rr_models_to_ht(models: List[Dict[str, Any]]) -> hl.Table:
    """
    Convert rr model dictionaries to a Hail Table keyed by cutoff, fold, and context.

    :param models: Model dictionaries (:func:`fit_rr_models`).
    :return: rr model Table.
    """
    schema = hl.tstruct(
        cutoff=hl.tstr,
        fold=hl.tstr,
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

    return hl.Table.parallelize(models, schema=schema).key_by(
        "cutoff", "fold", "context"
    )


# ---------------------------------------------------------------------------
# rr per window and context
# ---------------------------------------------------------------------------


def _rr_model_expr(
    ht: hl.Table,
    model: hl.Struct,
    masked_features: Optional[List[str]] = None,
    mask_autosomes: bool = False,
) -> hl.expr.Float64Expression:
    """
    Get the rr expression of one context model for each window of a feature Table.

    rr = sigmoid(intercept + sum(coefs * (x - mean) / sd)) / ave, and is missing when
    any (unmasked) model feature is missing for the window. Masked features are
    pinned to z = 0 (the genome-average window) on chrX/chrY windows (on all windows
    if ``mask_autosomes`` is True).

    :param ht: Feature Table keyed by ``interval``.
    :param model: Model row of the rr model Table.
    :param masked_features: Optional base feature names to mask. Default is None.
    :param mask_autosomes: Whether to also mask features on autosomal windows (used
        for the masking ablation). Default is False.
    :return: rr expression.
    """
    mask_expr = hl.bool(True) if mask_autosomes else ~ht.interval.start.in_autosome()
    eta = hl.float64(model.intercept)
    for feature, mean, sd, coef in zip(
        model.features, model.mean, model.sd, model.coefs
    ):
        z = (ht[feature] - mean) / sd
        if masked_features and feature.rsplit("_", 1)[0] in masked_features:
            z = hl.if_else(mask_expr, 0.0, z)
        eta = eta + coef * z

    return (1.0 / (1.0 + hl.exp(-eta))) / model.ave


def compute_rr_ht(
    features_ht: hl.Table,
    models: List[hl.Struct],
    masked_features: Optional[List[str]] = None,
    out_of_fold: bool = False,
) -> hl.Table:
    """
    Compute rr per 1kb window and trinucleotide context.

    Contexts without a model, and windows missing a model feature for a context, get
    an rr of 1.0. As in the legacy rr Tables, windows are kept only if at least one
    modeled context has an rr for the window.

    :param features_ht: Feature Table (:func:`import_genomic_features`).
    :param models: Rows of the rr model Table for one feature selection cutoff.
    :param masked_features: Optional base feature names to mask on chrX/chrY
        windows (see :func:`_rr_model_expr`). Default is None.
    :param out_of_fold: Whether to use the chromosome-split models, applying to each
        window the model trained without the window's fold. If False, the models
        trained on all sites are used. Default is False.
    :return: rr Table keyed by ``interval`` with an ``adj_r`` dict.
    """
    by_context = {}
    for model in models:
        by_context.setdefault(model.context, {})[model.fold] = model

    fold_expr = _chrom_fold_expr(features_ht.interval.start.contig)
    rr_exprs = {}
    for context, ctx_models in by_context.items():
        if out_of_fold:
            if not all(f in ctx_models for f in RR_FOLDS):
                continue
            rr_exprs[context] = hl.if_else(
                fold_expr == RR_FOLDS[0],
                _rr_model_expr(features_ht, ctx_models[RR_FOLDS[0]], masked_features),
                _rr_model_expr(features_ht, ctx_models[RR_FOLDS[1]], masked_features),
            )
        elif MODEL_FOLD_ALL in ctx_models:
            rr_exprs[context] = _rr_model_expr(
                features_ht, ctx_models[MODEL_FOLD_ALL], masked_features
            )

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


def clip_rr_ht(rr_ht: hl.Table, lower: float, upper: float) -> hl.Table:
    """
    Clip rr values to ``[lower, upper]``.

    :param rr_ht: rr Table (:func:`compute_rr_ht`).
    :param lower: Lower bound.
    :param upper: Upper bound.
    :return: rr Table with clipped ``adj_r`` values.
    """
    return rr_ht.annotate(
        adj_r=rr_ht.adj_r.map_values(lambda v: hl.max(hl.min(v, upper), lower))
    )


def _adj_freq_idx(ht: hl.Table) -> int:
    """
    Get the index of the adj frequency element in a per-variant expected Table.

    :param ht: Per-variant expected Table with an ``exomes_freq_meta`` global.
    :return: Index of the adj frequency element.
    """
    return hl.eval(ht.exomes_freq_meta).index(ADJ_FREQ_META)


def _filter_to_syn_mane_select(per_variant_ht: hl.Table) -> hl.Table:
    """
    Filter the per-variant expected Table to synonymous MANE Select rows.

    :param per_variant_ht: Per-variant expected Table from the constraint pipeline.
    :return: Filtered Table.
    """
    ht = per_variant_ht.filter(
        (per_variant_ht.annotation == SYNONYMOUS_CSQ) & per_variant_ht.mane_select
    )

    return ht.naive_coalesce(200)


def compute_rr_normalization_factors(
    per_variant_ht: hl.Table, rr_ht: hl.Table
) -> Dict[str, float]:
    """
    Compute the per-context rr normalization factors on exome synonymous sites.

    The factor of a context is the expected-weighted mean rr over synonymous MANE
    Select sites of the per-variant expected Table (unadjusted adj
    ``expected_variants``) that fall in an rr window:
    ``sum(expected * rr) / sum(expected)``.

    :param per_variant_ht: Per-variant expected Table from the constraint pipeline.
    :param rr_ht: rr Table (:func:`compute_rr_ht`, optionally clipped).
    :return: Dictionary of context to normalization factor.
    """
    adj_idx = _adj_freq_idx(per_variant_ht)
    ht = _filter_to_syn_mane_select(per_variant_ht)
    rr_expr = rr_ht[ht.locus].adj_r.get(ht.context)
    exp_expr = ht.expected_variants[adj_idx]
    sums = ht.aggregate(
        hl.agg.filter(
            hl.is_defined(rr_expr) & hl.is_defined(exp_expr),
            hl.agg.group_by(
                ht.context,
                hl.struct(num=hl.agg.sum(exp_expr * rr_expr), den=hl.agg.sum(exp_expr)),
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


def run_mask_ablation(
    features_ht: hl.Table, models: List[hl.Struct], masked_features: List[str]
) -> pd.DataFrame:
    """
    Measure on autosomal windows what masking features costs.

    For each modeled context, compares rr computed with all features to rr with
    ``masked_features`` pinned to the genome-average window, on autosomal windows
    where both are defined.

    :param features_ht: Feature Table (:func:`import_genomic_features`).
    :param models: Rows of the rr model Table (fold ``MODEL_FOLD_ALL``) for one
        feature selection cutoff.
    :param masked_features: Base feature names to mask.
    :return: DataFrame with ``context``, ``n_windows``, ``corr``, ``sd_retained``,
        and ``median_rel_diff`` (median of masked / full - 1).
    """
    ht = features_ht.filter(features_ht.interval.start.in_autosome())
    models = [m for m in models if m.fold == MODEL_FOLD_ALL]
    ht = ht.select(
        **{
            m.context: hl.struct(
                full=_rr_model_expr(ht, m),
                masked=_rr_model_expr(ht, m, masked_features, mask_autosomes=True),
            )
            for m in models
        }
    )
    res = ht.aggregate(
        hl.struct(
            **{
                m.context: hl.agg.filter(
                    hl.is_defined(ht[m.context].full)
                    & hl.is_defined(ht[m.context].masked),
                    hl.struct(
                        n=hl.agg.count(),
                        corr=hl.agg.corr(ht[m.context].full, ht[m.context].masked),
                        full_stats=hl.agg.stats(ht[m.context].full),
                        masked_stats=hl.agg.stats(ht[m.context].masked),
                        median_rel_diff=hl.agg.approx_median(
                            ht[m.context].masked / ht[m.context].full - 1
                        ),
                    ),
                )
                for m in models
            }
        )
    )

    return pd.DataFrame(
        [
            {
                "context": c,
                "n_windows": r.n,
                "corr": r.corr,
                "sd_retained": (
                    r.masked_stats.stdev / r.full_stats.stdev
                    if r.full_stats.stdev
                    else np.nan
                ),
                "median_rel_diff": r.median_rel_diff,
            }
            for c, r in res.items()
        ]
    )


# ---------------------------------------------------------------------------
# Feature selection cutoff evaluation
# ---------------------------------------------------------------------------


def evaluate_rr_cutoffs(
    per_variant_ht: hl.Table,
    rr_hts: Dict[str, hl.Table],
    cv_metrics_df: Optional[pd.DataFrame] = None,
    n_contexts_modeled: Optional[Dict[str, int]] = None,
    loeuf_ht: Optional[hl.Table] = None,
    syn_z_outlier_threshold: float = 8.0,
    min_expected: float = 0.0,
) -> pd.DataFrame:
    """
    Evaluate synonymous calibration of the rr Tables of each feature selection cutoff.

    Synonymous MANE Select rows of the per-variant expected Table are aggregated per
    transcript (adj observed, unadjusted adj expected, and expected multiplied by
    the rr of each cutoff). For each cutoff (and the unadjusted expected), reports:

        - ``syn_mean_z2``: mean squared raw z-score (dispersion).
        - ``n_syn_outliers``: transcripts with |raw z| > ``syn_z_outlier_threshold``.
        - ``oe_decile_{i}``: synonymous o/e (sum obs / sum exp) in each LOEUF decile,
          when ``loeuf_ht`` is provided.
        - ``n_contexts_modeled``, ``heldout_mean_log_loss``, ``heldout_mean_auc``
          when provided.

    Out-of-fold rr Tables should be used to avoid evaluating genes on the same DNMs
    used to train the rr models.

    :param per_variant_ht: Per-variant expected Table from the constraint pipeline.
    :param rr_hts: Dictionary of cutoff token to rr Table.
    :param cv_metrics_df: Optional held-out metrics (:func:`fit_rr_models`). Default
        is None.
    :param n_contexts_modeled: Optional dictionary of cutoff token to number of
        modeled contexts. Default is None.
    :param loeuf_ht: Optional release-format constraint Table with ``transcript``
        and ``lof.oe_ci.upper_bin_decile``. Default is None.
    :param syn_z_outlier_threshold: Absolute raw z-score above which a transcript is
        a synonymous outlier. Default is 8.0.
    :param min_expected: Minimum unadjusted expected count for a transcript to be
        included. Default is 0.0 (transcripts with expected > 0).
    :return: DataFrame with one row per cutoff (and "unadjusted").
    """
    adj_idx = _adj_freq_idx(per_variant_ht)
    ht = _filter_to_syn_mane_select(per_variant_ht)
    obs_expr = hl.or_else(ht.calibrate_mu.observed_variants[adj_idx], 0)
    exp_expr = ht.expected_variants[adj_idx]
    labels = ["unadjusted", *rr_hts]
    exp_fields = {"exp_0": exp_expr}
    for i, rr_ht in enumerate(rr_hts.values(), start=1):
        rr_expr = hl.or_else(rr_ht[ht.locus].adj_r.get(ht.context), 1.0)
        exp_fields[f"exp_{i}"] = exp_expr * rr_expr

    ht = ht.select("transcript", obs=obs_expr, **exp_fields)
    ht = ht.group_by("transcript").aggregate(
        obs=hl.agg.sum(ht.obs), **{f: hl.agg.sum(ht[f]) for f in exp_fields}
    )
    ht = ht.filter((ht.exp_0 > 0) & (ht.exp_0 >= min_expected))
    if loeuf_ht is not None:
        decile_ht = loeuf_ht.key_by("transcript").select(
            decile=loeuf_ht.lof.oe_ci.upper_bin_decile
        )
        decile_ht = decile_ht.distinct()
        ht = ht.annotate(decile=decile_ht[ht.transcript].decile)
    ht = ht.checkpoint(new_temp_file("rr_cutoff_evaluation", "ht"))

    rows = []
    for i, label in enumerate(labels):
        exp_f = ht[f"exp_{i}"]
        z_expr = calculate_raw_z_score(ht.obs, exp_f)
        agg = {
            "syn_mean_z2": hl.agg.mean(z_expr**2),
            "n_syn_outliers": hl.agg.count_where(
                hl.abs(z_expr) > syn_z_outlier_threshold
            ),
            "n_transcripts": hl.agg.count(),
        }
        if loeuf_ht is not None:
            agg["oe_by_decile"] = hl.agg.group_by(
                ht.decile, hl.agg.sum(ht.obs) / hl.agg.sum(exp_f)
            )
        res = ht.aggregate(hl.struct(**agg))
        row = {
            "cutoff": label,
            "n_transcripts": res.n_transcripts,
            "syn_mean_z2": res.syn_mean_z2,
            "n_syn_outliers": res.n_syn_outliers,
        }
        if loeuf_ht is not None:
            row.update(
                {
                    f"oe_decile_{d}": v
                    for d, v in sorted(
                        res.oe_by_decile.items(), key=lambda x: (x[0] is None, x[0])
                    )
                    if d is not None
                }
            )
        if n_contexts_modeled is not None:
            row["n_contexts_modeled"] = n_contexts_modeled.get(label)
        if cv_metrics_df is not None and not cv_metrics_df.empty:
            m = cv_metrics_df[cv_metrics_df["cutoff"] == label]
            row["heldout_mean_log_loss"] = m["log_loss"].mean() if len(m) else np.nan
            row["heldout_mean_auc"] = m["auc"].mean() if len(m) else np.nan
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Regional genomic features for chrX/chrY windows
# ---------------------------------------------------------------------------


def _read_tsv(path: str, **kwargs) -> pd.DataFrame:
    """
    Read a (optionally compressed) tab-delimited file from a local or GCS path.

    :param path: Local or ``gs://`` path.
    :param kwargs: Keyword arguments passed to ``pandas.read_csv``.
    :return: DataFrame, or an iterator of DataFrames if ``chunksize`` is given.
    """
    if path.startswith("gs://"):
        # hl.hadoop_open decompresses .gz/.bgz files based on the extension.
        return pd.read_csv(hl.hadoop_open(path, "r"), sep="\t", **kwargs)

    return pd.read_csv(path, sep="\t", compression="infer", **kwargs)


def load_bed_intervals(path: str, contig: Optional[str] = None) -> pd.DataFrame:
    """
    Read BED intervals as ``chrom``, ``start``, ``end``, optionally for one contig.

    :param path: BED path.
    :param contig: Optional contig to keep. Default is None.
    :return: DataFrame of intervals sorted by chrom and start.
    """
    df = _read_tsv(path, header=None, usecols=[0, 1, 2], low_memory=False)
    df.columns = ["chrom", "start", "end"]
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = np.where(
        df["chrom"].str.startswith("chr"), df["chrom"], "chr" + df["chrom"]
    )
    if contig is not None:
        df = df[df["chrom"] == contig]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).astype({"start": np.int64, "end": np.int64})
    df = df[df["end"] > df["start"]]

    return df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping (or touching) intervals per chromosome.

    :param df: DataFrame with ``chrom``, ``start``, and ``end``.
    :return: Merged, sorted, non-overlapping intervals.
    """
    out = []
    for chrom, sub in df.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        starts = sub["start"].to_numpy(dtype=np.int64)
        ends = np.maximum.accumulate(sub["end"].to_numpy(dtype=np.int64))
        # A new merged interval starts where the start is past every previous end.
        new = np.ones(len(starts), dtype=bool)
        new[1:] = starts[1:] > ends[:-1]
        out.append(
            pd.DataFrame(
                {
                    "chrom": chrom,
                    "start": starts[new],
                    "end": np.maximum.reduceat(ends, np.flatnonzero(new)),
                }
            )
        )

    if not out:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    return pd.concat(out, ignore_index=True)


def _covered_bases(
    positions: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> np.ndarray:
    """
    Count bases in ``[0, position)`` covered by merged, sorted intervals.

    :param positions: Query positions (0-based, exclusive).
    :param starts: Merged interval starts (sorted).
    :param ends: Merged interval ends.
    :return: Number of covered bases before each position.
    """
    cum = np.concatenate([[0], np.cumsum(ends - starts)])
    k = np.searchsorted(starts, positions, side="left")
    prev = np.maximum(k - 1, 0)
    partial = np.clip(np.minimum(positions, ends[prev]) - starts[prev], 0, None)

    return np.where(k > 0, cum[prev] + partial, 0)


def compute_coverage_feature(
    windows: pd.DataFrame, intervals: pd.DataFrame
) -> pd.Series:
    """
    Get the percent of each window covered by merged track intervals.

    Equivalent to ``bedtools coverage`` (bases covered by at least one interval).
    Returns NaN for every window when the track has no intervals on the windows'
    contig (a missing track, rather than a real 0% coverage).

    :param windows: Windows of one contig with ``start`` and ``end``.
    :param intervals: Merged track intervals of the same contig.
    :return: Percent coverage per window.
    """
    if intervals.empty:
        return pd.Series(np.nan, index=windows.index)

    starts = intervals["start"].to_numpy(dtype=np.int64)
    ends = intervals["end"].to_numpy(dtype=np.int64)
    w_start = windows["start"].to_numpy(dtype=np.int64)
    w_end = windows["end"].to_numpy(dtype=np.int64)
    covered = _covered_bases(w_end, starts, ends) - _covered_bases(
        w_start, starts, ends
    )

    return pd.Series(
        covered / np.clip(w_end - w_start, 1, None) * 100.0, index=windows.index
    )


def compute_distance_feature(
    mids: pd.DataFrame, intervals: pd.DataFrame, chrom_size: int
) -> pd.Series:
    """
    Get the distance from each window midpoint to the nearest interval.

    The distance is 0 inside an interval and is expressed as a percent of the
    chromosome length (``bedtools closest -d`` followed by the Gnocchi
    normalization). NaN when the track has no intervals on the contig.

    :param mids: Window midpoints of one contig with ``start``.
    :param intervals: Track intervals of the same contig.
    :param chrom_size: Chromosome length.
    :return: Normalized distance per window.
    """
    if intervals.empty:
        return pd.Series(np.nan, index=mids.index)

    starts = intervals["start"].to_numpy(dtype=np.int64)
    ends = intervals["end"].to_numpy(dtype=np.int64)
    pos = mids["start"].to_numpy(dtype=np.int64)
    dist = np.maximum(
        np.maximum(starts[None, :] - pos[:, None], pos[:, None] - ends[None, :]), 0
    ).min(axis=1)

    return pd.Series(dist / chrom_size * 100.0, index=mids.index)


def fetch_reference_sequence(
    contig: str, start: int, end: int, cache_dir: str, chunk: int = 1_000_000
) -> str:
    """
    Fetch hg38 reference sequence from the UCSC REST API, with a local cache.

    :param contig: Contig name.
    :param start: 0-based start.
    :param end: 0-based exclusive end.
    :param cache_dir: Local directory used to cache fetched sequence chunks.
    :param chunk: Chunk size of each request. Default is 1,000,000.
    :return: Upper-case sequence.
    """
    os.makedirs(cache_dir, exist_ok=True)
    parts = []
    for s in range(start, end, chunk):
        e = min(s + chunk, end)
        cache = os.path.join(cache_dir, f"{contig}_{s}_{e}.txt")
        if not os.path.exists(cache):
            url = f"{UCSC_SEQUENCE_API}?genome=hg38;chrom={contig};start={s};end={e}"
            with urllib.request.urlopen(url, timeout=300) as resp:
                seq = json.load(resp)["dna"].upper()
            with open(cache, "w") as fh:
                fh.write(seq)
        with open(cache) as fh:
            parts.append(fh.read())

    return "".join(parts)


def compute_gc_feature(
    windows: pd.DataFrame, gc_cumsum: np.ndarray, offset: int
) -> pd.Series:
    """
    Get the percent G+C of each window.

    N bases count toward the denominator (``hgGcPercent -doGaps``), matching the
    published feature values.

    :param windows: Windows of one contig with ``start`` and ``end``.
    :param gc_cumsum: Cumulative G+C count of the contig sequence starting at
        ``offset`` (with a leading 0).
    :param offset: Start position of the fetched sequence.
    :return: Percent G+C per window.
    """
    s = windows["start"].to_numpy(dtype=np.int64) - offset
    e = windows["end"].to_numpy(dtype=np.int64) - offset

    return pd.Series(
        (gc_cumsum[e] - gc_cumsum[s]) / (e - s) * 100.0, index=windows.index
    )


def load_feature_windows(
    tracks_dir: str,
    contig: str,
    element_ids: Optional[set] = None,
    mids: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Load the window (or midpoint) BEDs of each scale for one contig.

    :param tracks_dir: Genomic feature source track directory.
    :param contig: Contig to load.
    :param element_ids: Optional 1kb window element IDs to keep. Default is None.
    :param mids: Whether to load the midpoint BEDs. Default is False.
    :return: Dictionary of scale to DataFrame (``chrom``, ``start``, ``end``,
        ``element_id``) indexed by element ID.
    """
    beds = GENOMIC_FEATURE_MID_BEDS if mids else GENOMIC_FEATURE_WINDOW_BEDS
    out = {}
    for scale, name in beds.items():
        keep = []
        for chunk in _read_tsv(
            f"{tracks_dir}/{name}",
            header=None,
            names=["chrom", "start", "end", "element_id"],
            chunksize=1_000_000,
        ):
            sub = chunk[chunk["chrom"] == contig]
            if element_ids is not None:
                sub = sub[sub["element_id"].isin(element_ids)]
            keep.append(sub)
        out[scale] = pd.concat(keep).set_index("element_id", drop=False)

    return out


def build_feature_matrix(
    tracks_dir: str,
    contig: str,
    seq_cache_dir: str,
    masked_features: List[str],
    element_ids: Optional[set] = None,
) -> pd.DataFrame:
    """
    Build the 1kb regional genomic feature matrix of one contig.

    Computes the coverage features (``GENOMIC_FEATURE_COVERAGE_TRACKS``), distance
    features (``GENOMIC_FEATURE_DISTANCE_TRACKS``), and ``GC_content`` at each scale.
    ``masked_features`` (features without source data for the contig) are written
    as NaN, as are features whose track has no intervals on the contig.

    :param tracks_dir: Genomic feature source track directory.
    :param contig: Contig to build.
    :param seq_cache_dir: Local directory used to cache reference sequence.
    :param masked_features: Base feature names to write as NaN.
    :param element_ids: Optional 1kb window element IDs to build. Default is None
        (all windows of the contig).
    :return: DataFrame with ``element_id`` and ``GENOMIC_FEATURE_COLUMNS``.
    """
    chrom_sizes = _read_tsv(f"{tracks_dir}/hg38.chrom.sizes", header=None)
    chrom_size = int(dict(zip(chrom_sizes[0], chrom_sizes[1]))[contig])

    windows = load_feature_windows(tracks_dir, contig, element_ids)
    mids = load_feature_windows(tracks_dir, contig, element_ids, mids=True)
    ids = windows["1k"]["element_id"]
    frame = pd.DataFrame({"element_id": ids.to_numpy()}, index=ids.to_numpy())
    logger.info("Building features for %d %s windows...", len(frame), contig)

    # Fetch the reference sequence spanned by all window scales once.
    seq_start = int(min(w["start"].min() for w in windows.values()))
    seq_end = int(max(w["end"].max() for w in windows.values()))
    seq = fetch_reference_sequence(contig, seq_start, seq_end, seq_cache_dir)
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    gc_cumsum = np.concatenate(
        [[0], np.cumsum((arr == ord("G")) | (arr == ord("C")), dtype=np.int32)]
    )
    del arr, seq

    coverage_tracks = {
        f: merge_intervals(load_bed_intervals(f"{tracks_dir}/{t}", contig))
        for f, t in GENOMIC_FEATURE_COVERAGE_TRACKS.items()
    }
    distance_tracks = {
        f: load_bed_intervals(f"{tracks_dir}/{t}", contig)
        for f, t in GENOMIC_FEATURE_DISTANCE_TRACKS.items()
    }
    for scale in GENOMIC_FEATURE_SCALES:
        w = windows[scale]
        for feature, intervals in coverage_tracks.items():
            frame[f"{feature}_{scale}"] = compute_coverage_feature(w, intervals)
        for feature, intervals in distance_tracks.items():
            frame[f"{feature}_{scale}"] = compute_distance_feature(
                mids[scale], intervals, chrom_size
            )
        frame[f"GC_content_{scale}"] = compute_gc_feature(w, gc_cumsum, seq_start)

    for col in GENOMIC_FEATURE_COLUMNS:
        if col not in frame.columns or col.rsplit("_", 1)[0] in masked_features:
            frame[col] = np.nan

    return frame[["element_id", *GENOMIC_FEATURE_COLUMNS]].reset_index(drop=True)


def validate_feature_builder(
    tracks_dir: str,
    published_features_path: str,
    seq_cache_dir: str,
    contig: str = "chr21",
    n_windows: int = 3000,
) -> pd.DataFrame:
    """
    Rebuild autosomal windows and compare them to the published feature matrix.

    Only the features computed by :func:`build_feature_matrix` are compared.

    :param tracks_dir: Genomic feature source track directory.
    :param published_features_path: Path of the published feature matrix.
    :param seq_cache_dir: Local directory used to cache reference sequence.
    :param contig: Autosome to validate on. Default is "chr21".
    :param n_windows: Number of published windows of ``contig`` to rebuild. Default
        is 3000.
    :return: DataFrame with ``column``, ``n``, ``r``, ``median_abs_diff``,
        ``max_abs_diff``, and ``passed`` (r > 0.999 and median |diff| < 0.05).
    """
    published = []
    for chunk in _read_tsv(published_features_path, chunksize=200_000):
        sub = chunk[chunk["element_id"].str.startswith(f"{contig}-")]
        published.append(sub)
        if sum(len(p) for p in published) >= n_windows:
            break
    published = pd.concat(published).head(n_windows)

    built = build_feature_matrix(
        tracks_dir, contig, seq_cache_dir, [], set(published["element_id"])
    )
    merged = published.merge(built, on="element_id", suffixes=("_pub", "_built"))
    computed = [
        *GENOMIC_FEATURE_COVERAGE_TRACKS,
        *GENOMIC_FEATURE_DISTANCE_TRACKS,
        "GC_content",
    ]
    rows = []
    for col in [f"{f}_{s}" for s in GENOMIC_FEATURE_SCALES for f in computed]:
        a = merged[f"{col}_built"].to_numpy(dtype=float)
        b = merged[f"{col}_pub"].to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        d = np.abs(a[ok] - b[ok])
        r = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 1 else np.nan
        rows.append(
            {
                "column": col,
                "n": int(ok.sum()),
                "r": r,
                "median_abs_diff": float(np.median(d)) if len(d) else np.nan,
                "max_abs_diff": float(d.max()) if len(d) else np.nan,
                "passed": bool(r > 0.999 and np.median(d) < 0.05) if len(d) else False,
            }
        )

    return pd.DataFrame(rows)
