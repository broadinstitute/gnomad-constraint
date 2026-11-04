"""Compare rr training backgrounds (what the DNMs are contrasted with), locally.

Inputs:

- Ruchit Panchal's rr evaluation export (``regional_rate/rr_eval_export/`` under
  gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/): per (transcript, 1kb window, context)
  synonymous possible (``n_sites``), observed and expected counts, the window
  features, and the transcript table. The cells file must be decompressed.
- The Satterstrom 2026 ASD and DD DNM files (``get_dnm_cohort_resources``):
  synonymous SNVs (``isSyn``), deduplicated by variant.
- UCSC hg38.2bit, for the trinucleotide context of each DNM.
- The v4.1.2 release TSV, for QC gene flags in the evaluator.

Backgrounds are built from the exported cells, so they cover the synonymous sites of
Ensembl MANE Select transcripts, and DNMs are restricted to the same windows and
contexts. Every model feature is a window feature, so a background given as a weight
per (window, context) is the same as sampling its sites. Weights are 5% of the site
counts, as in the pipeline. The backgrounds:

- ``absent``: possible synonymous variants not seen in gnomAD v4 exomes (the current
  rule, but per variant; the pipeline keeps one row per locus).
- ``absent_context_matched``: the same, rescaled so each context has the same total
  as ``all_sites`` (sampling the absent sites by context, in proportion to all sites).
- ``absent_equal_per_context``: the same, rescaled so every context has the same total
  (even sampling by context).
- ``absent_reweighted``: the same, each weighted by 1 / (1 - q_exp), where q_exp =
  exp / n_sites is the mutation model's expected chance of being seen in gnomAD.
- ``all_sites``: every possible synonymous variant.
- ``polymorphisms``: variants seen in gnomAD, scaled like ``absent_context_matched``
  (a stand-in for untransmitted parental variants).

Each background gets ridge-penalized logistic fits per context on all features (as
``RR_RIDGE_CUTOFF``), with one per-DNM penalty (see ``PENALTIES``) chosen on held-out
DNMs (odd/even chromosome split). Each rr is also scored at a range of strengths
(rr^alpha), with alpha chosen on one chromosome half and scored on the other. rr is applied multiplicatively (as in the pipeline) or
saturation-aware: expected = n_sites * (1 - (1 - q)^rr), with q = exp / n_sites, then
normalized per context like the multiplicative rr.
"""

import argparse
import logging
import struct
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rr_shrinkage_experiments import (
    CONTEXTS,
    FOLDS,
    Evaluator,
    allowed_columns,
    chrom_fold,
)
from sklearn.linear_model import LogisticRegression

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

BACKGROUND_FRACTION = 0.05
PENALTIES = [3.0, 10.0, 30.0, 100.0]
"""Per-DNM ridge penalties: summed log loss + penalty * n_dnm / 2 * sum(coef^2). The
penalty is scaled by the number of DNMs, not of rows, so backgrounds of different
sizes get the same shrinkage (in a rare-outcome fit, the information about each
coefficient grows with the number of cases, not controls). 10 per DNM is about the
pipeline's 0.3 per site with its 5% absent background."""
STRENGTHS = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
CLASSES = ["A/T", "C non-CpG", "CpG"]
BACKGROUNDS = [
    "absent",
    "absent_context_matched",
    "absent_equal_per_context",
    "absent_reweighted",
    "all_sites",
    "polymorphisms",
]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


class TwoBit:
    """Minimal random-access reader for UCSC .2bit files."""

    def __init__(self, path: str):
        """
        Read the sequence index.

        :param path: Path to the .2bit file.
        """
        self.f = open(path, "rb")  # pylint: disable=consider-using-with
        sig, _, n_seq, _ = struct.unpack("<IIII", self.f.read(16))
        assert sig == 0x1A412743, "not a little-endian .2bit file"
        self.offsets = {}
        for _ in range(n_seq):
            name = self.f.read(ord(self.f.read(1))).decode()
            self.offsets[name] = struct.unpack("<I", self.f.read(4))[0]
        self.records = {}

    def _record(self, name: str) -> tuple:
        """Read and cache a sequence's size, N blocks and DNA offset."""
        if name not in self.records:
            self.f.seek(self.offsets[name])
            size, n_n = struct.unpack("<II", self.f.read(8))
            n_starts = np.frombuffer(self.f.read(4 * n_n), dtype="<u4")
            n_sizes = np.frombuffer(self.f.read(4 * n_n), dtype="<u4")
            (n_mask,) = struct.unpack("<I", self.f.read(4))
            self.f.seek(8 * n_mask + 4, 1)
            self.records[name] = (size, n_starts, n_sizes, self.f.tell())
        return self.records[name]

    def get(self, name: str, start: int, end: int) -> str:
        """
        Get the uppercase sequence of ``[start, end)`` (0-based).

        :param name: Sequence name, e.g. "chr1".
        :param start: 0-based start.
        :param end: 0-based exclusive end.
        :return: Sequence.
        """
        _, n_starts, n_sizes, dna = self._record(name)
        first, last = start // 4, (end - 1) // 4
        self.f.seek(dna + first)
        packed = np.frombuffer(self.f.read(last - first + 1), dtype=np.uint8)
        codes = np.stack([(packed >> s) & 3 for s in (6, 4, 2, 0)], axis=1).ravel()
        seq = list("TCAG"[c] for c in codes[start - 4 * first : end - 4 * first])
        hit = (n_starts < end) & (n_starts + n_sizes > start)
        for s, n in zip(n_starts[hit], n_sizes[hit]):
            for i in range(max(s, start), min(s + n, end)):
                seq[i - start] = "N"
        return "".join(seq)


def load_dnms(paths: List[str], twobit_path: str) -> pd.DataFrame:
    """
    Load synonymous SNV DNMs with their window and strand-collapsed context.

    :param paths: Satterstrom DNM files.
    :param twobit_path: hg38.2bit path.
    :return: One row per distinct DNM with ``win``, ``context`` and ``contig``.
    """
    df = pd.concat(
        [
            pd.read_csv(
                p, sep=",", quotechar='"', usecols=["Variant", "isSyn"], dtype=str
            )
            for p in paths
        ]
    )
    df = df[df.isSyn.str.upper() == "TRUE"].drop_duplicates("Variant")
    parts = df.Variant.str.split(":", expand=True)
    df = df.assign(
        contig=parts[0], pos=parts[1].astype(int), ref=parts[2], alt=parts[3]
    )
    df = df[(df.ref.str.len() == 1) & (df.alt.str.len() == 1)]
    tb = TwoBit(twobit_path)
    tri = [tb.get(c, p - 2, p + 1) for c, p in zip(df.contig, df.pos)]
    df["tri"] = tri
    ok = df.tri.str[1] == df.ref
    logger.info(
        "%d synonymous SNV DNMs; %d with a reference mismatch dropped.",
        len(df),
        (~ok).sum(),
    )
    df = df[ok]
    df["context"] = np.where(
        df.ref.isin(["G", "T"]), df.tri.str.translate(COMPLEMENT).str[::-1], df.tri
    )
    start = df.pos // 1000 * 1000
    df["win"] = df.contig + "-" + start.astype(str) + "-" + (start + 1000).astype(str)
    return df[["contig", "pos", "win", "context"]]


def cell_background(ev: Evaluator) -> pd.DataFrame:
    """
    Collapse the exported cells to one row per (window, context) with the backgrounds.

    Overlapping transcripts share sites, so counts take the maximum over transcripts.

    :param ev: Evaluator with the exported cells.
    :return: Rows with ``n``, ``obs``, ``exp`` and one weight column per background.
    """
    agg = ev.cells.groupby(["win", "context"], as_index=False)[
        ["n_sites", "obs", "exp"]
    ].max()
    agg = agg.rename(columns={"n_sites": "n"})
    agg = agg[agg.n > 0]
    f = BACKGROUND_FRACTION
    absent = (agg.n - agg.obs).clip(lower=0)
    agg["absent"] = f * absent
    ctx_all = agg.groupby("context").n.transform("sum")
    ctx_absent = absent.groupby(agg.context).transform("sum")
    agg["absent_context_matched"] = f * absent * ctx_all / ctx_absent
    agg["absent_equal_per_context"] = (
        f * absent * ctx_absent.groupby(agg.context).first().mean() / ctx_absent
    )
    q_exp = (agg.exp / agg.n).clip(upper=0.95)
    agg["absent_reweighted"] = f * absent / (1 - q_exp)
    agg["all_sites"] = f * agg.n
    agg["polymorphisms"] = (
        f * agg.obs * ctx_all / agg.obs.groupby(agg.context).transform("sum")
    )
    logger.info(
        "Fraction of possible synonymous variants seen in gnomAD: %s",
        (
            agg.groupby(class_of(agg.context)).obs.sum()
            / agg.groupby(class_of(agg.context)).n.sum()
        )
        .round(3)
        .to_dict(),
    )
    return agg


def class_of(context: pd.Series) -> np.ndarray:
    """
    Get the context class of each trinucleotide context.

    :param context: Strand-collapsed contexts (middle base A or C).
    :return: "A/T", "C non-CpG" or "CpG".
    """
    return np.where(
        context.str[1] == "A",
        "A/T",
        np.where(context.str[2] == "G", "CpG", "C non-CpG"),
    )


def fit_weighted(x: np.ndarray, y: np.ndarray, w: np.ndarray, penalty: float) -> dict:
    """
    Fit a ridge logistic model on weighted rows with z-scored features.

    :param x: Features.
    :param y: Outcome (1 = DNM).
    :param w: Row weights.
    :param penalty: Per-DNM penalty (see ``PENALTIES``).
    :return: Model with ``mean``, ``sd``, ``coef`` and ``intercept``.
    """
    mean = np.average(x, axis=0, weights=w)
    sd = np.sqrt(np.average((x - mean) ** 2, axis=0, weights=w))
    sd[sd == 0] = 1.0
    model = LogisticRegression(C=1.0 / (y.sum() * penalty), max_iter=5000)
    model.fit((x - mean) / sd, y, sample_weight=w)
    return {
        "mean": mean,
        "sd": sd,
        "coef": model.coef_[0],
        "intercept": model.intercept_[0],
    }


def log_loss(model: dict, x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Get the summed weighted log loss of a model."""
    eta = model["intercept"] + ((x - model["mean"]) / model["sd"]) @ model["coef"]
    p = np.clip(1 / (1 + np.exp(-eta)), 1e-12, 1 - 1e-12)
    return float(-np.sum(w * (y * np.log(p) + (1 - y) * np.log(1 - p))))


def fit_context(
    context: str,
    dnm: pd.DataFrame,
    bg: pd.DataFrame,
    weight: str,
    penalties: List[float],
) -> Tuple[str, Dict]:
    """
    Fit one context's models for each fold and penalty, with held-out log loss.

    :param context: Trinucleotide context.
    :param dnm: DNM rows of the context with features and ``fold``.
    :param bg: Background rows of the context with features, ``fold`` and the weight
        column.
    :param weight: Background weight column.
    :param penalties: Penalties to fit.
    :return: Context and {(fold, penalty): (model, held-out loss)}; fold "all" has no loss.
    """
    cols = allowed_columns(context)
    d = dnm.dropna(subset=cols)
    b = bg[bg[weight] > 0].dropna(subset=cols)
    x = np.vstack([d[cols].to_numpy(), b[cols].to_numpy()])
    y = np.concatenate([np.ones(len(d)), np.zeros(len(b))])
    w = np.concatenate([np.ones(len(d)), b[weight].to_numpy()])
    fold = np.concatenate([d.fold.to_numpy(), b.fold.to_numpy()])
    out = {}
    for pen in penalties:
        for f in (*FOLDS, "all"):
            train = fold != f
            model = fit_weighted(x[train], y[train], w[train], pen)
            model["cols"] = cols
            loss = (
                None if f == "all" else log_loss(model, x[~train], y[~train], w[~train])
            )
            out[(f, pen)] = (model, loss)
    return context, out


def rr_matrix(
    windows: pd.DataFrame, models: Dict[str, Dict], penalty: float, oof: bool
) -> np.ndarray:
    """
    Build a (window x context) rr matrix, out of fold or from all-chromosome models.

    :param windows: Window features with ``fold``.
    :param models: Context -> fit_context output.
    :param penalty: Chosen penalty.
    :param oof: Whether to use the other-fold models.
    :return: rr matrix (NaN where a feature is missing).
    """
    rr = np.full((len(windows), len(CONTEXTS)), np.nan)
    for context, fits in models.items():
        for f in FOLDS if oof else ("all",):
            rows = np.where(windows.fold == f)[0] if oof else np.arange(len(windows))
            m = fits[(f, penalty)][0]
            z = (windows.iloc[rows][m["cols"]].to_numpy() - m["mean"]) / m["sd"]
            eta0 = m["intercept"]
            p = 1 / (1 + np.exp(-(eta0 + z @ m["coef"])))
            rr[rows, CONTEXTS.index(context)] = p / (1 / (1 + np.exp(-eta0)))
    return rr


def saturating_expected(ev: Evaluator, rr: np.ndarray) -> np.ndarray:
    """
    Get per-transcript expected counts with rr applied before saturation.

    Each cell's expected becomes n_sites * (1 - (1 - q)^rr), q = exp / n_sites, then
    is rescaled per context so the context totals are unchanged.

    :param ev: Evaluator.
    :param rr: (window x context) rr matrix.
    :return: Expected per transcript.
    """
    ok = ~np.isnan(ev.cell_win)
    r = np.ones(len(ev.exp_cell))
    r[ok] = rr[ev.cell_win[ok].astype(int), ev.cell_ctx[ok]]
    r = np.where(np.isnan(r), 1.0, r)
    n = ev.cells.n_sites.to_numpy().astype(float)
    q = np.where(n > 0, np.clip(ev.exp_cell / np.maximum(n, 1), 0, 0.999), 0)
    e = np.where(n > 0, n * (1 - np.power(1 - q, r)), ev.exp_cell)
    num = np.bincount(ev.cell_ctx, weights=e, minlength=32)
    den = np.bincount(ev.cell_ctx, weights=ev.exp_cell, minlength=32)
    e = e * (den / num)[ev.cell_ctx]
    return np.bincount(ev.tx_codes, weights=e)


def score_expected(ev: Evaluator, exp: np.ndarray) -> dict:
    """
    Score per-transcript expected counts with the report's metrics.

    :param ev: Evaluator.
    :param exp: Expected per transcript.
    :return: Core noise, outliers and calibration slope.
    """
    z = (ev.obs - exp) / np.sqrt(exp)
    m = ev.slope_set
    x, y, w = np.log(exp[m] / ev.exp0[m]), np.log(ev.obs[m] / ev.exp0[m]), ev.exp0[m]
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    return {
        "core_noise": ev.noise(exp, ev.core),
        "outliers": int((np.abs(z) > 8).sum()),
        "outliers_no_qc_flag": int(((np.abs(z) > 8) & ev.no_qc_flag).sum()),
        "slope": float(np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2)),
    }


def class_slopes(
    ev: Evaluator, rr: np.ndarray, saturating: bool, mask: np.ndarray = None
) -> dict:
    """
    Get the calibration slope within each context class.

    Per transcript, obs and expected are summed over the class's cells only, and the
    slope is fit over unflagged transcripts with at least 10 expected in the class.

    :param ev: Evaluator.
    :param rr: (window x context) rr matrix.
    :param saturating: Whether rr is applied saturation-aware.
    :param mask: Optional transcripts to use (default all unflagged).
    :return: Class -> slope.
    """
    mask = ev.unflagged if mask is None else mask
    cls = class_of(ev.cells.context)
    ok = ~np.isnan(ev.cell_win)
    r = np.ones(len(ev.exp_cell))
    r[ok] = rr[ev.cell_win[ok].astype(int), ev.cell_ctx[ok]]
    r = np.where(np.isnan(r), 1.0, r)
    if saturating:
        n = ev.cells.n_sites.to_numpy().astype(float)
        q = np.where(n > 0, np.clip(ev.exp_cell / np.maximum(n, 1), 0, 0.999), 0)
        e = np.where(n > 0, n * (1 - np.power(1 - q, r)), ev.exp_cell)
    else:
        e = ev.exp_cell * r
    num = np.bincount(ev.cell_ctx, weights=e, minlength=32)
    den = np.bincount(ev.cell_ctx, weights=ev.exp_cell, minlength=32)
    e = e * (den / num)[ev.cell_ctx]
    out = {}
    obs_cell = ev.cells.obs.to_numpy()
    for c in CLASSES:
        k = cls == c
        o = np.bincount(ev.tx_codes[k], weights=obs_cell[k], minlength=len(ev.tx_names))
        e0 = np.bincount(
            ev.tx_codes[k], weights=ev.exp_cell[k], minlength=len(ev.tx_names)
        )
        e1 = np.bincount(ev.tx_codes[k], weights=e[k], minlength=len(ev.tx_names))
        m = mask & (e0 >= 10) & (o > 0)
        x, y, w = np.log(e1[m] / e0[m]), np.log(o[m] / e0[m]), e0[m]
        xm, ym = np.average(x, weights=w), np.average(y, weights=w)
        out[c] = float(np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2))
    return out


def strength_sweep(
    ev: Evaluator, rr_all: np.ndarray, rr_oof: np.ndarray, saturating: bool
) -> pd.DataFrame:
    """
    Score rr^alpha over ``STRENGTHS``, for all-chromosome and out-of-fold rr.

    :param ev: Evaluator.
    :param rr_all: All-chromosome rr matrix.
    :param rr_oof: Out-of-fold rr matrix.
    :param saturating: Whether rr is applied saturation-aware.
    :return: One row per strength.
    """
    rows = []
    odd, even = ev.core & (ev.tx_fold == "odd"), ev.core & (ev.tx_fold == "even")
    for a in STRENGTHS:
        if saturating:
            e_all = saturating_expected(ev, np.power(rr_all, a))
            e_oof = saturating_expected(ev, np.power(rr_oof, a))
        else:
            e_all, e_oof = ev.expected(rr_all, a), ev.expected(rr_oof, a)
        sc = score_expected(ev, e_all)
        rows.append(
            {
                "alpha": a,
                **sc,
                "oof_noise": ev.noise(e_oof, ev.core),
                "oof_odd": ev.noise(e_oof, odd),
                "oof_even": ev.noise(e_oof, even),
            }
        )
    return pd.DataFrame(rows)


def crossfit(sweep: pd.DataFrame, ev: Evaluator) -> Tuple[float, float, float]:
    """
    Choose alpha on one chromosome half and score it on the other.

    :param sweep: Output of :func:`strength_sweep`.
    :param ev: Evaluator (for the half sizes).
    :return: Alpha chosen on odd, alpha chosen on even, and the cross-fit noise.
    """
    a_odd = sweep.loc[sweep.oof_odd.idxmin()]
    a_even = sweep.loc[sweep.oof_even.idxmin()]
    n_odd = (ev.core & (ev.tx_fold == "odd")).sum()
    n_even = (ev.core & (ev.tx_fold == "even")).sum()
    noise = (a_odd.oof_even * n_even + a_even.oof_odd * n_odd) / (n_odd + n_even)
    return float(a_odd.alpha), float(a_even.alpha), float(noise)


def class_strength(rr: np.ndarray, alphas: dict) -> np.ndarray:
    """
    Raise each context class's rr columns to its own strength.

    :param rr: (window x context) rr matrix.
    :param alphas: Class -> strength.
    :return: Adjusted rr matrix.
    """
    out = rr.copy()
    for i, c in enumerate(CONTEXTS):
        out[:, i] = np.power(rr[:, i], alphas[class_of(pd.Series([c]))[0]])
    return out


def per_class_strength(ev: Evaluator, rr_all: np.ndarray, rr_oof: np.ndarray) -> dict:
    """
    Score rr with one strength per context class, set to that class's calibration slope.

    Out of fold, the strengths are measured on one chromosome half's transcripts and
    scored on the other half.

    :param ev: Evaluator.
    :param rr_all: All-chromosome rr matrix.
    :param rr_oof: Out-of-fold rr matrix.
    :return: Metrics like a strength-sweep row, plus the all-chromosome strengths.
    """
    alphas = {
        c: float(np.clip(a, 0.1, 3.0))
        for c, a in class_slopes(ev, rr_all, False).items()
    }
    sc = score_expected(ev, ev.expected(class_strength(rr_all, alphas)))
    halves = {f: ev.core & (ev.tx_fold == f) for f in FOLDS}
    total, n = 0.0, 0
    for f, other in (("odd", "even"), ("even", "odd")):
        a = {
            c: float(np.clip(v, 0.1, 3.0))
            for c, v in class_slopes(
                ev, rr_oof, False, ev.unflagged & (ev.tx_fold == f)
            ).items()
        }
        e = ev.expected(class_strength(rr_oof, a))
        total += ev.noise(e, halves[other]) * halves[other].sum()
        n += halves[other].sum()
    return {**sc, "oof_noise": total / n, "alphas": alphas}


def main(args: argparse.Namespace) -> None:
    """Fit each background, score it, and write the tables behind the figure."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    windows = ev.windows
    feat_cols = sorted({c for ctx in CONTEXTS for c in allowed_columns(ctx)})
    agg = cell_background(ev).merge(windows[["win", "fold", *feat_cols]], on="win")
    dnm = load_dnms(args.dnm_files, args.twobit)
    keys = set(zip(agg.win, agg.context))
    in_cells = np.array([k in keys for k in zip(dnm.win, dnm.context)])
    logger.info(
        "%d of %d DNMs fall in an exported (window, context) cell.",
        in_cells.sum(),
        len(dnm),
    )
    dnm = dnm[in_cells].merge(windows[["win", "fold", *feat_cols]], on="win")
    logger.info("DNMs per possible synonymous variant: %.4f", len(dnm) / agg.n.sum())

    dnm_by_ctx = {c: g for c, g in dnm.groupby("context")}
    bg_by_ctx = {c: g for c, g in agg.groupby("context")}
    rows, slopes, tuning, sweeps, coefs, rrs = [], [], [], [], [], {}
    none = ev.score(None)
    rows.append(
        {
            "option": "No correction",
            "background": "",
            "application": "",
            "penalty": np.nan,
            "core_noise": none["core_noise"],
            "oof_noise": none["core_noise"],
            "outliers": none["outliers"],
            "outliers_no_qc_flag": none["outliers_no_qc_flag"],
            "slope": np.nan,
        }
    )
    for bg_name in BACKGROUNDS:
        fits = dict(
            Parallel(n_jobs=args.n_jobs)(
                delayed(fit_context)(
                    c,
                    dnm_by_ctx[c][["fold", *allowed_columns(c)]],
                    bg_by_ctx[c][["fold", bg_name, *allowed_columns(c)]],
                    bg_name,
                    PENALTIES,
                )
                for c in CONTEXTS
            )
        )
        loss = {
            p: sum(fits[c][(f, p)][1] for c in fits for f in FOLDS) for p in PENALTIES
        }
        for p, v in loss.items():
            tuning.append({"background": bg_name, "penalty": p, "heldout_log_loss": v})
        pen = min(loss, key=loss.get)
        rr_all = rr_matrix(windows, fits, pen, oof=False)
        rr_oof = rr_matrix(windows, fits, pen, oof=True)
        rrs[bg_name] = rr_all
        for context, f in fits.items():
            m = f[("all", pen)][0]
            for col, b in zip(m["cols"], m["coef"]):
                coefs.append(
                    {
                        "background": bg_name,
                        "context": context,
                        "class": class_of(pd.Series([context]))[0],
                        "feature": col.rsplit("_", 1)[0],
                        "scale": col.rsplit("_", 1)[1],
                        "coef": b,
                    }
                )
        variants = [
            ("multiplicative", rr_all, rr_oof, False),
            ("saturation-aware", rr_all, rr_oof, True),
        ]
        if bg_name in ("absent", "all_sites"):
            cpg = [CONTEXTS.index(c) for c in CONTEXTS if c[1] == "C" and c[2] == "G"]
            no_cpg_all, no_cpg_oof = rr_all.copy(), rr_oof.copy()
            no_cpg_all[:, cpg] = np.nan
            no_cpg_oof[:, cpg] = np.nan
            variants.append(("no rr at CpG", no_cpg_all, no_cpg_oof, False))
        if bg_name in ("absent", "all_sites"):
            pc = per_class_strength(ev, rr_all, rr_oof)
            rows.append(
                {
                    "option": f"{bg_name}, per-class strength",
                    "background": bg_name,
                    "application": "per-class strength",
                    "penalty": pen,
                    "core_noise": pc["core_noise"],
                    "oof_noise": pc["oof_noise"],
                    "outliers": pc["outliers"],
                    "outliers_no_qc_flag": pc["outliers_no_qc_flag"],
                    "slope": pc["slope"],
                    "class_alphas": str(
                        {k: round(v, 2) for k, v in pc["alphas"].items()}
                    ),
                }
            )
            logger.info("%s", rows[-1])
            for c, sl in class_slopes(
                ev, class_strength(rr_all, pc["alphas"]), False
            ).items():
                slopes.append(
                    {
                        "background": bg_name,
                        "application": "per-class strength",
                        "class": c,
                        "slope": sl,
                    }
                )
        for application, r_all, r_oof, sat in variants:
            sw = strength_sweep(ev, r_all, r_oof, sat).assign(
                background=bg_name, application=application
            )
            sweeps.append(sw)
            at1 = sw[sw.alpha == 1.0].iloc[0]
            a_odd, a_even, cf = crossfit(sw, ev)
            row = {
                "option": f"{bg_name}, {application}",
                "background": bg_name,
                "application": application,
                "penalty": pen,
                "core_noise": at1.core_noise,
                "oof_noise": at1.oof_noise,
                "outliers": at1.outliers,
                "outliers_no_qc_flag": at1.outliers_no_qc_flag,
                "slope": at1.slope,
                "alpha_odd": a_odd,
                "alpha_even": a_even,
                "crossfit_noise": cf,
            }
            rows.append(row)
            logger.info("%s", row)
            for c, s in class_slopes(ev, r_all, sat).items():
                slopes.append(
                    {
                        "background": bg_name,
                        "application": application,
                        "class": c,
                        "slope": s,
                    }
                )

    pre = args.out_prefix
    pd.DataFrame(rows).to_csv(f"{pre}.metrics.tsv", sep="\t", index=False)
    pd.DataFrame(slopes).to_csv(f"{pre}.class_slopes.tsv", sep="\t", index=False)
    pd.DataFrame(tuning).to_csv(f"{pre}.penalty_tuning.tsv", sep="\t", index=False)
    pd.concat(sweeps).to_csv(f"{pre}.strength.tsv", sep="\t", index=False)
    pd.DataFrame(coefs).to_csv(f"{pre}.coefficients.tsv", sep="\t", index=False)

    # Diagnostics: what the absent-from-gnomAD rule removes, and how it changes rr.
    win_idx = pd.Series(np.arange(len(windows)), index=windows.win)
    wi = win_idx.reindex(agg.win).to_numpy()
    ci = agg.context.map({c: i for i, c in enumerate(CONTEXTS)}).to_numpy()
    agg["cls"] = class_of(agg.context)
    for name in ("all_sites", "absent"):
        r = rrs[name][wi, ci]
        # Normalize per context over possible sites so rr averages 1.
        s = pd.Series(r * agg.n).groupby(agg.context).transform("sum") / pd.Series(
            np.where(np.isnan(r), 0, agg.n)
        ).groupby(agg.context).transform("sum")
        agg[f"rr_{name}"] = r / s.to_numpy()
    d = agg.dropna(subset=["rr_all_sites", "rr_absent"])
    q_rows, pair_rows = [], []
    for c, g in d.groupby("cls"):
        g = g.assign(b=pd.qcut(g.rr_all_sites.rank(method="first"), 10, labels=False))
        q_bar = g.obs.sum() / g.n.sum()
        for b, h in g.groupby("b"):
            q_rows.append(
                {
                    "class": c,
                    "bin": b,
                    "rr_all_sites": np.average(h.rr_all_sites, weights=h.n),
                    "frac_seen": h.obs.sum() / h.n.sum(),
                    "class_frac_seen": q_bar,
                }
            )
        g = g.assign(
            b=pd.qcut(np.log(g.rr_all_sites).rank(method="first"), 20, labels=False)
        )
        for b, h in g.groupby("b"):
            pair_rows.append(
                {
                    "class": c,
                    "bin": b,
                    "log_rr_all_sites": np.average(np.log(h.rr_all_sites), weights=h.n),
                    "log_rr_absent": np.average(np.log(h.rr_absent), weights=h.n),
                    "class_frac_seen": q_bar,
                }
            )
    pd.DataFrame(q_rows).to_csv(f"{pre}.frac_seen_by_rr.tsv", sep="\t", index=False)
    pd.DataFrame(pair_rows).to_csv(f"{pre}.rr_absent_vs_all.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print(
        pd.DataFrame(slopes)
        .pivot_table(
            index=["background", "application"], columns="class", values="slope"
        )
        .round(2)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, help="Decompressed syn_cells TSV.")
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--dnm-files", nargs="+", required=True)
    parser.add_argument("--twobit", required=True)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
