"""Export what the per-context rr logistic regressions look like, for the report figures.

Each context's model is fit on the Satterstrom DNMs against each training background
(the local approximation of ``rr_background_experiments.py``: ridge on all features,
penalty 10 per DNM). Every site is scored out of fold: by the model fit on the other
half of the chromosomes. The score s = β·z is the model's combined regional score, in
log-odds units, so rr ≈ exp(s).

Writes:

- ``<prefix>.points.tsv``: a sample of DNM sites (y = 1) and background sites (y = 0,
  drawn in proportion to the current background's weights) with their score under the
  current background's model;
- ``<prefix>.bins.tsv``: per context, background and score decile (deciles of possible
  sites under the current background's model): DNMs, possible sites, background weight,
  the observed and fitted share of DNMs among DNMs plus background, the mean rr, the
  DNM rate per possible site relative to the context's, and the observed share if DNMs
  and background had equal total weight;
- ``<prefix>.contexts.tsv``: per context, the DNM count and the intercept of the
  current background's curve redrawn with equal total weight for DNMs and background.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rr_background_experiments import cell_background, class_of, fit_context, load_dnms
from rr_shrinkage_experiments import CONTEXTS, FOLDS, Evaluator, allowed_columns
from scipy.optimize import brentq

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

PENALTY = 10.0
BACKGROUNDS = [
    "absent",
    "absent_context_matched",
    "absent_equal_per_context",
    "all_sites",
]
REFERENCE = "absent"
N_BINS = 10
N_POINTS = 300


def sigmoid(x: np.ndarray) -> np.ndarray:
    """Get the logistic function of x."""
    return 1 / (1 + np.exp(-x))


def oof_eta(fits: dict, x: np.ndarray, fold: np.ndarray) -> tuple:
    """
    Score rows with the model from the other chromosome fold.

    :param fits: fit_context output for one penalty.
    :param x: Features.
    :param fold: Fold of each row.
    :return: (score β·z, intercept β₀) per row.
    """
    s = np.full(len(x), np.nan)
    b0 = np.full(len(x), np.nan)
    for f in FOLDS:
        m = fits[(f, PENALTY)][0]
        k = fold == f
        s[k] = ((x[k] - m["mean"]) / m["sd"]) @ m["coef"]
        b0[k] = m["intercept"]
    return s, b0


def weighted_quantiles(v: np.ndarray, w: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Get weighted quantiles of v."""
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return np.interp(q * cw[-1], cw, v[o])


def summarize_context(
    context: str, dnm: pd.DataFrame, cells: pd.DataFrame, seed: int
) -> tuple:
    """
    Fit one context's models for every background and summarize them.

    :param context: Trinucleotide context.
    :param dnm: DNM rows of the context with features and ``fold``.
    :param cells: Cells of the context with features, ``fold``, ``n`` and the
        background weight columns.
    :param seed: Sampling seed.
    :return: (bins rows, points rows, context row).
    """
    cols = allowed_columns(context)
    d = dnm.dropna(subset=cols)
    a = cells.dropna(subset=cols)
    xd, xa = d[cols].to_numpy(), a[cols].to_numpy()
    fd, fa = d.fold.to_numpy(), a.fold.to_numpy()
    n = a.n.to_numpy()
    scores = {}
    for bg in BACKGROUNDS:
        _, fits = fit_context(context, d, a, bg, [PENALTY])
        scores[bg] = (oof_eta(fits, xd, fd), oof_eta(fits, xa, fa))
    (s_d, _), (s_a, _) = scores[REFERENCE]

    edges = weighted_quantiles(s_a, n, np.linspace(0, 1, N_BINS + 1))
    bin_d = np.clip(np.searchsorted(edges, s_d, side="right") - 1, 0, N_BINS - 1)
    bin_a = np.clip(np.searchsorted(edges, s_a, side="right") - 1, 0, N_BINS - 1)
    n_bin = np.bincount(bin_a, weights=n, minlength=N_BINS)
    d_bin = np.bincount(bin_d, minlength=N_BINS).astype(float)
    s_mid = np.bincount(bin_a, weights=n * s_a, minlength=N_BINS) / n_bin
    dnm_rel = (d_bin / n_bin) / (d_bin.sum() / n_bin.sum())

    rows = []
    for bg in BACKGROUNDS:
        (sd, b0d), (sa, b0a) = scores[bg]
        w = a[bg].to_numpy()
        pd_, pa = sigmoid(b0d + sd), sigmoid(b0a + sa)
        rr = pa / sigmoid(b0a)
        rr = rr / (np.sum(n * rr) / n.sum())
        w_bin = np.bincount(bin_a, weights=w, minlength=N_BINS)
        fit_bin = np.bincount(bin_d, weights=pd_, minlength=N_BINS) + np.bincount(
            bin_a, weights=w * pa, minlength=N_BINS
        )
        for i in range(N_BINS):
            rows.append(
                {
                    "context": context,
                    "background": bg,
                    "bin": i,
                    "score": s_mid[i],
                    "dnms": d_bin[i],
                    "possible": n_bin[i],
                    "background_weight": w_bin[i],
                    "observed_p": d_bin[i] / (d_bin[i] + w_bin[i]),
                    "fitted_p": fit_bin[i] / (d_bin[i] + w_bin[i]),
                    "rr": np.bincount(bin_a, weights=n * rr, minlength=N_BINS)[i]
                    / n_bin[i],
                    "dnm_rate_relative": dnm_rel[i],
                }
            )
    bins = pd.DataFrame(rows)

    # The current background's curve with equal total weight for DNMs and background:
    # the same slope, a new intercept.
    w = a[REFERENCE].to_numpy()
    k = len(d) / w.sum()
    a_bal = brentq(
        lambda v: np.sum(1 - sigmoid(v + s_d)) - k * np.sum(w * sigmoid(v + s_a)),
        -20,
        20,
    )
    ref = bins.background == REFERENCE
    w_bin = bins.loc[ref, "background_weight"].to_numpy()
    bins.loc[ref, "observed_p_equal_weight"] = d_bin / (d_bin + k * w_bin)

    rng = np.random.default_rng(seed)
    pick_d = rng.choice(len(d), size=min(N_POINTS, len(d)), replace=False)
    pick_a = rng.choice(len(a), size=N_POINTS, replace=True, p=w / w.sum())
    points = pd.concat(
        [
            pd.DataFrame({"context": context, "y": 1, "score": s_d[pick_d]}),
            pd.DataFrame({"context": context, "y": 0, "score": s_a[pick_a]}),
        ]
    )
    ctx = {
        "context": context,
        "class": class_of(pd.Series([context]))[0],
        "dnms": len(d),
        "intercept_equal_weight": a_bal,
        "score_p05": weighted_quantiles(s_a, n, np.array([0.05]))[0],
        "score_p95": weighted_quantiles(s_a, n, np.array([0.95]))[0],
    }
    return bins, points, ctx


def main(args: argparse.Namespace) -> None:
    """Fit, score and export the per-context logistic fits."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    feat_cols = sorted({c for ctx in CONTEXTS for c in allowed_columns(ctx)})
    agg = cell_background(ev).merge(ev.windows[["win", "fold", *feat_cols]], on="win")
    dnm = load_dnms(args.dnm_files, args.twobit)
    dnm = dnm.merge(agg[["win", "context"]], on=["win", "context"]).merge(
        ev.windows[["win", "fold", *feat_cols]], on="win"
    )
    keep = ["fold", "n", *BACKGROUNDS]
    out = Parallel(n_jobs=args.n_jobs)(
        delayed(summarize_context)(
            c,
            dnm[dnm.context == c][["fold", *allowed_columns(c)]],
            agg[agg.context == c][[*keep, *allowed_columns(c)]],
            i,
        )
        for i, c in enumerate(CONTEXTS)
    )
    bins = pd.concat([o[0] for o in out])
    bins.to_csv(f"{args.out_prefix}.bins.tsv", sep="\t", index=False)
    pd.concat([o[1] for o in out]).to_csv(
        f"{args.out_prefix}.points.tsv", sep="\t", index=False
    )
    ctx = pd.DataFrame([o[2] for o in out])
    ctx.to_csv(f"{args.out_prefix}.contexts.tsv", sep="\t", index=False)
    logger.info("Contexts:\n%s", ctx.round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, help="Decompressed syn_cells TSV.")
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--dnm-files", nargs="+", required=True)
    parser.add_argument("--twobit", required=True)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
