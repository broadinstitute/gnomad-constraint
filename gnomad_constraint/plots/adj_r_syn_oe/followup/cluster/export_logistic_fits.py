"""Export the per-context rr logistic fits of a pipeline run, for the report figures.

Reads a run's rr training Table and model Table, scores every training site with the
chromosome-split model that didn't see it (``score`` = β·z, in log-odds), and writes
the same tables as ``rr_logistic_fits.py``: ``<prefix>.points.tsv`` (a sample of DNM
and background sites with their scores), ``<prefix>.bins.tsv`` (DNMs and background
sites per score decile, with the share of DNMs if both had equal total weight) and
``<prefix>.contexts.tsv`` (per context: DNMs, the intercept of the curve redrawn
with equal weight, the 5th and 95th percentiles of the background scores, and the
mean score of each class with its standard error), and ``<prefix>.hist.tsv`` (per
context and class, the share of sites in each of ``N_HIST`` score bins).
"""

import argparse

import hail as hl
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from gnomad_constraint.resources.constants import CPG_CONTEXTS
from gnomad_constraint.utils.regional_rate import _chrom_fold, rr_training_ht_to_pandas

N_BINS = 10
N_POINTS = 300
N_HIST = 60


def sigmoid(x: np.ndarray) -> np.ndarray:
    """Get the logistic function of x."""
    return 1 / (1 + np.exp(-x))


def class_of(context: str) -> str:
    """Get the context class of a strand-collapsed trinucleotide context."""
    if context[1] == "A":
        return "A/T"
    return "CpG" if context in CPG_CONTEXTS else "C non-CpG"


def main(args: argparse.Namespace) -> None:
    """Score the training sites out of fold and write the tables."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    df = rr_training_ht_to_pandas(hl.read_table(args.training_ht))
    df["fold"] = df.locus.str.split(":").str[0].map(_chrom_fold)
    models = [
        m
        for m in hl.read_table(args.models_ht).collect()
        if m.cutoff == args.cutoff and m.fold != "all"
    ]
    rng = np.random.default_rng(0)
    bins, points, contexts, hists = [], [], [], []
    for context in sorted({m.context for m in models}):
        ctx = df[df.context == context]
        parts = []
        for m in (m for m in models if m.context == context):
            d = ctx[ctx.fold == m.fold].dropna(subset=list(m.features))
            z = (d[list(m.features)].to_numpy() - np.array(m.mean)) / np.array(m.sd)
            parts.append(d[["y"]].assign(score=z @ np.array(m.coefs)))
        d = pd.concat(parts)
        s1, s0 = d[d.y == 1].score.to_numpy(), d[d.y == 0].score.to_numpy()
        k = len(s1) / len(s0)
        a_bal = brentq(
            lambda v: np.sum(1 - sigmoid(v + s1)) - k * np.sum(sigmoid(v + s0)), -20, 20
        )
        edges = np.quantile(s0, np.linspace(0, 1, N_BINS + 1))
        b1 = np.clip(np.searchsorted(edges, s1, side="right") - 1, 0, N_BINS - 1)
        b0 = np.clip(np.searchsorted(edges, s0, side="right") - 1, 0, N_BINS - 1)
        n1 = np.bincount(b1, minlength=N_BINS)
        n0 = np.bincount(b0, minlength=N_BINS)
        mid = np.bincount(b0, weights=s0, minlength=N_BINS) / n0
        for i in range(N_BINS):
            bins.append(
                {
                    "context": context,
                    "background": "absent",
                    "bin": i,
                    "score": mid[i],
                    "dnms": n1[i],
                    "background_sites": n0[i],
                    "observed_p_equal_weight": n1[i] / (n1[i] + k * n0[i]),
                }
            )
        hist_edges = np.linspace(*np.quantile(s0, [0.005, 0.995]), N_HIST + 1)
        for y, s in ((1, s1), (0, s0)):
            pick = rng.choice(len(s), size=min(N_POINTS, len(s)), replace=False)
            points.append(pd.DataFrame({"context": context, "y": y, "score": s[pick]}))
            clipped = np.clip(s, hist_edges[0], hist_edges[-1])
            hists.append(
                pd.DataFrame(
                    {
                        "context": context,
                        "y": y,
                        "lo": hist_edges[:-1],
                        "hi": hist_edges[1:],
                        "share": np.histogram(clipped, hist_edges)[0] / len(s),
                    }
                )
            )
        contexts.append(
            {
                "context": context,
                "class": class_of(context),
                "dnms": len(s1),
                "intercept_equal_weight": a_bal,
                "score_p05": np.quantile(s0, 0.05),
                "score_p95": np.quantile(s0, 0.95),
                "mean_dnm": s1.mean(),
                "se_dnm": s1.std() / np.sqrt(len(s1)),
                "mean_background": s0.mean(),
                "se_background": s0.std() / np.sqrt(len(s0)),
                "sd_background": s0.std(),
            }
        )
    for name, table in (
        ("bins", pd.DataFrame(bins)),
        ("points", pd.concat(points)),
        ("contexts", pd.DataFrame(contexts)),
        ("hist", pd.concat(hists)),
    ):
        with hl.hadoop_open(f"{args.out_prefix}.{name}.tsv", "w") as f:
            table.to_csv(f, sep="\t", index=False)
    print(pd.DataFrame(contexts).round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-ht", required=True)
    parser.add_argument("--models-ht", required=True)
    parser.add_argument("--cutoff", default="ridge")
    parser.add_argument("--out-prefix", required=True, help="GCS path prefix.")
    main(parser.parse_args())
