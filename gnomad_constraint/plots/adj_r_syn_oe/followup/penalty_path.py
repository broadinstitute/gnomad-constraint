"""Trace what the ridge penalty does to the rr models, for the report's penalty figure.

The local approximation of ``rr_background_experiments.py`` (Satterstrom DNMs against
the absent-from-gnomAD background built from Ruchit's exported MANE Select cells; a
penalty per DNM). For each penalty on a grid, every context's model is fit on each
chromosome half and scored on the other. Writes:

- ``<prefix>.loss.tsv``: held-out DNM log loss per penalty, summed over contexts and
  both halves (how the penalty is chosen);
- ``<prefix>.coefs.tsv``: one example context's weights per penalty (all chromosomes);
- ``<prefix>.rr_range.tsv``: that context's rr at the 5th and 95th percentile of
  possible sites per penalty (what the penalty does to rr's strength).
"""

import argparse
import logging

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rr_background_experiments import cell_background, fit_context, load_dnms
from rr_logistic_fits import weighted_quantiles
from rr_shrinkage_experiments import CONTEXTS, FOLDS, Evaluator, allowed_columns

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

PENALTIES = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000]
"""Per-DNM penalties traced (0.03 is close to no penalty)."""
BACKGROUND = "absent"


def main(args: argparse.Namespace) -> None:
    """Fit every context at every penalty and write the tables."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    feat_cols = sorted({c for ctx in CONTEXTS for c in allowed_columns(ctx)})
    agg = cell_background(ev).merge(ev.windows[["win", "fold", *feat_cols]], on="win")
    dnm = load_dnms(args.dnm_files, args.twobit)
    keys = set(zip(agg.win, agg.context))
    dnm = dnm[[k in keys for k in zip(dnm.win, dnm.context)]]
    dnm = dnm.merge(ev.windows[["win", "fold", *feat_cols]], on="win")
    fits = dict(
        Parallel(n_jobs=args.n_jobs)(
            delayed(fit_context)(
                c,
                dnm[dnm.context == c][["fold", *allowed_columns(c)]],
                agg[agg.context == c][["fold", BACKGROUND, *allowed_columns(c)]],
                BACKGROUND,
                PENALTIES,
            )
            for c in CONTEXTS
        )
    )
    loss = pd.DataFrame(
        [
            {
                "penalty": p,
                "heldout_log_loss": sum(
                    fits[c][(f, p)][1] for c in fits for f in FOLDS
                ),
            }
            for p in PENALTIES
        ]
    )
    loss.to_csv(f"{args.out_prefix}.loss.tsv", sep="\t", index=False)
    logger.info("Held-out log loss:\n%s", loss.to_string(index=False))

    c = args.context
    bg = agg[(agg.context == c) & (agg[BACKGROUND] > 0)].dropna(
        subset=allowed_columns(c)
    )
    coefs, ranges = [], []
    for p in PENALTIES:
        m = fits[c][("all", p)][0]
        for col, b in zip(m["cols"], m["coef"]):
            coefs.append({"penalty": p, "column": col, "weight": b})
        z = (bg[m["cols"]].to_numpy() - m["mean"]) / m["sd"]
        eta0 = m["intercept"]
        rr = (1 / (1 + np.exp(-(eta0 + z @ m["coef"])))) / (1 / (1 + np.exp(-eta0)))
        lo, hi = weighted_quantiles(rr, bg.n.to_numpy(), np.array([0.05, 0.95]))
        ranges.append({"penalty": p, "rr_p05": lo, "rr_p95": hi})
    pd.DataFrame(coefs).to_csv(f"{args.out_prefix}.coefs.tsv", sep="\t", index=False)
    pd.DataFrame(ranges).to_csv(
        f"{args.out_prefix}.rr_range.tsv", sep="\t", index=False
    )
    logger.info(
        "%s rr range:\n%s", c, pd.DataFrame(ranges).round(3).to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, help="Decompressed syn_cells TSV.")
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--dnm-files", nargs="+", required=True)
    parser.add_argument("--twobit", required=True)
    parser.add_argument("--context", default="CCT")
    parser.add_argument("--n-jobs", type=int, default=6)
    parser.add_argument("--out-prefix", default="penalty_path")
    main(parser.parse_args())
