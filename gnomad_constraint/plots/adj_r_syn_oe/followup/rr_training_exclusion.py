"""Test removing poorly covered or mapped genes from rr training, locally.

Same data, background (absent from gnomAD, the current rule) and ridge fit as
``rr_background_experiments.py``. Training DNMs and background cells in the removed
genes are dropped before fitting; the evaluation is unchanged (rr is still applied to
every gene, and the scores use the same transcripts). A (window, context) cell is
removed if any transcript covering it is removed. Removal sets:

- ``none``: train on everything (the current pipeline).
- ``qc_flagged``: genes with a release QC gene flag (``low_exome_coverage`` or
  ``low_exome_mapping_quality``).
- ``qc_flagged_or_segdup``: also genes with segdup fraction >= 0.5 or AN90 < 0.5 (gene
  quality metrics only, not outlier status, so the outcome isn't used).

For each set it reports how many DNMs and background sites are dropped and their
DNM obs/exp (DNMs against the context-rate expectation), which shows whether the DNM
cohorts miss calls there.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rr_background_experiments import (
    PENALTIES,
    cell_background,
    class_slopes,
    crossfit,
    fit_context,
    load_dnms,
    per_class_strength,
    rr_matrix,
    strength_sweep,
)
from rr_shrinkage_experiments import CONTEXTS, FOLDS, Evaluator, allowed_columns

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

BACKGROUND = "absent"


def removed_transcripts(release_tsv: str) -> dict:
    """
    Get the transcripts in each removal set.

    :param release_tsv: v4.1.2 release TSV.
    :return: Removal set -> set of transcripts.
    """
    rel = pd.read_csv(
        release_tsv,
        sep="\t",
        compression="gzip",
        na_values=["NA"],
        usecols=[
            "transcript",
            "gene_flags",
            "gene_quality_metrics.exome_prop_segdup",
            "gene_quality_metrics.exome_prop_bp_AN90",
        ],
    )
    rel = rel[rel.transcript.str.startswith("ENST")].drop_duplicates("transcript")
    qc = rel.gene_flags.fillna("[]").str.strip("[]").str.replace('"', "").str.len() > 0
    segdup = (rel["gene_quality_metrics.exome_prop_segdup"] >= 0.5) | (
        rel["gene_quality_metrics.exome_prop_bp_AN90"] < 0.5
    )
    return {
        "none": set(),
        "qc_flagged": set(rel.transcript[qc]),
        "qc_flagged_or_segdup": set(rel.transcript[qc | segdup]),
    }


def main(args: argparse.Namespace) -> None:
    """Fit and score the rr with each removal set."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    windows = ev.windows
    feat_cols = sorted({c for ctx in CONTEXTS for c in allowed_columns(ctx)})
    agg = cell_background(ev).merge(windows[["win", "fold", *feat_cols]], on="win")
    dnm = load_dnms(args.dnm_files, args.twobit)
    keys = set(zip(agg.win, agg.context))
    dnm = dnm[[k in keys for k in zip(dnm.win, dnm.context)]].merge(
        windows[["win", "fold", *feat_cols]], on="win"
    )
    # Context-rate DNM expectation per (window, context), for the obs/exp of
    # removed cells.
    rate = dnm.groupby("context").size() / agg.groupby("context").n.sum()
    agg["dnm_exp"] = agg.n * agg.context.map(rate)
    dnm_counts = dnm.groupby(["win", "context"]).size().rename("dnm")
    agg = agg.join(dnm_counts, on=["win", "context"]).fillna({"dnm": 0})

    rows, slopes = [], []
    for name, tx in removed_transcripts(args.release_tsv).items():
        cells = ev.cells[ev.cells.transcript.isin(tx)]
        tx = set(cells.transcript)
        drop = set(zip(cells.win, cells.context))
        bad = np.array([k in drop for k in zip(agg.win, agg.context)])
        bad_dnm = np.array([k in drop for k in zip(dnm.win, dnm.context)])
        d = agg[bad]
        logger.info(
            "%s: %d exported transcripts; drops %d of %d DNMs and %.1f%% of the background; "
            "DNM obs/exp in the dropped cells %.2f.",
            name,
            len(tx),
            bad_dnm.sum(),
            len(dnm),
            100 * d[BACKGROUND].sum() / agg[BACKGROUND].sum(),
            d.dnm.sum() / d.dnm_exp.sum() if len(d) else np.nan,
        )
        train_bg, train_dnm = agg[~bad], dnm[~bad_dnm]
        dnm_by_ctx = {c: g for c, g in train_dnm.groupby("context")}
        bg_by_ctx = {c: g for c, g in train_bg.groupby("context")}
        fits = dict(
            Parallel(n_jobs=args.n_jobs)(
                delayed(fit_context)(
                    c,
                    dnm_by_ctx[c][["fold", *allowed_columns(c)]],
                    bg_by_ctx[c][["fold", BACKGROUND, *allowed_columns(c)]],
                    BACKGROUND,
                    PENALTIES,
                )
                for c in CONTEXTS
            )
        )
        loss = {
            p: sum(fits[c][(f, p)][1] for c in fits for f in FOLDS) for p in PENALTIES
        }
        pen = min(loss, key=loss.get)
        rr_all = rr_matrix(windows, fits, pen, oof=False)
        rr_oof = rr_matrix(windows, fits, pen, oof=True)
        sw = strength_sweep(ev, rr_all, rr_oof, False)
        at1 = sw[sw.alpha == 1.0].iloc[0]
        pc = per_class_strength(ev, rr_all, rr_oof)
        row = {
            "removed": name,
            "n_transcripts": len(tx),
            "dnms_dropped": int(bad_dnm.sum()),
            "background_dropped_pct": 100 * d[BACKGROUND].sum() / agg[BACKGROUND].sum(),
            "dropped_dnm_oe": d.dnm.sum() / d.dnm_exp.sum() if len(d) else np.nan,
            "penalty": pen,
            "core_noise": at1.core_noise,
            "oof_noise": at1.oof_noise,
            "outliers": at1.outliers,
            "outliers_no_qc_flag": at1.outliers_no_qc_flag,
            "slope": at1.slope,
            "crossfit_noise": crossfit(sw, ev)[2],
            "per_class_oof_noise": pc["oof_noise"],
        }
        rows.append(row)
        logger.info("%s", row)
        for c, s in class_slopes(ev, rr_all, False).items():
            slopes.append({"removed": name, "class": c, "slope": s})
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.metrics.tsv", sep="\t", index=False)
    pd.DataFrame(slopes).to_csv(
        f"{args.out_prefix}.class_slopes.tsv", sep="\t", index=False
    )
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print(
        pd.DataFrame(slopes)
        .pivot_table(index="removed", columns="class", values="slope")
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
