"""Compare the ridge rr with Ruchit's recommended cutoff + strength options, locally.

Fits rr models with Ruchit's code (``gnomad_constraint.utils.regional_rate`` on branch
``jg/rr_ridge_fit``) on his exported training sites, and scores all-chromosome and
out-of-fold rr on his evaluation export with two core sets: his current one
(unflagged, expected >= 20), and the one his 2026-10-27 follow-up recommends, which
doesn't select on the outcome (expected >= 20, release QC gene flags removed).
"""

import argparse

import numpy as np
import pandas as pd
from rr_shrinkage_experiments import CONTEXTS, Evaluator

from gnomad_constraint.utils import regional_rate as rr

OPTIONS = [
    # (label, DNM set, cutoff token, strength alpha)
    ("Old DNMs, Bonferroni, α 1 (reference)", "17k", "bonferroni", 1.0),
    ("New DNMs, 1e-2, α 0.625 (Ruchit's recommendation)", "Satterstrom", "1e-2", 0.625),
    ("New DNMs, 1e-2, α 0.75", "Satterstrom", "1e-2", 0.75),
    ("Old DNMs, ridge, α 1", "17k", "ridge", 1.0),
    ("New DNMs, ridge, α 1", "Satterstrom", "ridge", 1.0),
]


def rr_matrix(ev: Evaluator, models: list, cutoff: str, oof: bool) -> np.ndarray:
    """
    Build a (window x context) rr matrix from fitted models.

    :param ev: Evaluator holding the window features.
    :param models: Model dictionaries from ``fit_rr_models``.
    :param cutoff: Cutoff token to use.
    :param oof: Whether to use the chromosome-split models (out of fold).
    :return: rr matrix (NaN = no rr).
    """
    mat = np.full((len(ev.windows), len(CONTEXTS)), np.nan)
    for m in models:
        if m["cutoff"] != cutoff:
            continue
        if oof and m["fold"] in ("odd", "even"):
            idx = np.where(ev.windows.fold == m["fold"])[0]
        elif not oof and m["fold"] == "all":
            idx = np.arange(len(ev.windows))
        else:
            continue
        w = ev.windows.iloc[idx]
        mat[idx, CONTEXTS.index(m["context"])] = rr.predict_rr_model(m, w) / m["ave"]
    return mat


def main(args: argparse.Namespace) -> None:
    """Fit, score and write one row per option and model type."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    core_b = ev.no_qc_flag & (ev.exp0 >= 20)
    models = {}
    for name, sub in (("17k", "v1_ddd_asd"), ("Satterstrom", "satterstrom")):
        df = pd.read_csv(
            f"{args.training_dir}/{sub}/rr_training.tsv.bgz",
            sep="\t",
            compression="gzip",
            na_values=["NA"],
        )
        df = df.sort_values(["y", "row_idx"], ascending=[False, True], kind="stable")
        screen = pd.read_csv(
            f"{args.training_dir}/{sub}/gnomad.v4.1.1.rr_feature_screen.tsv", sep="\t"
        )
        cutoffs = sorted({c for _, s, c, _ in OPTIONS if s == name})
        models[name], _ = rr.fit_rr_models(
            df.reset_index(drop=True), screen, cutoffs, cv_chrom_split=True
        )

    rows = [
        {
            "option": "No correction",
            "models": "",
            "core_noise": ev.score(None)["core_noise"],
            "core_noise_no_outcome_selection": ev.noise(ev.exp0, core_b),
            "outliers": 270,
        }
    ]
    for label, dnm_set, cutoff, alpha in OPTIONS:
        for oof in (False, True):
            mat = rr_matrix(ev, models[dnm_set], cutoff, oof)
            sc = ev.score(mat, alpha)
            rows.append(
                {
                    "option": label,
                    "models": "out-of-fold" if oof else "all chromosomes",
                    "core_noise": sc["core_noise"],
                    "core_noise_no_outcome_selection": ev.noise(
                        ev.expected(mat, alpha), core_b
                    ),
                    "outliers": sc["outliers"],
                    "outliers_no_qc_flag": sc["outliers_no_qc_flag"],
                    "slope": sc["slope"],
                }
            )
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print("core set without outcome selection:", int(core_b.sum()), "transcripts")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--training-dir", required=True)
    parser.add_argument("--out", required=True)
    main(parser.parse_args())
