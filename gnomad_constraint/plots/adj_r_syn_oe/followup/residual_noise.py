"""Characterize the per-gene synonymous noise left after the ridge rr.

Uses Ruchit Panchal's rr evaluation export (cells with an ``rr_recommended`` column
written by ``rr_shrinkage_experiments.py --experiments recommended``), the window
features, and the v4.1.2 release TSV (gene quality metrics and flags).

Questions:

1. Does the leftover noise grow with a gene's expected count (a gene-wide rate
   error, which rr could address) or not (site-level effects rr can't see)?
2. How much of it can gene-level covariates predict out of fold: the 13 regional
   features (nonlinearly), the release gene quality metrics, or both?
3. Which genes stay noisiest?

Genes with any release flag (QC gene flags or constraint outlier flags) are left out:
the release already handles them.
"""

import argparse
import logging
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

QC_COLS = [
    "gene_quality_metrics.exome_prop_bp_AN90",
    "gene_quality_metrics.exome_mean_AS_MQ",
    "gene_quality_metrics.exome_prop_segdup",
    "gene_quality_metrics.exome_prop_LCR",
]
GENE_COLS = ["cds_length", "num_coding_exons"]


def load(args: argparse.Namespace) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build one row per transcript with counts, regional feature means and covariates.

    :param args: Parsed arguments.
    :return: Per-transcript DataFrame and the regional feature column names.
    """
    cells = pd.read_csv(
        args.cells,
        sep="\t",
        usecols=["transcript", "win", "exp", "obs", "n_sites", "rr_recommended"],
    )
    cells["rr_recommended"] = cells.rr_recommended.fillna(1.0)
    cells["exp_rr"] = cells.exp * cells.rr_recommended
    p = (cells.exp / cells.n_sites).clip(upper=1)
    cells["floor_var"] = cells.n_sites * p * (1 - p)

    tx = cells.groupby("transcript").agg(
        obs=("obs", "sum"),
        exp=("exp", "sum"),
        exp_rr=("exp_rr", "sum"),
        floor_var=("floor_var", "sum"),
    )
    # Expected-weighted mean of each regional feature over the transcript's windows.
    tw = cells.groupby(["transcript", "win"], as_index=False).exp.sum()
    windows = pd.read_csv(args.windows, sep="\t", compression="gzip")
    feat_cols = [c for c in windows.columns if c not in ("win", "contig")]
    tw = tw.merge(windows[["win", *feat_cols]], on="win", how="left")
    num = tw[feat_cols].multiply(tw.exp, axis=0).assign(transcript=tw.transcript)
    den = (
        tw[feat_cols].notna().multiply(tw.exp, axis=0).assign(transcript=tw.transcript)
    )
    tx = tx.join(
        num.groupby("transcript")[feat_cols].sum()
        / den.groupby("transcript")[feat_cols].sum()
    )

    tinfo = pd.read_csv(args.transcripts, sep="\t", compression="gzip")
    tx = tx.join(tinfo.set_index("transcript")[["chrom", "is_xy", "decile"]])
    rel = pd.read_csv(
        args.release_tsv,
        sep="\t",
        compression="gzip",
        na_values=["NA"],
        usecols=[
            "gene",
            "gene_id",
            "transcript",
            "gene_flags",
            "constraint_flags",
            *QC_COLS,
            *GENE_COLS,
        ],
    )
    rel = rel[rel.transcript.str.startswith("ENST")].drop_duplicates("transcript")
    tx = tx.join(rel.set_index("transcript"))
    tx["qc_flagged"] = tx.gene_flags.fillna("[]").str.strip("[]").str.len() > 0
    if args.gtex:
        tx = tx.join(read_gtex(args.gtex), on="gene_id")
    tx["fold"] = np.where(
        pd.to_numeric(tx.chrom.str.replace("chr", "").replace({"X": "1", "Y": "2"})) % 2
        == 1,
        "odd",
        "even",
    )
    return tx, feat_cols


EXPR_COLS = ["expr_mean_log2tpm", "expr_max_log2tpm", "expr_n_tissues_tpm1", "expr_tau"]


def read_gtex(path: str) -> pd.DataFrame:
    """
    Summarize GTEx v8 median TPM per tissue into gene-level expression covariates.

    :param path: GTEx gene median TPM GCT file.
    :return: DataFrame indexed by unversioned gene ID with ``EXPR_COLS``.
    """
    g = pd.read_csv(path, sep="\t", skiprows=2, compression="gzip")
    g.index = g.Name.str.split(".").str[0]
    tpm = g.drop(columns=["Name", "Description"])
    x = np.log2(tpm + 1)
    xmax = x.max(axis=1)
    out = pd.DataFrame(
        {
            "expr_mean_log2tpm": x.mean(axis=1),
            "expr_max_log2tpm": xmax,
            "expr_n_tissues_tpm1": (tpm > 1).sum(axis=1),
            # Tissue specificity: 0 = expressed evenly, 1 = in one tissue only.
            "expr_tau": (1 - x.div(xmax.replace(0, np.nan), axis=0)).sum(axis=1)
            / (x.shape[1] - 1),
        }
    )
    return out[~out.index.duplicated()]


def noise(obs: np.ndarray, exp: np.ndarray) -> float:
    """
    Get level-free mean z^2.

    :param obs: Observed counts.
    :param exp: Expected counts.
    :return: Mean z^2 after rescaling expected so total obs/exp = 1.
    """
    exp = exp * obs.sum() / exp.sum()
    return float(np.mean((obs - exp) ** 2 / exp))


def oof_multiplier(d: pd.DataFrame, cols: list) -> np.ndarray:
    """
    Predict a per-gene multiplicative correction to expected, out of fold.

    A gradient-boosted regression of log((obs + 0.5) / exp_rr) on ``cols``, weighted
    by expected, trained on the other chromosome half.

    :param d: Per-transcript rows.
    :param cols: Covariate columns.
    :return: Multiplier per row.
    """
    y = np.log((d.obs + 0.5) / d.exp_rr)
    pred = np.zeros(len(d))
    for fold in ("odd", "even"):
        tr, te = (d.fold != fold).to_numpy(), (d.fold == fold).to_numpy()
        model = HistGradientBoostingRegressor(
            max_iter=300,
            learning_rate=0.05,
            max_leaf_nodes=15,
            min_samples_leaf=50,
            l2_regularization=1.0,
            random_state=0,
        )
        model.fit(d.loc[tr, cols], y[tr], sample_weight=d.exp_rr[tr])
        pred[te] = model.predict(d.loc[te, cols])
    return np.exp(pred)


def main(args: argparse.Namespace) -> None:
    """Run the three analyses and write their tables."""
    tx, feat_cols = load(args)
    any_flag = tx.qc_flagged | tx.constraint_flags.fillna("").str.contains("outlier")
    d = tx[~any_flag & (tx.exp >= 20)].copy()
    logger.info("%d transcripts without any release flag and expected >= 20.", len(d))

    # 1. Noise by gene size: chance floor, no correction, ridge rr.
    s0 = d.obs.sum() / d.exp.sum()
    s1 = d.obs.sum() / d.exp_rr.sum()
    d["z2_none"] = (d.obs - d.exp * s0) ** 2 / (d.exp * s0)
    d["z2_rr"] = (d.obs - d.exp_rr * s1) ** 2 / (d.exp_rr * s1)
    d["z_rr"] = (d.obs - d.exp_rr * s1) / np.sqrt(d.exp_rr * s1)
    d["floor"] = d.floor_var / d.exp
    bins = [20, 50, 100, 200, 400, 800, np.inf]
    d["size_bin"] = pd.cut(d.exp, bins, right=False)
    by_size = (
        d.groupby("size_bin", observed=True)
        .agg(
            n=("obs", "size"),
            median_exp=("exp", "median"),
            chance=("floor", "mean"),
            none=("z2_none", "mean"),
            ridge_rr=("z2_rr", "mean"),
        )
        .reset_index()
    )
    by_size["size_bin"] = by_size.size_bin.astype(str)
    by_size.to_csv(f"{args.out_prefix}.by_size.tsv", sep="\t", index=False)
    for col in ("z2_none", "z2_rr"):
        b, a = np.polyfit(d.exp, d[col], 1)
        logger.info(
            "%s = %.2f + %.5f * exp (chance floor %.2f)", col, a, b, d.floor.mean()
        )
    excess = d.z2_rr - d.floor
    logger.info(
        "Leftover excess noise from genes with z > 0 (more variants than expected): %.0f%%",
        100 * excess[d.z_rr > 0].sum() / excess.sum(),
    )
    ranked = excess.sort_values(ascending=False)
    for frac in (0.01, 0.05):
        k = int(len(ranked) * frac)
        logger.info(
            "Top %.0f%% of genes carry %.0f%% of the leftover excess noise.",
            100 * frac,
            100 * ranked.head(k).sum() / ranked.sum(),
        )
    histone = d.gene.fillna("").str.match(r"^(H1-|H2AC|H2BC|H3C|H4C|H2AW|H3-|H4-)")
    logger.info(
        "Histone cluster genes: %d genes, %.1f%% of the leftover excess noise.",
        histone.sum(),
        100 * excess[histone].sum() / excess.sum(),
    )

    # 2. How much can gene-level covariates predict, out of fold?
    rows = [
        {
            "covariates": "none (ridge rr only)",
            "noise": noise(d.obs.to_numpy(), d.exp_rr.to_numpy()),
        },
        {
            "covariates": "no correction",
            "noise": noise(d.obs.to_numpy(), d.exp.to_numpy()),
        },
    ]
    sets = {
        "13 regional features, nonlinear": feat_cols,
        "gene quality metrics": QC_COLS,
        "gene size": GENE_COLS,
    }
    if args.gtex:
        sets["expression (GTEx)"] = EXPR_COLS
    sets["all of the above"] = [c for cols in sets.values() for c in cols]
    for label, cols in sets.items():
        mult = oof_multiplier(d, cols)
        rows.append(
            {
                "covariates": label,
                "noise": noise(d.obs.to_numpy(), (d.exp_rr * mult).to_numpy()),
            }
        )
        logger.info("%s: %s", label, rows[-1])
    pd.DataFrame(rows).to_csv(
        f"{args.out_prefix}.covariates.tsv", sep="\t", index=False
    )

    # 3. Noisiest unflagged genes.
    top = d.reindex(d.z_rr.abs().sort_values(ascending=False).index).head(40)
    top = top.assign(oe_rr=top.obs / (top.exp_rr * s1))[
        [
            "gene",
            "chrom",
            "obs",
            "exp_rr",
            "oe_rr",
            "z_rr",
            "decile",
            *QC_COLS,
            "constraint_flags",
        ]
    ]
    top.to_csv(f"{args.out_prefix}.top_genes.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(by_size.round(3).to_string(index=False))
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print(top.round(3).head(25).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--gtex", help="GTEx v8 gene median TPM GCT (optional).")
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
