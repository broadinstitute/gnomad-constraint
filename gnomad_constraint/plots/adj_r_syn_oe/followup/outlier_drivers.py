"""Test what drives the synonymous outliers left after the ridge rr.

The outliers are the transcripts without a QC gene flag or histone gene symbol with
|z| > 8 under the ridge rr (``--rr-col`` of the cells: ``rr_recommended`` from
``rr_shrinkage_experiments.py --experiments recommended``, or the full run's
``rr_ridge_qcexcl`` from ``cluster/annotate_cells_rr.py``).

Questions:

1. Do de novo mutations (DNMs) show the same excess? DNMs measure the mutation rate
   directly, in different samples, so a shared excess means the mutation model is
   too low for these genes rather than a gnomAD artifact.
2. Do DNMs track synonymous obs/exp across all genes, not just the outliers?
3. Are the outliers in the regions with elevated maternal C>G DNM rates reported by
   Jónsson et al. 2017 (Nature 549:519-522): the first 20 Mb of 8p, 2p16-22, 7p, 9p,
   16p and 16q?
4. Is the excess the same in every genetic ancestry group and downsampling?
5. Is it spread across the gene or concentrated, which contexts carry it, and how do
   the outliers differ from genes of the same size?

Expected counts are rescaled so that genes without a flag have obs/exp = 1, for both
gnomAD and DNMs.
"""

import argparse
import json
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

HISTONE = r"^(H1-|H2AC|H2BC|H3C|H4C|H2AW|H2BW|H3-|H4-)"
HISTONE_CLUSTER_6P22 = ("chr6", 25_000_000, 28_500_000)

MATERNAL_HOTSPOTS = {
    "8p, first 20 Mb": ("chr8", 0, 20_000_000),
    "2p16-22": ("chr2", 31_800_000, 61_000_000),
    "7p": ("chr7", 0, 60_100_000),
    "9p": ("chr9", 0, 43_000_000),
    "16p": ("chr16", 0, 36_800_000),
    "16q": ("chr16", 36_800_000, 90_400_000),
}
"""Regions with elevated maternal C>G DNM rates (Jónsson et al. 2017), in GRCh38
using UCSC hg38 cytoband boundaries."""

REPEAT_FAMILIES = r"^(MUC\d|OR\d|PSG\d)"
GROUPS = ["maternal hotspot region", "repeat, segdup or low coverage", "other"]

DOWNSAMPLINGS = [
    10,
    100,
    500,
    1000,
    2000,
    2884,
    5000,
    10000,
    13068,
    16740,
    19850,
    20000,
    22362,
    26572,
    30198,
    43129,
    50000,
    100000,
    200000,
    500000,
    556006,
]
"""Sample size at each index of the release ``gen_anc_obs`` / ``gen_anc_exp`` arrays:
the standard v4 downsamplings plus each group's full size, sorted (the union used by
``gnomad_constraint.utils.constraint``)."""
STANDARD_DOWNSAMPLINGS = {
    500,
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
    200000,
    500000,
}
ANCESTRIES = ["global", "afr", "amr", "eas", "nfe", "sas"]

QC_COLS = [
    "gene_quality_metrics.exome_prop_bp_AN90",
    "gene_quality_metrics.exome_mean_AS_MQ",
    "gene_quality_metrics.exome_prop_segdup",
]


def poisson_ci(k: float, e: float) -> tuple:
    """
    Get a 95% CI for obs/exp with a Poisson count.

    :param k: Observed count.
    :param e: Expected count.
    :return: Lower and upper bound of obs/exp.
    """
    lo, hi = stats.poisson.interval(0.95, k)
    return lo / e, hi / e


def load_cells(path: str, rr_col: str = "rr_recommended") -> pd.DataFrame:
    """
    Read the per (transcript, window, context) cells with the ridge rr.

    :param path: Cells TSV with ``rr_col``.
    :param rr_col: The rr column to use.
    :return: Cells with ``exp_rr`` and a context class.
    """
    cells = pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "transcript",
            "win",
            "context",
            "exp",
            "obs",
            "n_sites",
            rr_col,
        ],
    )
    cells["rr"] = cells[rr_col].fillna(1.0)
    cells["exp_rr"] = cells.exp * cells.rr
    ctx = cells.context
    cells["cls"] = np.where(
        ctx.str[1] == "A", "A/T", np.where(ctx.str[2] == "G", "CpG", "C non-CpG")
    )
    return cells


def load_release(path: str) -> pd.DataFrame:
    """
    Read the release TSV, one row per Ensembl transcript.

    :param path: v4.1.2 constraint metrics TSV.
    :return: Release rows indexed by transcript.
    """
    rel = pd.read_csv(
        path, sep="\t", compression="gzip", na_values=["NA"], low_memory=False
    )
    rel = rel[rel.transcript.str.startswith("ENST")].drop_duplicates("transcript")
    return rel.set_index("transcript")


def add_dnms(cells: pd.DataFrame, training: str) -> pd.DataFrame:
    """
    Add observed and expected synonymous DNMs to each cell.

    Expected DNMs are n_sites times the exome-wide DNM rate of the context, times rr.
    A DNM in a window and context shared by overlapping transcripts is split by
    n_sites.

    :param cells: Cells from :func:`load_cells`.
    :param training: rr training TSV (rows with y = 1 are DNMs).
    :return: Cells with ``dnm`` and ``dnm_exp`` columns.
    """
    dn = pd.read_csv(
        training, sep="\t", compression="gzip", usecols=["locus", "context", "y"]
    )
    dn = dn[dn.y == 1]
    chrom = dn.locus.str.split(":").str[0]
    start = dn.locus.str.split(":").str[1].astype(int) // 1000 * 1000
    dn["win"] = chrom + "-" + start.astype(str) + "-" + (start + 1000).astype(str)
    counts = dn.groupby(["win", "context"]).size().rename("dnm_wc").reset_index()
    c = cells.merge(counts, on=["win", "context"], how="left").fillna({"dnm_wc": 0})
    logger.info(
        "%d DNMs, %.1f%% in an exported cell.",
        len(dn),
        100 * c.drop_duplicates(["win", "context"]).dnm_wc.sum() / len(dn),
    )
    c["dnm"] = (
        c.dnm_wc * c.n_sites / c.groupby(["win", "context"]).n_sites.transform("sum")
    )
    sites = c.groupby(["win", "context"]).agg(
        n=("n_sites", "max"), d=("dnm_wc", "first")
    )
    rate = (
        sites.groupby(level="context").d.sum() / sites.groupby(level="context").n.sum()
    )
    c["dnm_exp"] = c.n_sites * c.context.map(rate) * c.rr
    return c


def build_genes(cells: pd.DataFrame, rel: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate cells to transcripts and add release annotations and groups.

    :param cells: Cells with DNM columns.
    :param rel: Release rows.
    :return: One row per transcript.
    """
    tx = cells.groupby("transcript").agg(
        obs=("obs", "sum"),
        exp=("exp", "sum"),
        exp_rr=("exp_rr", "sum"),
        dnm=("dnm", "sum"),
        dnm_exp=("dnm_exp", "sum"),
    )
    cols = [
        "gene",
        "gene_id",
        "chromosome",
        "start_position",
        "end_position",
        "gene_flags",
        "constraint_flags",
        "mis.z_raw",
        *QC_COLS,
    ]
    tx = tx.join(rel[cols])
    tx["mid"] = (tx.start_position + tx.end_position) / 2
    tx["qc"] = tx.gene_flags.fillna("[]").str.strip("[]").str.len() > 0
    tx["histone"] = tx.gene.fillna("").str.match(HISTONE)
    c, s, e = HISTONE_CLUSTER_6P22
    tx["histone_6p22"] = tx.histone & (tx.chromosome == c) & tx.mid.between(s, e)
    tx["z"] = (tx.obs - tx.exp_rr) / np.sqrt(tx.exp_rr)
    tx["outlier"] = ~tx.qc & ~tx.histone & (tx.z.abs() > 8)
    tx["control"] = ~tx.qc & ~tx.histone & ~tx.outlier & (tx.exp_rr >= 20)
    ctrl = tx[tx.control]
    tx["exp_g"] = tx.exp_rr * ctrl.obs.sum() / ctrl.exp_rr.sum()
    tx["exp_none"] = tx.exp * ctrl.obs.sum() / ctrl.exp.sum()
    tx["dnm_exp"] = tx.dnm_exp * ctrl.dnm.sum() / ctrl.dnm_exp.sum()
    tx["oe"] = tx.obs / tx.exp_g
    tx["dnm_oe"] = tx.dnm / tx.dnm_exp

    tx["region"] = ""
    for name, (c, s, e) in MATERNAL_HOTSPOTS.items():
        tx.loc[
            (tx.chromosome == c) & tx.mid.between(s, e, inclusive="left"), "region"
        ] = name
    repeat = (
        (tx["gene_quality_metrics.exome_prop_segdup"] >= 0.5)
        | (tx["gene_quality_metrics.exome_prop_bp_AN90"] < 0.5)
        | tx.gene.fillna("").str.match(REPEAT_FAMILIES)
    )
    tx["group"] = np.where(
        repeat, GROUPS[1], np.where(tx.region != "", GROUPS[0], GROUPS[2])
    )
    logger.info(
        "%d outliers (%d with too many variants).",
        tx.outlier.sum(),
        (tx.outlier & (tx.z > 0)).sum(),
    )
    return tx


def size_matched(tx: pd.DataFrame, k: int = 50) -> Dict[str, pd.Index]:
    """
    Pick the ``k`` controls closest in log expected count to each excess outlier.

    :param tx: Transcript rows.
    :param k: Controls per outlier.
    :return: Outlier transcript -> matched control transcripts.
    """
    ctrl = tx[tx.control]
    lc = np.log(ctrl.exp_g.to_numpy())
    return {
        t: ctrl.index[np.argsort(np.abs(lc - np.log(r.exp_g)))[:k]]
        for t, r in tx[tx.outlier & (tx.z > 0)].iterrows()
    }


def pooled(d: pd.DataFrame, label: str) -> dict:
    """
    Pool gnomAD and DNM obs/exp over a set of transcripts.

    :param d: Transcript rows.
    :param label: Row label.
    :return: Summary row.
    """
    lo, hi = poisson_ci(d.dnm.sum(), d.dnm_exp.sum())
    return {
        "set": label,
        "n": len(d),
        "gnomad_oe": d.obs.sum() / d.exp_g.sum(),
        "gnomad_oe_no_rr": (
            d.obs.sum() / d.exp_none.sum() if "exp_none" in d else np.nan
        ),
        "dnm_obs": d.dnm.sum(),
        "dnm_exp": d.dnm_exp.sum(),
        "dnm_oe": d.dnm.sum() / d.dnm_exp.sum(),
        "dnm_lo": lo,
        "dnm_hi": hi,
    }


def dnm_slope(d: pd.DataFrame) -> tuple:
    """
    Fit DNM count ~ Poisson(dnm_exp * (1 + b * (gnomAD obs/exp - 1))).

    :param d: Transcript rows.
    :return: Slope b and its standard error.
    """

    def nll(b):
        mu = d.dnm_exp * np.clip(1 + b * (d.oe - 1), 1e-6, None)
        return float((mu - d.dnm * np.log(mu)).sum())

    b = minimize_scalar(nll, bounds=(-2, 3), method="bounded").x
    h = 1e-3
    return b, 1 / np.sqrt((nll(b + h) - 2 * nll(b) + nll(b - h)) / h**2)


def neighbourhood_oe(tx: pd.DataFrame, width: float = 3e6) -> pd.Series:
    """
    Get pooled obs/exp of control genes within ``width`` of each gene, excluding it.

    :param tx: Transcript rows.
    :param width: Half-width in bp.
    :return: Neighbourhood obs/exp per transcript.
    """
    out = pd.Series(np.nan, index=tx.index)
    for _, g in tx.groupby("chromosome"):
        c = g[g.control]
        mid, obs, exp = c.mid.to_numpy(), c.obs.to_numpy(), c.exp_g.to_numpy()
        for t, m in zip(g.index, g.mid):
            near = (np.abs(mid - m) <= width) & (c.index != t)
            if near.any():
                out[t] = obs[near].sum() / exp[near].sum()
    return out


def within_gene(cells: pd.DataFrame, tx: pd.DataFrame) -> pd.DataFrame:
    """
    Measure how evenly each gene's excess is spread over its 1-kb windows.

    :param cells: Cells.
    :param tx: Transcript rows.
    :return: Dispersion (chi-square per window after scaling expected to the gene's own
        obs/exp; about 1 if the excess is a gene-wide factor) and the fraction of
        windows above expected.
    """
    w = cells.groupby(["transcript", "win"])[["obs", "exp_rr"]].sum().reset_index()
    w = w.join((tx.obs / tx.exp_rr).rename("k"), on="transcript")
    w = w[w.exp_rr > 0]
    w["chi"] = (w.obs - w.exp_rr * w.k) ** 2 / (w.exp_rr * w.k)
    w["up"] = w.obs > w.exp_rr
    g = w.groupby("transcript").agg(
        chi=("chi", "sum"), n_windows=("chi", "size"), frac_windows_up=("up", "mean")
    )
    g["dispersion"] = g.chi / (g.n_windows - 1).clip(lower=1)
    return g[["dispersion", "n_windows", "frac_windows_up"]]


def parse_arrays(rel: pd.DataFrame, index: pd.Index, col: str) -> np.ndarray:
    """
    Parse a release array column into a matrix padded with NaN.

    :param rel: Release rows.
    :param index: Transcripts to keep, in order.
    :param col: Column holding JSON arrays.
    :return: Matrix (transcripts x downsamplings).
    """
    m = np.full((len(index), len(DOWNSAMPLINGS)), np.nan)
    for i, v in enumerate(rel[col].reindex(index)):
        if isinstance(v, str):
            a = json.loads(v)
            m[i, : len(a)] = a
    return m


def by_downsampling(rel: pd.DataFrame, tx: pd.DataFrame) -> pd.DataFrame:
    """
    Get the excess outliers' pooled obs/exp, relative to controls, per downsampling.

    Release expected counts carry no rr, so each gene's expected is scaled by its mean
    ridge rr.

    :param rel: Release rows.
    :param tx: Transcript rows.
    :return: Rows of ancestry, sample size and relative obs/exp.
    """
    rows = []
    rr = (tx.exp_rr / tx.exp).to_numpy()[:, None]
    ex, c = (tx.outlier & (tx.z > 0)).to_numpy(), tx.control.to_numpy()
    for anc in ANCESTRIES:
        o = parse_arrays(rel, tx.index, f"syn.gen_anc_obs.{anc}")
        e = parse_arrays(rel, tx.index, f"syn.gen_anc_exp.{anc}") * rr
        for i, n in enumerate(DOWNSAMPLINGS):
            if (
                n not in STANDARD_DOWNSAMPLINGS
                or np.isnan(e[c, i]).all()
                or np.nansum(e[ex, i]) == 0
            ):
                continue
            ratio = (np.nansum(o[ex, i]) / np.nansum(e[ex, i])) / (
                np.nansum(o[c, i]) / np.nansum(e[c, i])
            )
            lo, hi = poisson_ci(
                np.nansum(o[ex, i]),
                np.nansum(e[ex, i]) * np.nansum(o[c, i]) / np.nansum(e[c, i]),
            )
            rows.append(
                {
                    "ancestry": anc,
                    "n_samples": n,
                    "rel_oe": ratio,
                    "lo": lo,
                    "hi": hi,
                    "outlier_exp": np.nansum(e[ex, i]),
                }
            )
    return pd.DataFrame(rows)


def covariates(
    cells: pd.DataFrame,
    tx: pd.DataFrame,
    match: Dict[str, pd.Index],
    windows: Optional[str],
    gtex: Optional[str],
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Compare the excess outliers with size-matched controls on gene covariates.

    :param cells: Cells.
    :param tx: Transcript rows.
    :param match: Matched controls per outlier.
    :param windows: Window feature TSV (optional).
    :param gtex: GTEx v8 gene median TPM GCT (optional).
    :return: One row per covariate (mean percentile of the outliers among their matched
        controls, 0.5 = no difference, and a Wilcoxon p-value), and the neighbourhood
        obs/exp per transcript.
    """
    d = tx.copy()
    d["rr_mean"] = d.exp_rr / d.exp
    d["nb_oe_3Mb"] = neighbourhood_oe(d)
    cols: List[str] = [*QC_COLS, "rr_mean", "nb_oe_3Mb"]
    if windows:
        win = pd.read_csv(windows, sep="\t", compression="gzip")
        feats = [f for f in win.columns if f.endswith("_100k")]
        tw = (
            cells.groupby(["transcript", "win"], as_index=False)
            .exp.sum()
            .merge(win[["win", *feats]], on="win", how="left")
        )
        num = tw[feats].multiply(tw.exp, axis=0).groupby(tw.transcript).sum()
        den = tw[feats].notna().multiply(tw.exp, axis=0).groupby(tw.transcript).sum()
        d = d.join(num / den)
        cols += feats
    if gtex:
        g = pd.read_csv(gtex, sep="\t", skiprows=2, compression="gzip")
        g.index = g.Name.str.split(".").str[0]
        g = g[~g.index.duplicated()].drop(columns=["Name", "Description"])
        d = d.join(
            np.log2(g + 1)
            .agg(["mean"], axis=1)
            .rename(columns={"mean": "expr_mean_log2tpm"})
            .assign(testis_log2tpm=np.log2(g.Testis + 1)),
            on="gene_id",
        )
        cols += ["expr_mean_log2tpm", "testis_log2tpm"]
    rows = []
    for col in cols:
        pct = []
        for t, ref_idx in match.items():
            ref = d.loc[ref_idx, col].dropna()
            if pd.notna(d.at[t, col]) and len(ref) >= 10:
                pct.append(
                    (ref < d.at[t, col]).mean() + 0.5 * (ref == d.at[t, col]).mean()
                )
        pct = np.array(pct)
        rows.append(
            {
                "covariate": col,
                "n": len(pct),
                "mean_percentile": pct.mean(),
                "outlier_median": d.loc[list(match), col].median(),
                "control_median": d.loc[
                    np.concatenate(list(match.values())), col
                ].median(),
                "p": stats.wilcoxon(pct - 0.5).pvalue,
            }
        )
    return pd.DataFrame(rows).sort_values("p"), d.nb_oe_3Mb


def main(args: argparse.Namespace) -> None:
    """Run the analyses and write their tables."""
    rel = load_release(args.release_tsv)
    cells = add_dnms(load_cells(args.cells, args.rr_col), args.dnm_training)
    tx = build_genes(cells, rel)
    ex = tx.outlier & (tx.z > 0)
    match = size_matched(tx)
    prefix = args.out_prefix

    # 1. DNMs in the outlier groups.
    rows = [pooled(tx[ex & (tx.group == g)], g) for g in GROUPS]
    rows += [
        pooled(tx[ex], "all excess outliers"),
        pooled(tx[tx.histone_6p22], "histone genes, chr6p22 cluster"),
        pooled(
            tx.loc[np.unique(np.concatenate(list(match.values())))],
            "size-matched controls",
        ),
    ]
    groups = pd.DataFrame(rows)
    groups.to_csv(f"{prefix}.groups.tsv", sep="\t", index=False)

    # Robustness: non-CpG contexts only (no methylation in the DNM expectation).
    nc = (
        cells[cells.cls != "CpG"]
        .groupby("transcript")[["obs", "exp_rr", "dnm", "dnm_exp"]]
        .sum()
    )
    nc = nc.join(tx[["control", "group", "histone_6p22"]]).assign(outlier=ex)
    ctrl = nc[nc.control]
    nc["exp_g"] = nc.exp_rr * ctrl.obs.sum() / ctrl.exp_rr.sum()
    nc["dnm_exp"] = nc.dnm_exp * ctrl.dnm.sum() / ctrl.dnm_exp.sum()
    for label, m in (
        ("excess outliers", nc.outlier),
        ("histone chr6p22", nc.histone_6p22),
    ):
        r = pooled(nc[m], label)
        logger.info(
            "Non-CpG only, %s: gnomAD %.2f, DNM %.2f [%.2f-%.2f]",
            label,
            r["gnomad_oe"],
            r["dnm_oe"],
            r["dnm_lo"],
            r["dnm_hi"],
        )

    # 2. DNMs across all genes, by gnomAD obs/exp decile.
    a = tx[(tx.control | tx.outlier) & (tx.exp_g >= 100)].copy()
    a["decile"] = pd.qcut(a.oe, 10, labels=False)
    dec = pd.DataFrame([pooled(d, str(i)) for i, d in a.groupby("decile")])
    dec.to_csv(f"{prefix}.dnm_by_oe_decile.tsv", sep="\t", index=False)
    b, se = dnm_slope(a)
    w = a.exp_g
    v_total = np.average((a.oe - 1) ** 2, weights=w)
    v_chance = 0.7 * np.average(1 / a.exp_g, weights=w)
    logger.info(
        "DNM slope on gnomAD obs/exp (expected >= 100): %.2f +/- %.2f; %.2f if all spread above chance were mutation rate.",
        b,
        se,
        1 - v_chance / v_total,
    )

    # 3. Maternal hotspot regions.
    in_region = (ex & (tx.region != "")).sum()
    frac = (tx.control & (tx.region != "") & (tx.exp_g >= 300)).sum() / (
        tx.control & (tx.exp_g >= 300)
    ).sum()
    logger.info(
        "%d of %d excess outliers in maternal hotspot regions; %.1f expected (binomial p = %.1e).",
        in_region,
        ex.sum(),
        ex.sum() * frac,
        stats.binomtest(in_region, ex.sum(), frac, alternative="greater").pvalue,
    )
    regions = pd.DataFrame(
        [pooled(tx[tx.control & (tx.region == r)], r) for r in MATERNAL_HOTSPOTS]
        + [
            pooled(tx[tx.control & (tx.region != "")], "all six regions"),
            pooled(tx[tx.control & (tx.region == "")], "rest of the genome"),
        ]
    )
    regions.to_csv(f"{prefix}.regions.tsv", sep="\t", index=False)

    # 4. Downsamplings and genetic ancestry groups.
    by_downsampling(rel, tx).to_csv(f"{prefix}.downsampling.tsv", sep="\t", index=False)

    # 5. Within-gene spread, contexts, covariates.
    tx = tx.join(within_gene(cells, tx))
    cls = cells.groupby(["transcript", "cls"])[["obs", "exp_rr"]].sum().unstack()
    base = cls["obs"][tx.control].sum() / cls["exp_rr"][tx.control].sum()
    tx = tx.join((cls["obs"] / cls["exp_rr"] / base).add_prefix("oe_"))
    logger.info(
        "Excess outliers, obs/exp by context class: %s",
        (cls["obs"][ex].sum() / cls["exp_rr"][ex].sum() / base).round(2).to_dict(),
    )
    cov, nb = covariates(cells, tx, match, args.windows, args.gtex)
    cov.to_csv(f"{prefix}.covariates.tsv", sep="\t", index=False)
    tx["nb_oe_3Mb"] = nb
    m = list(match)
    allm = np.concatenate(list(match.values()))
    for col in ("dispersion", "frac_windows_up"):
        logger.info(
            "%s: outliers %.2f, size-matched controls %.2f",
            col,
            tx.loc[m, col].median(),
            tx.loc[allm, col].median(),
        )

    genes = tx[tx.outlier].sort_values("z", ascending=False)[
        [
            "gene",
            "chromosome",
            "start_position",
            "group",
            "region",
            "obs",
            "exp_rr",
            "exp_g",
            "oe",
            "z",
            "dnm",
            "dnm_exp",
            "dnm_oe",
            "nb_oe_3Mb",
            "dispersion",
            "frac_windows_up",
            "oe_A/T",
            "oe_C non-CpG",
            "oe_CpG",
            *QC_COLS,
            "constraint_flags",
            "mis.z_raw",
        ]
    ]
    genes.to_csv(f"{prefix}.genes.tsv", sep="\t")
    tx[~tx.qc][
        ["gene", "chromosome", "mid", "z", "oe", "outlier", "group", "histone_6p22"]
    ].to_csv(f"{prefix}.all_genes.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(groups.round(3).to_string(index=False))
    print(regions.round(3).to_string(index=False))
    print(cov.round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, help="Cells TSV with --rr-col.")
    parser.add_argument(
        "--rr-col",
        default="rr_recommended",
        help="rr column of the cells, e.g. rr_ridge_qcexcl in the full run's cells.",
    )
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument(
        "--dnm-training", required=True, help="rr_training.tsv.bgz of the DNM set."
    )
    parser.add_argument(
        "--windows", help="Window feature TSV (optional, for covariates)."
    )
    parser.add_argument("--gtex", help="GTEx v8 gene median TPM GCT (optional).")
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
