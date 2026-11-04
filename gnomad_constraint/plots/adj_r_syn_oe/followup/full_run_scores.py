"""Score the rr tables from the full pipeline run on Ruchit's evaluation export.

Inputs:

- The exported synonymous cells annotated with each rr (``cluster/annotate_cells_rr.py``,
  decompressed): columns ``rr_<name>`` for the legacy tables and for each full-run
  table (``<cutoff>_<run>`` and ``<cutoff>_<run>_oof``, run = ``full`` or ``qcexcl``).
- The transcripts table of the export, the v4.1.2 release TSV, the Kar et al. 2026
  s_het table, and the ``rr_cv_metrics.tsv`` of each run.

Each version is scored like the rest of the report: per-gene synonymous noise (level-
free mean z²) on Ruchit's core set, on held-out chromosomes, and on a core set that
doesn't select on the outcome; outliers by release flag; the calibration slope overall
and by context class; noise at each strength; synonymous obs/exp by LOEUF and s_het
decile; and noise on chrX/Y. The held-out DNM fit comes from the cv metrics.
"""

import argparse
import logging
from typing import Dict, Optional

import numpy as np
import pandas as pd
from remap_rr_to_shet import GENE_FAMILIES, read_release_flags, read_shet
from rr_shrinkage_experiments import CONTEXTS, Evaluator

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

STRENGTHS = [0.25, 0.5, 0.625, 0.75, 1.0, 1.25, 1.5]
CLASSES = ["A/T", "C non-CpG", "CpG"]
RUNS = {"full": "all training genes", "qcexcl": "QC-flagged genes left out"}
FITS = {
    "bonferroni": ("bonferroni", 1.0, "Current method (Bonferroni screen, no penalty)"),
    "1e-2_a0.625": ("1e-2", 0.625, "Ruchit's latest pick (1e-2 cutoff, α 0.625)"),
    "ridge": ("ridge", 1.0, "Ridge"),
}


def class_of(context: pd.Series) -> np.ndarray:
    """
    Get the context class of each strand-collapsed trinucleotide context.

    :param context: Contexts (middle base A or C).
    :return: "A/T", "C non-CpG" or "CpG".
    """
    return np.where(
        context.str[1] == "A",
        "A/T",
        np.where(context.str[2] == "G", "CpG", "C non-CpG"),
    )


def cell_rr(ev: Evaluator, col: Optional[str], alpha: float = 1.0) -> np.ndarray:
    """
    Get each cell's rr, raised to ``alpha`` and renormalized per context if alpha != 1.

    :param ev: Evaluator holding the cells.
    :param col: rr column, or None for no correction.
    :param alpha: Strength.
    :return: rr per cell (1 where missing).
    """
    if col is None:
        return np.ones(len(ev.exp_cell))
    r = ev.cells[col].to_numpy(dtype=float)
    r = np.where(np.isnan(r), 1.0, r)
    if alpha != 1.0:
        r = np.power(r, alpha)
        num = np.bincount(ev.cell_ctx, weights=ev.exp_cell * r, minlength=32)
        den = np.bincount(ev.cell_ctx, weights=ev.exp_cell, minlength=32)
        r = r / (num / den)[ev.cell_ctx]
    return r


def tx_expected(ev: Evaluator, r: np.ndarray) -> np.ndarray:
    """Get per-transcript expected counts for per-cell rr ``r``."""
    return np.bincount(ev.tx_codes, weights=ev.exp_cell * r, minlength=len(ev.tx_names))


def slope(obs: np.ndarray, e0: np.ndarray, e1: np.ndarray, mask: np.ndarray) -> float:
    """
    Get the calibration slope: log obs/e0 on log e1/e0, weighted by e0.

    :param obs: Observed per transcript.
    :param e0: Uncorrected expected.
    :param e1: rr-corrected expected.
    :param mask: Transcripts to use.
    :return: Slope (1 = right strength).
    """
    m = mask & (obs > 0) & (e0 > 0)
    x, y, w = np.log(e1[m] / e0[m]), np.log(obs[m] / e0[m]), e0[m]
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2))


def score(
    ev: Evaluator,
    r: np.ndarray,
    r_oof: Optional[np.ndarray],
    masks: Dict[str, np.ndarray],
) -> dict:
    """
    Score one version.

    :param ev: Evaluator.
    :param r: All-chromosome per-cell rr.
    :param r_oof: Held-out-chromosome per-cell rr, or None.
    :param masks: Named transcript masks (``outcome_free``, ``xy``).
    :return: Metrics.
    """
    exp = tx_expected(ev, r)
    z = (ev.obs - exp) / np.sqrt(exp)
    out = {
        "core_noise": ev.noise(exp, ev.core),
        "oof_noise": (
            ev.noise(tx_expected(ev, r_oof), ev.core) if r_oof is not None else np.nan
        ),
        "outcome_free_noise": ev.noise(exp, masks["outcome_free"]),
        "xy_noise": ev.noise(exp, masks["xy"]),
        "outliers": int((np.abs(z) > 8).sum()),
        "outliers_no_qc_flag": int(((np.abs(z) > 8) & ev.no_qc_flag).sum()),
        "slope": slope(ev.obs, ev.exp0, exp, ev.slope_set),
    }
    cls = class_of(ev.cells.context)
    obs_cell = ev.cells.obs.to_numpy()
    for c in CLASSES:
        k = cls == c
        n = len(ev.tx_names)
        o = np.bincount(ev.tx_codes[k], weights=obs_cell[k], minlength=n)
        e0 = np.bincount(ev.tx_codes[k], weights=ev.exp_cell[k], minlength=n)
        e1 = np.bincount(ev.tx_codes[k], weights=(ev.exp_cell * r)[k], minlength=n)
        out[f"slope_{c}"] = slope(o, e0, e1, ev.unflagged & (e0 >= 10))
    return out


def outlier_categories(ev: Evaluator, r: np.ndarray, rel: pd.DataFrame) -> pd.DataFrame:
    """
    Count synonymous outliers by direction and release flag (histone genes split out).

    :param ev: Evaluator.
    :param r: Per-cell rr.
    :param rel: Release rows indexed by transcript with ``gene`` and ``gene_flags``.
    :return: Counts by direction and category.
    """
    exp = tx_expected(ev, r)
    z = (ev.obs - exp) / np.sqrt(exp)
    info = rel.reindex(ev.tx_names)
    flags = info.gene_flags.map(lambda f: f if isinstance(f, list) else [])
    histone = info.gene.fillna("").str.match(GENE_FAMILIES["histone genes"]).to_numpy()
    cat = np.where(
        flags.map(lambda f: "low_exome_mapping_quality" in f),
        "low_exome_mapping_quality",
        np.where(
            flags.map(lambda f: "low_exome_coverage" in f),
            "low_exome_coverage",
            np.where(histone, "histone genes", "unflagged"),
        ),
    )
    d = pd.DataFrame({"z": z, "category": cat})
    d = d[d.z.abs() > 8].assign(
        direction=np.where(d.z[d.z.abs() > 8] > 0, "excess", "deficit")
    )
    return d.groupby(["direction", "category"]).size().rename("n").reset_index()


def deciles(
    ev: Evaluator,
    versions: Dict[str, np.ndarray],
    tx_path: str,
    rel: pd.DataFrame,
    shet_path: str,
) -> pd.DataFrame:
    """
    Sum synonymous obs and expected by LOEUF and s_het decile, as in remap_rr_to_shet.py.

    :param ev: Evaluator.
    :param versions: Version name -> per-cell rr.
    :param tx_path: Transcripts table of the export.
    :param rel: Release rows with ``gene_id`` and ``gene_flags``.
    :param shet_path: Kar et al. 2026 s_het xlsx.
    :return: One row per binning, decile and version.
    """
    tx = pd.read_csv(tx_path, sep="\t", compression="gzip", na_values=["NA"])
    tx = tx[tx.is_xy.astype(str).str.lower().eq("false") & tx.decile.notna()]
    keep = rel[rel.gene_flags.map(len) == 0].reset_index()[["transcript", "gene_id"]]
    tx = tx.merge(keep, on="transcript").merge(read_shet(shet_path), on="gene_id")
    tx["shet"] = pd.qcut(
        tx.shet.rank(method="first", ascending=False), 10, labels=False
    )
    tx["loeuf"] = tx.decile.astype(int)
    idx = (
        pd.Series(np.arange(len(ev.tx_names)), index=ev.tx_names)
        .reindex(tx.transcript)
        .to_numpy()
    )
    rows = []
    for name, r in versions.items():
        exp = tx_expected(ev, r)[idx]
        obs = ev.obs[idx]
        for binning in ("loeuf", "shet"):
            g = (
                pd.DataFrame({"bin": tx[binning].to_numpy(), "obs": obs, "exp": exp})
                .groupby("bin")
                .sum()
            )
            for dec, row in g.iterrows():
                rows.append(
                    {
                        "binning": binning,
                        "decile": dec,
                        "version": name,
                        "n": len(tx),
                        "obs": row.obs,
                        "exp": row.exp,
                        "oe": row.obs / row.exp,
                    }
                )
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    """Score every version and write the tables."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    rel = read_release_flags(args.release_tsv).set_index("transcript")
    tx = (
        pd.read_csv(args.transcripts, sep="\t", compression="gzip")
        .set_index("transcript")
        .reindex(ev.tx_names)
    )
    masks = {
        "outcome_free": ev.no_qc_flag & (ev.exp0 >= 20),
        "xy": tx.is_xy.astype(str).str.lower().eq("true").to_numpy(),
    }

    # syn_adj_r normalized per context over the exported cells, like the pipeline's rr.
    r = ev.cells.rr_syn_adj_r.fillna(1.0).to_numpy()
    norm = np.bincount(
        ev.cell_ctx, weights=ev.exp_cell * r, minlength=32
    ) / np.bincount(ev.cell_ctx, weights=ev.exp_cell, minlength=32)
    ev.cells["rr_syn_adj_r_norm"] = r / norm[ev.cell_ctx]

    versions = {
        "none": ("No correction", None, None, 1.0),
        "adj_r": ("Gnocchi adj_r", "rr_adj_r", None, 1.0),
        "syn_adj_r_norm": (
            "Legacy syn_adj_r (normalized)",
            "rr_syn_adj_r_norm",
            None,
            1.0,
        ),
    }
    for run in RUNS:
        for key, (cutoff, alpha, label) in FITS.items():
            versions[f"{key}_{run}"] = (
                f"{label}, {RUNS[run]}",
                f"rr_{cutoff}_{run}",
                f"rr_{cutoff}_{run}_oof",
                alpha,
            )

    rows, strength, outl, cells_rr = [], [], [], {}
    for name, (label, col, col_oof, alpha) in versions.items():
        r = cell_rr(ev, col, alpha)
        r_oof = cell_rr(ev, col_oof, alpha) if col_oof else None
        cells_rr[name] = r
        sc = score(ev, r, r_oof, masks)
        rows.append({"version": name, "label": label, **sc})
        logger.info(
            "%s: %s",
            name,
            {k: round(v, 3) if isinstance(v, float) else v for k, v in sc.items()},
        )
        if col is not None:
            for a in STRENGTHS:
                ra = cell_rr(ev, col, alpha * a)
                strength.append(
                    {
                        "version": name,
                        "alpha": a,
                        "core_noise": ev.noise(tx_expected(ev, ra), ev.core),
                    }
                )
        outl.append(outlier_categories(ev, r, rel).assign(version=name))

    pre = args.out_prefix
    pd.DataFrame(rows).to_csv(f"{pre}.metrics.tsv", sep="\t", index=False)
    pd.DataFrame(strength).to_csv(f"{pre}.strength.tsv", sep="\t", index=False)
    pd.concat(outl).to_csv(f"{pre}.outliers_by_flag.tsv", sep="\t", index=False)
    deciles(ev, cells_rr, args.transcripts, rel, args.shet).to_csv(
        f"{pre}.deciles.tsv", sep="\t", index=False
    )

    fit = []
    for run, path in (
        ("full", args.cv_metrics_full),
        ("qcexcl", args.cv_metrics_qcexcl),
    ):
        cv = pd.read_csv(path, sep="\t")
        cv["ll_gain"] = (cv.null_log_loss - cv.log_loss) * cv.n_test
        g = (
            cv.groupby("cutoff")
            .agg(
                heldout_ll_gain=("ll_gain", "sum"),
                n_test=("n_test", "sum"),
                n_contexts=("context", "nunique"),
            )
            .reset_index()
        )
        fit.append(
            g.assign(
                run=run, gain_per_1000_dnm_sites=1000 * g.heldout_ll_gain / g.n_test
            )
        )
    pd.concat(fit).to_csv(f"{pre}.dnm_fit.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).drop(columns="label").round(3).to_string(index=False))
    print(pd.concat(fit).round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells", required=True, help="Decompressed annotated cells TSV."
    )
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--shet", required=True)
    parser.add_argument("--cv-metrics-full", required=True)
    parser.add_argument("--cv-metrics-qcexcl", required=True)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
