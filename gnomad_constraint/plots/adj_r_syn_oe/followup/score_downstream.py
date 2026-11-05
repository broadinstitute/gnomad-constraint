"""Score the downstream constraint runs from their per-transcript group counts.

Inputs are the ``<version>.group_oe.tsv.bgz`` files from
``cluster/export_group_oe.py``: per transcript, the adj observed, expected (with rr and
any per-variant rates) and unadjusted expected counts of each constraint group, from
the constraint pipeline itself rather than the cell-level approximation.

For each version, on Ensembl MANE Select transcripts off chrX/Y:

- synonymous noise (level-free mean z²) on Ruchit's core set (``unflagged`` in the
  evaluation export, expected >= 20) and on a set that doesn't select on the outcome
  (no release QC gene flag, expected >= 20);
- synonymous outliers (|z| > 8), all and without a QC gene flag;
- synonymous and missense obs/exp by LOEUF decile (v4.1.2 release) and by s_het
  decile (Kar et al. 2026);
- with ``--genes`` (``outlier_drivers.all_genes.tsv``): synonymous obs/exp, relative to
  the core set's, of the control genes in each maternal hotspot region (no QC flag, not
  a histone gene, not an outlier, expected >= 20), and of the 23 excess outliers in
  those regions, with how many of them are still outliers.

The "none" version is the pipeline run without rr; its expected counts equal the
unadjusted ones of every other version.
"""

import argparse
import logging
import os

import numpy as np
import pandas as pd
from outlier_drivers import HISTONE, MATERNAL_HOTSPOTS
from remap_rr_to_shet import read_release_flags, read_shet

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def noise(o: pd.Series, e: pd.Series) -> float:
    """Get the level-free mean z² (expected rescaled to the observed total)."""
    e = e * o.sum() / e.sum()
    return float(np.mean((o - e) ** 2 / e))


def load(path: str, tx: pd.DataFrame) -> pd.DataFrame:
    """Read one version's group counts, keeping Ensembl MANE Select autosomal rows."""
    d = pd.read_csv(path, sep="\t", compression="gzip")
    d = d[
        d.transcript.str.startswith("ENST")
        & d.mane_select.astype(str).str.lower().eq("true")
    ]
    d = d.merge(tx, on="transcript")
    return d[~d.is_xy]


def hotspots(d: pd.DataFrame, genes: pd.DataFrame, lvl: float) -> list:
    """
    Pool synonymous obs/exp over the maternal hotspot regions and their excess outliers.

    :param d: One version's rows (from ``load``), with each row's synonymous ``z``.
    :param genes: ``outlier_drivers.all_genes.tsv``.
    :param lvl: Synonymous obs/exp of the core set, which obs/exp is relative to.
    :return: One row per region and set.
    """
    g = d.merge(
        genes[["transcript", "chromosome", "mid", "outlier", "group"]], on="transcript"
    )
    control = (
        g.no_qc_flag.fillna(False).astype(bool)
        & ~g.gene.fillna("").str.match(HISTONE)
        & ~g.outlier
        & (g.csq_set_syn_unadj_expected >= 20)
    )
    hot = g.outlier & g.group.eq("maternal hotspot region")
    rows = []
    for name, (c, s, e) in MATERNAL_HOTSPOTS.items():
        r = g.chromosome.eq(c) & g.mid.between(s, e)
        for label, m in (("control genes", r & control), ("excess outliers", r & hot)):
            x = g[m]
            rows.append(
                {
                    "region": name,
                    "set": label,
                    "n": len(x),
                    "oe": x.csq_set_syn_observed_variants.sum()
                    / x.csq_set_syn_expected_variants.sum()
                    / lvl,
                    "still_outliers": int(x.z.gt(8).sum()),
                }
            )
    return rows


def main(args: argparse.Namespace) -> None:
    """Score every version."""
    tx = pd.read_csv(args.transcripts, sep="\t", compression="gzip")
    tx["unflagged"] = tx.unflagged.astype(str).str.lower().eq("true")
    tx["is_xy"] = tx.is_xy.astype(str).str.lower().eq("true")
    rel = read_release_flags(args.release_tsv)
    rel["no_qc_flag"] = rel.gene_flags.map(len).eq(0)
    tx = tx.merge(
        rel[["transcript", "gene_id", "no_qc_flag"]], on="transcript", how="left"
    )
    shet = read_shet(args.shet)
    tx = tx.merge(shet, on="gene_id", how="left")
    tx["shet_decile"] = pd.qcut(tx.shet.rank(method="first"), 10, labels=False)
    tx = tx.rename(columns={"decile": "loeuf_decile"})

    genes = pd.read_csv(args.genes, sep="\t") if args.genes else None

    rows, deciles, hot = [], [], []
    for v in args.versions:
        d = load(os.path.join(args.group_oe_dir, f"{v}.group_oe.tsv.bgz"), tx)
        o = d.csq_set_syn_observed_variants
        e = d.csq_set_syn_expected_variants
        e0 = d.csq_set_syn_unadj_expected
        core = d.unflagged & (e0 >= 20)
        free = d.no_qc_flag.fillna(False).astype(bool) & (e0 >= 20)
        lvl = o[core].sum() / e[core].sum()
        z = (o - e * lvl) / np.sqrt(e * lvl)
        rows.append(
            {
                "version": v,
                "n_core": int(core.sum()),
                "core_noise": noise(o[core], e[core]),
                "outcome_free_noise": noise(o[free], e[free]),
                "outliers": int((z.abs() > 8).sum()),
                "outliers_no_qc_flag": int(
                    ((z.abs() > 8) & d.no_qc_flag.fillna(False)).sum()
                ),
                "syn_oe": float(lvl),
                "mis_oe": float(
                    d.csq_set_mis_observed_variants[core].sum()
                    / d.csq_set_mis_expected_variants[core].sum()
                ),
            }
        )
        for binning in ("loeuf_decile", "shet_decile"):
            for csq in ("syn", "mis"):
                g = d.dropna(subset=[binning]).groupby(binning)
                oe = (
                    g[f"csq_set_{csq}_observed_variants"].sum()
                    / g[f"csq_set_{csq}_expected_variants"].sum()
                )
                for dec, val in oe.items():
                    deciles.append(
                        {
                            "version": v,
                            "binning": binning,
                            "csq": csq,
                            "decile": int(dec),
                            "oe": val,
                        }
                    )
        if genes is not None:
            hot += [{"version": v, **h} for h in hotspots(d.assign(z=z), genes, lvl)]
        logger.info(
            "%s: %s",
            v,
            {
                k: round(x, 3) if isinstance(x, float) else x
                for k, x in rows[-1].items()
            },
        )
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.scores.tsv", sep="\t", index=False)
    pd.DataFrame(deciles).to_csv(
        f"{args.out_prefix}.deciles.tsv", sep="\t", index=False
    )
    pd.set_option("display.width", 200)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    if hot:
        h = pd.DataFrame(hot)
        h.to_csv(f"{args.out_prefix}.hotspots.tsv", sep="\t", index=False)
        print(
            h.pivot(
                index=["region", "set"],
                columns="version",
                values=["oe", "still_outliers"],
            )
            .round(3)
            .to_string()
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-oe-dir", required=True)
    parser.add_argument("--versions", nargs="+", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--shet", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--genes", help="outlier_drivers.all_genes.tsv")
    main(parser.parse_args())
