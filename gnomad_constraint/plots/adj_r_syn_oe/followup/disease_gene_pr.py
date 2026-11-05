"""Precision-recall of LOEUF for NDD and HI genes, for each version of the full run.

The same tests as the gnomAD v4 flagship preprint (Guez et al. 2026), with its gene
lists and rules from
[atgu/gnomAD_v4_flagship_paper](https://github.com/atgu/gnomAD_v4_flagship_paper/tree/47ad465e2a80454066740b5c7952aff296f40961),
copied into ``data/ext_flagship/``:

- **NDD genes** (Fig. 4d and Extended Data Fig. 7b; ``Figure_4/figure_4.R`` and
  ``scripts/functions_figure4.R``): positives are the 598 NDD genes in ``ndd.txt``
  (DDG2P monoallelic, confirmed or probable, as of 2020-11-27) and negatives every
  other gene. Haploinsufficient genes are left out entirely (the severe 2016 curated
  list plus the 77 genes added in ``load_gene_lists_figure4``, ``hi_added_genes.txt``),
  as are the 44 genes in ``ces_genes_list.tsv``;
- **HI genes** (Extended Data Fig. 5b; ``Figure_3/figure_3.R`` and
  ``scripts/compute_loeuf_by_lof_filter.R``): positives are the 90 genes in
  ``hi_genes_union.tsv`` and negatives every other gene, with nothing left out.

For each, genes are ranked by LOEUF (lower = more constrained), and only genes with a
LOEUF in every version are used, so each version is scored on the same genes. Two gene
sets: all genes, and genes with fewer than 50 expected pLoF variants (the preprint's
Fig. 4d rule, applied here to both tests), using the run without rr's expected count so
the set doesn't change between versions.

The preprint's LOEUF is a newer LOEUF (its own pLoF filter and adj_r), so its AUPRC
values aren't directly comparable with these. One row per gene: the Ensembl MANE
Select transcript of each version's slim release TSV
(``cluster/export_release_slim.py``). The run without rr reproduces the released
v4.1.1 LOEUF exactly.

AUPRC is average precision (scikit-learn), close to but not the same as the PRROC
``auc.integral`` the preprint uses. Each version's difference from the run without
rr gets a 95% interval from a paired bootstrap over genes.

Writes ``<prefix>.auprc.tsv`` (AUPRC and the bootstrap interval of its difference
from no rr, per test, gene set and version) and ``<prefix>.curves.tsv`` (the curves,
one point per positive gene reached).
"""

import argparse
import logging
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

VERSIONS = [
    "none",
    "bonferroni_full",
    "bonferroni_qcexcl",
    "ridge_full",
    "ridge_qcexcl",
    "roulette_ridge",
    "homebuilt_ridge",
]
LOEUF = "lof.oe_ci.upper"
EXP_CUTOFF = 50


def read_list(path: str, header: bool = False) -> set:
    """Read a one-column gene list."""
    return set(pd.read_csv(path, sep="\t", header=0 if header else None).iloc[:, 0])


def load_versions(release_dir: str) -> pd.DataFrame:
    """Get each version's LOEUF on Ensembl MANE Select transcripts, one row per gene."""
    out = None
    for v in VERSIONS:
        d = pd.read_csv(
            os.path.join(release_dir, f"{v}.tsv.bgz"),
            sep="\t",
            compression="gzip",
            na_values="NA",
            usecols=["gene", "gene_id", "transcript", "mane_select", "lof.exp", LOEUF],
        )
        d = d[
            d.transcript.str.startswith("ENST")
            & d.mane_select.astype(str).str.lower().eq("true")
        ]
        d = d.rename(columns={LOEUF: v, "lof.exp": "exp_none"})
        if out is None:
            out = d[["gene", "gene_id", "transcript", "exp_none", v]]
        else:
            out = out.merge(d[["transcript", v]], on="transcript")
    return out[np.isfinite(out[VERSIONS]).all(axis=1)]


def score(
    g: pd.DataFrame,
    positive: pd.Series,
    test: str,
    rng: np.random.Generator,
    n_boot: int,
) -> tuple:
    """Score each version on all genes and on short genes."""
    rows, curves = [], []
    for name, m in (
        ("all genes", np.ones(len(g), dtype=bool)),
        (f"expected pLoF < {EXP_CUTOFF}", (g.exp_none < EXP_CUTOFF).to_numpy()),
    ):
        y = positive.to_numpy()[m]
        scores = {v: -g[v].to_numpy()[m] for v in VERSIONS}
        ap = {v: average_precision_score(y, s) for v, s in scores.items()}
        boot = {v: [] for v in VERSIONS}
        for _ in range(n_boot):
            i = rng.integers(0, len(y), len(y))
            base = average_precision_score(y[i], scores["none"][i])
            for v in VERSIONS:
                boot[v].append(average_precision_score(y[i], scores[v][i]) - base)
        for v in VERSIONS:
            lo, hi = np.percentile(boot[v], [2.5, 97.5])
            rows.append(
                {
                    "test": test,
                    "gene_set": name,
                    "version": v,
                    "n_genes": len(y),
                    "n_positive": int(y.sum()),
                    "auprc": ap[v],
                    "diff_vs_none": ap[v] - ap["none"],
                    "diff_lo": lo,
                    "diff_hi": hi,
                }
            )
            # Keep one point per recall level: the precision just as each positive
            # gene is reached (sklearn orders recall from 1 down to 0).
            p, r, _ = precision_recall_curve(y, scores[v])
            keep = np.r_[np.diff(r) != 0, True]
            p, r = p[keep], r[keep]
            curves.append(
                pd.DataFrame(
                    {
                        "test": test,
                        "gene_set": name,
                        "version": v,
                        "recall": r,
                        "precision": p,
                    }
                )
            )
    return rows, curves


def main(args: argparse.Namespace) -> None:
    """Score each version for both tests and write the tables."""
    d = args.lists
    g = load_versions(args.release_dir)
    rng = np.random.default_rng(args.seed)

    ndd = read_list(os.path.join(d, "ndd.txt"), header=True)
    hi_severe = read_list(os.path.join(d, "haploinsufficiency_severe_curated_2016.tsv"))
    hi_severe |= read_list(os.path.join(d, "hi_added_genes.txt"))
    ces = set(pd.read_csv(os.path.join(d, "ces_genes_list.tsv"), sep="\t").gene_id)
    g_ndd = g[~g.gene.isin(hi_severe) & ~g.gene_id.isin(ces)]

    hi = read_list(os.path.join(d, "hi_genes_union.tsv"), header=True)

    rows, curves = [], []
    for test, sub, pos in (
        ("NDD", g_ndd, g_ndd.gene.isin(ndd)),
        ("HI", g, g.gene.isin(hi)),
    ):
        logger.info("%s: %d genes, %d positive.", test, len(sub), pos.sum())
        r, c = score(sub, pos, test, rng, args.n_boot)
        rows += r
        curves += c
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}.auprc.tsv", sep="\t", index=False)
    pd.concat(curves).to_csv(f"{args.out_prefix}.curves.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    logger.info("\n%s", out.round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", default="data/cluster/release_slim")
    parser.add_argument("--lists", default="data/ext_flagship")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-prefix", default="disease_gene_pr")
    main(parser.parse_args())
