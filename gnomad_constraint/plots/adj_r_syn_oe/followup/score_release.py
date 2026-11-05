"""Score the LoF and missense metrics of the downstream constraint runs.

Inputs are the slim release TSVs from ``cluster/export_release_slim.py``. For each
version, on Ensembl MANE Select transcripts off chrX/Y, against the run without rr
("none"):

- LOEUF: Spearman correlation with no rr, genes changing LOEUF decile (by any and by
  2+ deciles), genes moving into or out of the most constrained decile, and genes
  crossing LOEUF 0.6;
- ClinGen haploinsufficiency (score 3) against all other genes: AUC of LOEUF, share in
  the most constrained decile and in the two most constrained, and median LOEUF; the
  same shares for autosomal recessive genes (score 30) and genes where dosage
  sensitivity is unlikely (score 40), which shouldn't look constrained;
- LoF and missense obs/exp by LOEUF decile of the run without rr (so the bins don't
  move between versions), missense-constrained genes (missense z >= 3.09) and pLI >=
  0.9;
- constraint flag counts.
"""

import argparse
import logging
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

LOEUF = "lof.oe_ci.upper"
DECILE = "lof.oe_ci.upper_bin_decile"


def load(path: str) -> pd.DataFrame:
    """Read one version's slim release TSV, keeping Ensembl MANE Select autosomal rows."""
    d = pd.read_csv(path, sep="\t", compression="gzip", na_values="NA")
    d = d[
        d.transcript.str.startswith("ENST")
        & d.mane_select.astype(str).str.lower().eq("true")
        & ~d.chromosome.isin(["chrX", "chrY"])
    ]
    return d.set_index("transcript")


def read_clingen(path: str) -> pd.Series:
    """Get the ClinGen haploinsufficiency score of each gene symbol."""
    c = pd.read_csv(
        path,
        sep="\t",
        skiprows=1,
        header=None,
        names=["chr", "start", "end", "gene", "score"],
    )
    return c.groupby("gene").score.max()


def auc(pos: pd.Series, neg: pd.Series) -> float:
    """Get the AUC of low LOEUF for telling ``pos`` genes from ``neg`` genes."""
    pos, neg = pos.dropna(), neg.dropna()
    return float(mannwhitneyu(-pos, -neg).statistic / (len(pos) * len(neg)))


def main(args: argparse.Namespace) -> None:
    """Score every version."""
    clingen = read_clingen(args.clingen)
    base = load(os.path.join(args.release_dir, "none.tsv.bgz"))
    base_dec = base[DECILE]
    rows, deciles, flags = [], [], []
    for v in args.versions:
        d = load(os.path.join(args.release_dir, f"{v}.tsv.bgz"))
        d = d.join(base[[LOEUF, DECILE]], rsuffix="_none")
        score = d.gene.map(clingen)
        hi, ar, unlikely = score.eq(3), score.eq(30), score.eq(40)
        both = d[[DECILE, f"{DECILE}_none"]].dropna()
        moved = both[DECILE] - both[f"{DECILE}_none"]
        top, top_none = d[DECILE].eq(0), d[f"{DECILE}_none"].eq(0)
        lt, lt_none = d[LOEUF].lt(0.6), d[f"{LOEUF}_none"].lt(0.6)
        row = {
            "version": v,
            "n": len(d),
            "loeuf_spearman_vs_none": spearmanr(
                d[LOEUF], d[f"{LOEUF}_none"], nan_policy="omit"
            ).statistic,
            "decile_changed": float(moved.ne(0).mean()),
            "decile_changed_2plus": float(moved.abs().ge(2).mean()),
            "into_top_decile": int((top & ~top_none).sum()),
            "out_of_top_decile": int((~top & top_none).sum()),
            "into_loeuf_lt_0.6": int((lt & ~lt_none).sum()),
            "out_of_loeuf_lt_0.6": int((~lt & lt_none).sum()),
            "n_clingen_hi": int(hi.sum()),
            "hi_auc": auc(d[LOEUF][hi], d[LOEUF][~hi]),
            "hi_top_decile": float(d[DECILE][hi].eq(0).mean()),
            "hi_top_2_deciles": float(d[DECILE][hi].le(1).mean()),
            "hi_median_loeuf": float(d[LOEUF][hi].median()),
            "ar_top_decile": float(d[DECILE][ar].eq(0).mean()),
            "unlikely_top_decile": float(d[DECILE][unlikely].eq(0).mean()),
            "lof_oe": float(d["lof.obs"].sum() / d["lof.exp"].sum()),
            "mis_oe": float(d["mis.obs"].sum() / d["mis.exp"].sum()),
            "mis_z_ge_3.09": int(d["mis.z_score"].ge(3.09).sum()),
            "pli_ge_0.9": int(d["lof.pLI"].ge(0.9).sum()),
        }
        rows.append(row)
        g = d.groupby(base_dec.reindex(d.index))
        for csq in ("lof", "mis", "syn"):
            oe = g[f"{csq}.obs"].sum() / g[f"{csq}.exp"].sum()
            for dec, val in oe.items():
                deciles.append(
                    {"version": v, "csq": csq, "decile_none": int(dec), "oe": val}
                )
        f = (
            d.constraint_flags.fillna("[]")
            .str.strip("[]")
            .str.replace('"', "")
            .str.split(",")
            .explode()
            .str.strip()
        )
        for flag, n in f[f.ne("")].value_counts().items():
            flags.append({"version": v, "flag": flag, "n": int(n)})
        logger.info(
            "%s: %s",
            v,
            {k: round(x, 4) if isinstance(x, float) else x for k, x in row.items()},
        )
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.scores.tsv", sep="\t", index=False)
    pd.DataFrame(deciles).to_csv(
        f"{args.out_prefix}.deciles.tsv", sep="\t", index=False
    )
    pd.DataFrame(flags).pivot(index="flag", columns="version", values="n").fillna(
        0
    ).astype(int)[args.versions].to_csv(f"{args.out_prefix}.flags.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).set_index("version").T.round(4).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True)
    parser.add_argument("--versions", nargs="+", required=True)
    parser.add_argument("--clingen", required=True)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
