"""Check the synonymous excess in the chr6p22 histone cluster.

1. Excess by downsampling and genetic ancestry group: pooled syn obs/exp of the
   cluster's genes relative to genes without a QC flag, from the release
   ``gen_anc_obs`` / ``gen_anc_exp`` arrays (no rr). A higher mutation rate gives the
   full ratio while few sites are seen and a smaller one as sites saturate; a fixed
   set of common artifact variants gives a ratio that keeps falling with sample size.
2. What a histone factor learned from DNMs would do to gnomAD: multiply expected by
   the cluster's DNM obs/exp, either directly or before saturation
   (n_sites * (1 - (1 - q)^r), q = exp / n_sites), overall or by CpG status.
"""

import argparse
import json

import numpy as np
import pandas as pd
from outlier_drivers import (
    DOWNSAMPLINGS,
    HISTONE,
    HISTONE_CLUSTER_6P22,
    STANDARD_DOWNSAMPLINGS,
)

ANCESTRIES = ["global", "afr", "amr", "eas", "nfe", "sas"]


def parse(rel: pd.DataFrame, col: str) -> np.ndarray:
    """
    Parse a release array column into a matrix padded with NaN.

    :param rel: Release rows.
    :param col: Column of JSON arrays.
    :return: Matrix (rows x downsamplings).
    """
    m = np.full((len(rel), len(DOWNSAMPLINGS)), np.nan)
    for i, v in enumerate(rel[col]):
        if isinstance(v, str):
            a = json.loads(v)
            m[i, : len(a)] = a
    return m


def main(args: argparse.Namespace) -> None:
    """Write the downsampling and DNM-factor tables."""
    rel = pd.read_csv(
        args.release_tsv,
        sep="\t",
        compression="gzip",
        na_values=["NA"],
        low_memory=False,
    )
    rel = rel[
        rel.transcript.str.startswith("ENST") & (rel.mane_select == True)
    ].drop_duplicates("transcript")
    mid = (rel.start_position + rel.end_position) / 2
    c, s, e = HISTONE_CLUSTER_6P22
    histone = rel.gene.fillna("").str.match(HISTONE)
    cluster = histone & (rel.chromosome == c) & mid.between(s, e)
    qc = rel.gene_flags.fillna("[]").str.strip("[]").str.len() > 0
    ctrl = ~qc & ~histone & (rel["syn.exp"] >= 20)

    rows = []
    for anc in ANCESTRIES:
        o, x = parse(rel, f"syn.gen_anc_obs.{anc}"), parse(
            rel, f"syn.gen_anc_exp.{anc}"
        )
        h, k = cluster.to_numpy(), ctrl.to_numpy()
        for i, n in enumerate(DOWNSAMPLINGS):
            if n not in STANDARD_DOWNSAMPLINGS or np.nansum(x[h, i]) == 0:
                continue
            base = np.nansum(o[k, i]) / np.nansum(x[k, i])
            rows.append(
                {
                    "ancestry": anc,
                    "n_samples": n,
                    "rel_oe": np.nansum(o[h, i]) / np.nansum(x[h, i]) / base,
                    "exp": np.nansum(x[h, i]),
                }
            )
    pd.DataFrame(rows).to_csv(
        f"{args.out_prefix}.downsampling.tsv", sep="\t", index=False
    )

    cells = pd.read_csv(args.cells, sep="\t")
    k = cells.obs.sum() / cells.exp.sum()
    d = cells[cells.transcript.isin(rel.transcript[cluster])].copy()
    cpg = ((d.context.str[1] == "C") & (d.context.str[2] == "G")).to_numpy()
    q = (d.exp / d.n_sites).clip(upper=0.999).to_numpy()
    out = []
    for label, r_non, r_cpg, sat in (
        ("no correction", 1.0, 1.0, False),
        (f"× {args.dnm_oe:g}, multiplied", args.dnm_oe, args.dnm_oe, False),
        (f"× {args.dnm_oe:g}, before saturation", args.dnm_oe, args.dnm_oe, True),
        (
            f"× {args.dnm_oe_noncpg:g} non-CpG, × {args.dnm_oe_cpg:g} CpG, before saturation",
            args.dnm_oe_noncpg,
            args.dnm_oe_cpg,
            True,
        ),
    ):
        r = np.where(cpg, r_cpg, r_non)
        exp = d.n_sites * (1 - (1 - q) ** r) * k if sat else d.exp * r * k
        out.append(
            {
                "option": label,
                "oe": d.obs.sum() / exp.sum(),
                "oe_noncpg": d.obs[~cpg].sum() / exp[~cpg].sum(),
                "oe_cpg": d.obs[cpg].sum() / exp[cpg].sum(),
            }
        )
    pd.DataFrame(out).to_csv(f"{args.out_prefix}.dnm_factor.tsv", sep="\t", index=False)
    print(pd.DataFrame(out).round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--cells", required=True, help="Decompressed syn_cells TSV.")
    parser.add_argument(
        "--dnm-oe",
        type=float,
        default=2.84,
        help="Cluster DNM obs/exp (outlier_drivers.py).",
    )
    parser.add_argument("--dnm-oe-noncpg", type=float, default=5.24)
    parser.add_argument("--dnm-oe-cpg", type=float, default=1.15)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
