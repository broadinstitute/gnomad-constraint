"""Count which synonymous outliers each version of the full run still leaves.

The outliers are the groups in ``outlier_drivers.all_genes.tsv`` (the 47 non-histone
genes with |z| > 8 and no QC flag under the ridge rr, split by explanation), plus the
histone genes. For each version, z comes from ``full_run_scores.py``'s per-transcript
tables (``full_run`` and ``rr_rates``), as in the outlier counts elsewhere. Writes
``<prefix>.tsv``: per group and version, the genes in the group and how many still
have |z| > 8, and the remaining genes.
"""

import argparse

import pandas as pd
from outlier_drivers import HISTONE

GROUPS = [
    ("maternal hotspot region", "maternal hotspot region"),
    ("repeat, segdup or low coverage", "repeat, segdup or low coverage"),
    ("other", "other"),
]
VERSIONS = ["none", "ridge_qcexcl", "roulette_ridge", "homebuilt_ridge"]


def main(args: argparse.Namespace) -> None:
    """Count the remaining outliers per group and version."""
    genes = pd.read_csv(args.genes, sep="\t")
    pt = pd.read_csv(args.per_transcript[0], sep="\t")
    for path in args.per_transcript[1:]:
        more = pd.read_csv(path, sep="\t")
        pt = pt.merge(
            more[[c for c in more if c not in pt or c == "transcript"]], on="transcript"
        )
    d = genes.merge(pt, on="transcript")
    histone = d.gene.fillna("").str.match(HISTONE)
    masks = {label: d.outlier & d.group.eq(g) for label, g in GROUPS}
    masks["histone genes"] = histone
    rows = []
    for v in VERSIONS:
        hit = d[f"z_{v}"].abs() > 8
        for label, m in masks.items():
            rows.append(
                {
                    "group": label,
                    "version": v,
                    "n_genes": int(m.sum()),
                    "still_outliers": int((m & hit).sum()),
                    "genes": ",".join(sorted(d.gene[m & hit].astype(str))),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    print(
        out.pivot(index="group", columns="version", values="still_outliers")[VERSIONS]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", default="outlier_drivers.all_genes.tsv")
    parser.add_argument(
        "--per-transcript",
        nargs="+",
        default=["full_run.per_transcript.tsv.gz", "rr_rates.per_transcript.tsv.gz"],
    )
    parser.add_argument("--out-prefix", default="outlier_rates")
    main(parser.parse_args())
