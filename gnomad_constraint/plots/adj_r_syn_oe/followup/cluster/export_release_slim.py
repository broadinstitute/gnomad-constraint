"""Export a slim copy of each downstream version's release TSV.

The release TSVs are about 550 MB each. This keeps the MANE Select and canonical rows
and the columns the report scores (obs, exp, obs/exp, z-score and upper interval of
syn, mis and LoF, with and without rr; LOEUF rank and decile; pLI; flags), writing
``<out>/<version>.tsv.bgz``. Versions whose slim copy exists are skipped.
"""

import argparse

import hail as hl

ROOT = "gs://gnomad-tmp/jg/rr_downstream"
PATH = "gnomad_v4.1.1/constraint/release/gnomad.v4.1.1.constraint_metrics.tsv.bgz"
ID_COLS = [
    "gene",
    "gene_id",
    "transcript",
    "canonical",
    "mane_select",
    "chromosome",
    "gene_flags",
    "constraint_flags",
]
METRICS = ["obs", "exp", "oe", "z_score", "oe_ci.upper"]
UNADJ = ["unadj.exp", "unadj.oe", "unadj.oe_ci.upper", "unadj.z_score"]
LOF_ONLY = ["oe_ci.upper_rank", "oe_ci.upper_bin_decile", "pLI"]


def main(args: argparse.Namespace) -> None:
    """Export each version's slim TSV."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    cols = ID_COLS + [
        f"{csq}.{m}"
        for csq in ("syn", "mis", "lof")
        for m in METRICS + UNADJ + (LOF_ONLY if csq == "lof" else [])
    ]
    for v in args.versions:
        out = f"{args.out}/{v}.tsv.bgz"
        if hl.hadoop_exists(out):
            print(f"{out} exists, skipping")
            continue
        ht = hl.import_table(f"{ROOT}/{v}/{PATH}", force_bgz=True)
        ht = ht.filter((ht.mane_select == "true") | (ht.canonical == "true"))
        ht.select(**{c: ht[c] for c in cols}).export(out)
        print(f"Wrote {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--versions", nargs="+", required=True)
    parser.add_argument("--out", default=f"{ROOT}/release_slim")
    main(parser.parse_args())
