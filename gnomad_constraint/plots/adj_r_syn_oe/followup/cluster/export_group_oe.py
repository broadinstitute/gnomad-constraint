"""Export per-transcript observed and expected counts from downstream constraint groups.

For each rr_downstream version given, reads its constraint group Table and writes one
row per transcript with the adj observed, expected and unadjusted expected counts of
each constraint group (syn, mis, lof, ...), to ``<out>/<version>.group_oe.tsv.bgz``.
These are what the report scores (synonymous noise, outliers, missense in
unconstrained genes) without needing the metrics step.
"""

import argparse

import hail as hl

ROOT = "gs://gnomad-tmp/jg/rr_downstream"
PATH = "gnomad_v4.1.1/constraint/apply_models/transcript_consequences/gnomad.v4.1.1.constraint_group.ht"


def main(args: argparse.Namespace) -> None:
    """Export each version's table."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    for v in args.versions:
        ht = hl.read_table(f"{ROOT}/{v}/{PATH}")
        if args.print_schema:
            print(v, "schema:", ht.row.dtype)
            print(v, "globals:", ht.globals.dtype)
            print(v, "row:", ht.head(1).collect())
            continue
        # One element per constraint group; take the adj element of each oe_info.
        meta = hl.eval(ht.constraint_group_meta)
        adj = hl.eval(ht.exomes_freq_meta).index({"group": "adj"})
        names = ["_".join(f"{k}_{val}" for k, val in m.items()) for m in meta]
        print(v, "groups:", names)
        cg = ht.constraint_groups
        ht = ht.select(
            *[
                f
                for f in ("gene", "transcript", "canonical", "mane_select")
                if f not in ht.key
            ],
            **{
                f"{n}_{f}": cg[i].oe_info[adj][f]
                for i, n in enumerate(names)
                for f in ("observed_variants", "expected_variants", "unadj_expected")
            },
            **{f"{n}_possible": cg[i].possible_variants for i, n in enumerate(names)},
        )
        ht.export(f"{args.out}/{v}.group_oe.tsv.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--versions", nargs="+", required=True)
    parser.add_argument("--out", default="gs://gnomad-tmp/jg/rr_downstream/group_oe")
    parser.add_argument("--print-schema", action="store_true")
    main(parser.parse_args())
