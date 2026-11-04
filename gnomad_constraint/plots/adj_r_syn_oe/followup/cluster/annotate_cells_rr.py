"""Add rr from the full-run rr Tables to Ruchit's exported synonymous cells.

Each cell is (transcript, 1kb window, context). A window's interval in every rr Table
starts at the window's BED start, so the locus at BED start + 1 is inside it. Cells
without a window get a missing rr (no correction).
"""

import argparse

import hail as hl

ROOTS = {
    "full": "gs://gnomad-tmp/jg/rr_full",
    "qcexcl": "gs://gnomad-tmp/jg/rr_full_qcexcl",
}
CUTOFFS = {"bonferroni": "bonferroni", "1e-2": "1e-2", "ridge": "ridge"}


def main(args: argparse.Namespace) -> None:
    """Annotate the cells and export them."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    ht = hl.import_table(
        args.cells,
        types={"exp": hl.tfloat64, "obs": hl.tint64, "n_sites": hl.tint64},
        force_bgz=True,
        min_partitions=64,
    )
    win = ht.win.split("-")
    ht = ht.annotate(
        _locus=hl.locus(win[0], hl.int(win[1]) + 1, reference_genome="GRCh38")
    )
    for name, path in args.rr_tables:
        rr_ht = hl.read_table(path).select("adj_r")
        ht = ht.annotate(**{f"rr_{name}": rr_ht[ht._locus].adj_r.get(ht.context)})
    ht.drop("_locus").export(args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True)
    parser.add_argument(
        "--rr-tables", nargs=2, action="append", metavar=("NAME", "PATH"), required=True
    )
    parser.add_argument("--out", required=True)
    main(parser.parse_args())
