"""Run ``full_run_scores.py`` on the Dataproc driver.

Copies the inputs from GCS to the driver's disk, decompresses the annotated cells, scores
every version and copies the ``<prefix>.*`` tables to ``--out``. ``--rr-names`` is
passed through (for example to score the background runs). Submit with a zip of
``full_run_scores.py``, ``remap_rr_to_shet.py`` and ``rr_shrinkage_experiments.py``.
"""

import argparse
import glob
import os
import subprocess
import sys
from types import SimpleNamespace

import full_run_scores

EXPORT = "gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export"


def main(args: argparse.Namespace) -> None:
    """Copy the inputs, score and copy the results back."""
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "openpyxl"])
    work = "/tmp/full_run_scores"
    os.makedirs(work, exist_ok=True)
    src = {
        "cells.tsv.bgz": args.cells,
        "windows.tsv.bgz": f"{EXPORT}/windows.tsv.bgz",
        "transcripts.tsv.bgz": f"{EXPORT}/transcripts.tsv.bgz",
        "release.tsv.bgz": f"{args.inputs}/gnomad.v4.1.2.constraint_metrics.tsv.bgz",
        "shet.xlsx": f"{args.inputs}/supp_table3.xlsx",
    }
    if not args.rr_names:
        src["cv_full.tsv"] = f"{args.inputs}/rr_cv_metrics.full.tsv"
        src["cv_qcexcl.tsv"] = f"{args.inputs}/rr_cv_metrics.qcexcl.tsv"
    for name, path in src.items():
        subprocess.check_call(["gsutil", "-q", "cp", path, f"{work}/{name}"])
    subprocess.check_call(
        f"gunzip -c {work}/cells.tsv.bgz > {work}/cells.tsv", shell=True
    )
    full_run_scores.main(
        SimpleNamespace(
            cells=f"{work}/cells.tsv",
            windows=f"{work}/windows.tsv.bgz",
            transcripts=f"{work}/transcripts.tsv.bgz",
            release_tsv=f"{work}/release.tsv.bgz",
            shet=f"{work}/shet.xlsx",
            cv_metrics_full=f"{work}/cv_full.tsv",
            cv_metrics_qcexcl=f"{work}/cv_qcexcl.tsv",
            out_prefix=f"{work}/{args.prefix}",
            rr_names=args.rr_names,
        )
    )
    outs = glob.glob(f"{work}/{args.prefix}.*")
    subprocess.check_call(["gsutil", "-q", "cp", *outs, f"{args.out}/"])
    print("Copied", outs, "to", args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells",
        default="gs://gnomad-tmp/jg/rr_downstream/syn_cells_rr_full_run.tsv.bgz",
    )
    parser.add_argument(
        "--inputs", default="gs://gnomad-tmp/jg/rr_downstream/full_run_scores/inputs"
    )
    parser.add_argument(
        "--out", default="gs://gnomad-tmp/jg/rr_downstream/full_run_scores/out"
    )
    parser.add_argument("--prefix", default="full_run")
    parser.add_argument("--rr-names", nargs="+")
    main(parser.parse_args())
