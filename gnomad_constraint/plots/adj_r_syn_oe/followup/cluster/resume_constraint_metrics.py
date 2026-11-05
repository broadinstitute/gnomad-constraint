"""Finish --compute-constraint-metrics from a checkpoint a failed run left in temp.

``compute_constraint_metrics`` checkpoints each step to a randomly named temp Table
(``constraint_metrics.oe_ci_z-*``, then ``.z_scores-*``, ``.percentile_bins-*`` and
``.pli-*``), so a rerun starts over, including the 2-3 hour confidence-interval step.
This picks up from such a checkpoint instead:

- ``--map``: print each checkpoint's synonymous expected count for its first few
  transcripts, to match it to the downstream version that wrote it (whose
  per-transcript counts give the same values).
- ``--resume``: run the remaining steps of ``compute_constraint_metrics`` from an
  ``oe_ci_z`` or ``z_scores`` checkpoint, with the pipeline's default settings, and
  write the constraint metrics Table of ``--output-root`` (never overwriting).
  ``constraint_pipeline.py --prepare-release --export-release-tsv`` then runs as usual.
"""

import argparse

import hail as hl

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.utils.constraint import (
    _compute_percentile_bins,
    _compute_pli_scores,
    _compute_z_scores,
    add_gencode_transcript_annotations,
)

VERSION = "4.1.1"
EXPECTED_VALUES = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}
"""The pipeline's defaults for --expectation-null, --expectation-rec, --expectation-li."""


def main(args: argparse.Namespace) -> None:
    """Map or resume."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    if args.map:
        for path in args.checkpoints:
            ht = hl.read_table(path)
            adj = hl.eval(ht.exomes_freq_meta).index({"group": "adj"})
            # The first few rows only: a full pass over ~10,000 partitions is slow.
            rows = (
                ht.select(
                    syn_exp=ht.constraint_groups[0].oe_info[adj].expected_variants
                )
                .head(5)
                .collect()
            )
            for r in rows:
                print(f"MAP {path} {r.transcript} {r.syn_exp}")
        return

    out = constraint_res.get_constraint_metrics_dataset(
        version=VERSION, output_root=args.output_root
    ).path
    if hl.hadoop_exists(f"{out}/_SUCCESS"):
        raise ValueError(f"{out} already exists!")
    # The checkpoints have ~10,000 partitions for ~310,000 rows; read fewer.
    ht = hl.read_table(args.checkpoint, _n_partitions=args.n_partitions)
    if args.stage == "oe_ci_z":
        ht = _compute_z_scores(ht)
        ht = ht.checkpoint(hl.utils.new_temp_file("constraint_metrics.z_scores", "ht"))
    ht = _compute_percentile_bins(ht, True)
    ht = ht.checkpoint(
        hl.utils.new_temp_file("constraint_metrics.percentile_bins", "ht")
    )
    ht = _compute_pli_scores(ht, EXPECTED_VALUES, 0.001)
    ht = ht.checkpoint(hl.utils.new_temp_file("constraint_metrics.pli", "ht"))
    gq = constraint_res.get_gene_quality_metrics_ht(VERSION).ht()
    ht = ht.annotate(**gq[ht.transcript])
    ht = add_gencode_transcript_annotations(ht, constraint_res.get_gencode_ht(VERSION))
    ht.write(out)
    print(f"Wrote {out} from {args.checkpoint}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", action="store_true")
    parser.add_argument("--checkpoints", nargs="+", default=[])
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--checkpoint")
    parser.add_argument("--stage", choices=["oe_ci_z", "z_scores"])
    parser.add_argument("--output-root")
    parser.add_argument("--n-partitions", type=int, default=200)
    main(parser.parse_args())
