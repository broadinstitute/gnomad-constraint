"""Script to build the 1kb regional genomic feature matrix for chrX and chrY windows.

The published Gnocchi feature matrix (genomic_features13_genome_1kb.txt) covers
chr1-22 only. This script computes the same feature columns for chrX (PAR and
non-PAR) and chrY windows from the Gnocchi source tracks (the extracted contents of
gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13.tar.gz). Port of the
chrX PAR feature builder in Ruchit10/gene_deserts (regional_corrections).

Only features whose source data covers chrX/chrY are computed (GC_content, LCR, SINE,
LINE, CpG_island, dist2telo, dist2cent). The other features (``--masked-features``)
are written as NaN (never 0) and are masked when computing rr
(regional_rate_pipeline.py --compute-rr --mask-features).

The script consists of the following steps:
    - --validate: rebuild autosomal windows and compare them to the published matrix.
    - --build-xy-features: build the chrX/chrY feature matrix.
"""

import argparse
import logging
import os
import tempfile

import hail as hl
import pandas as pd

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    GENOMIC_FEATURE_COLUMNS,
    RR_XY_MASKED_FEATURES,
    VERSIONS,
)
from gnomad_constraint.utils.regional_rate import (
    build_feature_matrix,
    load_feature_windows,
    validate_feature_builder,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("build_xy_features")
logger.setLevel(logging.INFO)


def main(args):
    """Build the chrX/chrY regional genomic feature matrix."""
    hl.init(
        log="/build_xy_features.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test = args.test
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_build_xy_features_resources(
        version,
        overwrite,
        test,
        directory_post_fix=args.directory_post_fix,
        path_post_fix=args.path_post_fix,
        output_root=args.output_root,
        tracks_dir=args.tracks_dir,
    )
    tracks_dir = args.tracks_dir or constraint_res.get_genomic_feature_tracks_dir(
        version=version, test=test, output_root=args.output_root
    )

    try:
        if args.validate:
            logger.info("Validating the feature builder on published windows...")
            df = validate_feature_builder(
                tracks_dir,
                constraint_res.get_published_genomic_features_path(),
                args.seq_cache_dir,
                contig=args.validate_contig,
                n_windows=args.validate_n_windows,
            )
            logger.info("Feature builder validation:\n%s", df.to_string(index=False))
            n_failed = int((~df["passed"]).sum())
            if n_failed:
                raise ValueError(
                    f"{n_failed} of {len(df)} feature columns failed validation!"
                )
            logger.info("All %d feature columns passed validation.", len(df))

        if args.build_xy_features:
            logger.info("Building the chrX/chrY regional genomic feature matrix...")
            res = resources.build_xy_features
            res.check_resource_existence()

            frames = []
            for contig in args.contigs:
                element_ids = None
                if test:
                    # Use a slice of windows from the middle of the contig.
                    ids = load_feature_windows(tracks_dir, contig)["1k"]["element_id"]
                    mid = len(ids) // 2
                    element_ids = set(ids.iloc[mid : mid + args.test_n_windows])
                frames.append(
                    build_feature_matrix(
                        tracks_dir,
                        contig,
                        args.seq_cache_dir,
                        args.masked_features,
                        element_ids,
                    )
                )
            df = pd.concat(frames, ignore_index=True)

            populated = [c for c in GENOMIC_FEATURE_COLUMNS if df[c].notna().any()]
            logger.info(
                "%d windows built; %d of %d feature columns populated.",
                len(df),
                len(populated),
                len(GENOMIC_FEATURE_COLUMNS),
            )
            # Write with pandas to keep full float precision (hl.Table.export
            # rounds floats). hl.hadoop_open compresses based on the extension.
            with hl.hadoop_open(res.xy_features_tsv, "w") as f:
                df.to_csv(f, sep="\t", index=False, na_rep="NaN")
            logger.info("Done building the chrX/chrY feature matrix.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path("build_xy_features", version=version)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--version",
        help=(
            "Which version of the resource Tables will be used. Default is"
            f" {CURRENT_VERSION}."
        ),
        type=str,
        default=CURRENT_VERSION,
    )
    parser.add_argument(
        "--directory-post-fix",
        help="Post-fix to append to the output directory path.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--path-post-fix",
        help="Post-fix to append to the output file path.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--output-root",
        help=(
            "Optional bucket path (e.g. gs://gnomad-tmp/rp) to use as the root of the"
            " outputs instead of the production root."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help="Build only a slice of windows of each contig for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-windows",
        help="Number of windows per contig to build with --test. Default is 2000.",
        type=int,
        default=2000,
    )
    parser.add_argument(
        "--tracks-dir",
        help=(
            "Directory (local or gs://) with the extracted genomic_features13.tar.gz"
            " source tracks. Default is the regional_rate genomic_features13_tracks"
            " directory under the output root."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--seq-cache-dir",
        help=(
            "Local directory used to cache hg38 reference sequence fetched from the"
            " UCSC REST API for GC content."
        ),
        type=str,
        default=os.path.join(tempfile.gettempdir(), "ucsc_hg38_sequence_cache"),
    )

    validate_args = parser.add_argument_group(
        "Validate args", "Arguments used for validating the feature builder."
    )
    validate_args.add_argument(
        "--validate",
        help=(
            "Rebuild autosomal windows and compare them to the published feature"
            " matrix. Raises an error if any computed feature column fails."
        ),
        action="store_true",
    )
    validate_args.add_argument(
        "--validate-contig",
        help="Autosome to validate on. Default is chr21.",
        type=str,
        default="chr21",
    )
    validate_args.add_argument(
        "--validate-n-windows",
        help="Number of published windows to rebuild. Default is 3000.",
        type=int,
        default=3000,
    )

    build_args = parser.add_argument_group(
        "Build chrX/chrY features args",
        "Arguments used for building the chrX/chrY feature matrix.",
    )
    build_args.add_argument(
        "--build-xy-features",
        help="Build the chrX/chrY regional genomic feature matrix.",
        action="store_true",
    )
    build_args.add_argument(
        "--contigs",
        help="Contigs to build. Default is chrX and chrY.",
        nargs="+",
        default=["chrX", "chrY"],
    )
    build_args.add_argument(
        "--masked-features",
        help=(
            "Features without source data covering chrX/chrY, written as NaN."
            f" Default is {RR_XY_MASKED_FEATURES}."
        ),
        nargs="+",
        default=RR_XY_MASKED_FEATURES,
    )

    args = parser.parse_args()
    main(args)
