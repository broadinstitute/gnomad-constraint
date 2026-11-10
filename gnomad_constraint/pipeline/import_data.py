"""
Script to import the external data used to compute the regional mutation rate (rr).

The imported Tables are written under ``imported_data/`` in the constraint output
directory.

The import script consists of the following steps:
    - --import-dnms: import the published de novo mutations (DNMs) as a Table of DNM
      sites.
    - --import-genomic-features: import the Gnocchi 1kb regional genomic feature
      matrix.
"""

import argparse
import logging

import hail as hl

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    GENOMIC_FEATURE_COLUMNS,
    VERSIONS,
)
from gnomad_constraint.utils.import_data import format_dnm_ht, format_genomic_features

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("import_data")
logger.setLevel(logging.INFO)


def main(args: argparse.Namespace) -> None:
    """Import the external data used to compute the regional mutation rate."""
    hl.init(
        log="/import_data.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_import_data_resources(version, overwrite, args.test)

    try:
        if args.import_dnms:
            logger.info("Importing the DNMs...")
            res = resources.import_dnms
            res.check_resource_existence()

            # The files are gzip (not bgzip) compressed and quote fields with commas.
            # TODO: The calls aren't QC-filtered. The ASD file has Confidence and Role
            #  columns and the DDD file has LCR and gnomad_non_neuro_AF columns that a
            #  filter could use.
            ht = hl.Table.union(
                *[
                    hl.import_table(path, delimiter=",", quote='"', force=True).select(
                        "Variant"
                    )
                    for path in [res.asd_dnm_path, res.ddd_dnm_path]
                ]
            )
            format_dnm_ht(ht).write(res.dnm_ht.path, overwrite=overwrite)
            logger.info("Done importing the DNMs.")

        if args.import_genomic_features:
            logger.info("Importing the regional genomic feature matrix...")
            res = resources.import_genomic_features
            res.check_resource_existence()

            ht = hl.import_table(
                res.published_features_path,
                types={c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
                missing=["NA", "nan", "NaN", ""],
                min_partitions=100,
            )
            ht = format_genomic_features(ht)
            ht = ht.checkpoint(res.genomic_features_ht.path, overwrite=overwrite)
            n_duplicated = ht.count() - ht.distinct().count()
            if n_duplicated:
                raise ValueError(f"{n_duplicated} windows are duplicated!")
            logger.info("Done importing the regional genomic feature matrix.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(constraint_res.get_logging_path("import_data", version=version))


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
        "--test",
        help="Whether to write the outputs to the testing bucket.",
        action="store_true",
    )
    parser.add_argument(
        "--import-dnms",
        help="Import the published DNMs as a Table of DNM sites.",
        action="store_true",
    )
    parser.add_argument(
        "--import-genomic-features",
        help="Import the Gnocchi 1kb regional genomic feature matrix.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
