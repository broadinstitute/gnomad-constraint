"""
Script to import external data used with the constraint pipeline.

The imported Tables are inputs to ``regional_rate_pipeline.py``, and are written under
``imported_data/`` in the constraint output directory.

The import script consists of the following steps:
    - --import-dnms: combine the de novo mutation (DNM) cohorts into one Table of SNV
      DNMs annotated with the annotated context Table's context fields and VEP most
      severe consequence.
    - --import-genomic-features: import the 1kb regional genomic feature matrix.
"""

import argparse
import logging
from typing import Union

import hail as hl
from gnomad.resources.resource_utils import TableResource

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    GENOMIC_FEATURE_COLUMNS,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import (
    filter_for_test,
    get_dnm_cohort_resources,
)
from gnomad_constraint.utils.import_data import (
    combine_dnms,
    format_dnm_cohort,
    format_genomic_features,
    summarize_genomic_features,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("import_data")
logger.setLevel(logging.INFO)


def import_dnm_cohort(cohort: str, resource: Union[str, TableResource]) -> hl.Table:
    """
    Import one DNM cohort (see ``get_dnm_cohort_resources``).

    :param cohort: Cohort name.
    :param resource: Path (text cohorts) or TableResource (Hail Table cohorts).
    :return: Imported cohort Table.
    """
    if isinstance(resource, TableResource):
        return resource.ht()
    if cohort.startswith("satterstrom_2026"):
        # The files are gzip (not bgzip) compressed and quote fields with commas.
        return hl.import_table(resource, delimiter=",", quote='"', force=True)

    return hl.import_table(resource, impute=True)


def main(args: argparse.Namespace) -> None:
    """Import external data used with the constraint pipeline."""
    hl.init(
        log="/import_data.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_import_data_resources(version, overwrite, test)

    try:
        if args.import_dnms:
            logger.info("Combining the DNM cohorts into the DNM Table...")
            res = resources.import_dnms
            res.check_resource_existence()

            cohort_hts = {
                cohort: format_dnm_cohort(cohort, import_dnm_cohort(cohort, resource))
                for cohort, resource in get_dnm_cohort_resources().items()
            }
            # The cohorts are small, so counting them is cheap. Comparing these
            # counts with the SNV DNMs kept per cohort shows a cohort that lost many
            # rows in the context join (e.g. positions on the wrong build).
            n_input = {cohort: ht.count() for cohort, ht in cohort_hts.items()}
            context_ht = res.annotated_context_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht, use_gene_list=test_gene_list)

            ht = combine_dnms(cohort_hts, context_ht)
            ht = ht.checkpoint(res.dnm_ht.path, overwrite=overwrite)

            # The DNM Table is small and checkpointed, so aggregating is cheap.
            stats = ht.aggregate(
                hl.struct(
                    n_sites=hl.agg.count(),
                    n_per_csq=hl.agg.counter(ht.most_severe_consequence),
                    n_per_cohort=hl.agg.explode(
                        lambda kv: hl.agg.group_by(kv[0], hl.agg.sum(kv[1])),
                        ht.cohort_counts.items(),
                    ),
                )
            )
            for cohort, n in n_input.items():
                logger.info(
                    "    %s: %d DNMs imported (after liftover), %d SNV DNMs kept.",
                    cohort,
                    n,
                    stats.n_per_cohort.get(cohort, 0),
                )
            logger.info("%d DNM sites after deduplication.", stats.n_sites)
            logger.info(
                "DNM sites per most severe consequence: %s", dict(stats.n_per_csq)
            )
            logger.info("Done creating the DNM Table.")

        if args.import_genomic_features:
            logger.info("Importing the regional genomic feature matrix...")
            res = resources.import_genomic_features
            res.check_resource_existence()

            path = res.published_features
            ht = hl.import_table(
                path,
                types={c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
                missing=["NA", "nan", "NaN", ""],
                min_partitions=100,
                force=path.endswith(".gz"),
            )
            ht = format_genomic_features(ht)
            ht = ht.checkpoint(res.genomic_features_ht.path, overwrite=overwrite)

            summary = summarize_genomic_features(ht)
            if summary["n_duplicated_windows"]:
                raise ValueError(
                    f"{summary['n_duplicated_windows']} windows are duplicated in the"
                    " imported feature matrix!"
                )
            for contig, s in sorted(summary["contigs"].items()):
                logger.info(
                    "    %s: %d windows; features missing in all windows: %s",
                    contig,
                    s["n_windows"],
                    s["missing_features"],
                )
            logger.info("Done importing the regional genomic feature matrix.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(constraint_res.get_logging_path("import_data", version=version))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
        help=(
            "Whether to filter the annotated context Table to only chromosome 20,"
            " chromosome X, and chromosome Y for testing, and write the outputs to"
            " the testing bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-gene-list",
        help=(
            "Whether to filter the annotated context Table to only six genes for"
            " testing (see --test-gene-list in constraint_pipeline.py)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--import-dnms",
        help=(
            "Combine the DNM cohorts into one Table of SNV DNMs annotated with their"
            " context and VEP most severe consequence."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--import-genomic-features",
        help="Import the 1kb regional genomic feature matrix.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
