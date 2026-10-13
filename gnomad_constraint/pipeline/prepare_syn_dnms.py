"""Script to create the synonymous de novo mutation (DNM) Table used to train rr models.

The synonymous DNM Table is the training input of the regional mutation rate (rr)
models (see regional_rate_pipeline.py). DNMs from multiple cohorts are lifted over to
GRCh38 (when needed), combined, filtered to synonymous DNMs using the VEP most severe
consequence of the constraint pipeline's annotated context Table, and deduplicated by
locus and alleles.

The script consists of the following step:
    - --prepare-syn-dnms: combine the DNM cohorts into the synonymous DNM Table.
"""

import argparse
import logging

import hail as hl

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import CURRENT_VERSION, VERSIONS
from gnomad_constraint.resources.resource_utils import (
    filter_for_test,
    get_dnm_cohort_resources,
)
from gnomad_constraint.utils.regional_rate import combine_syn_dnms, load_dnm_cohort

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("prepare_syn_dnms")
logger.setLevel(logging.INFO)


def main(args):
    """Create the synonymous DNM Table."""
    hl.init(
        log="/prepare_syn_dnms.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_prepare_syn_dnms_resources(
        version,
        overwrite,
        test,
        directory_post_fix=args.directory_post_fix,
        path_post_fix=args.path_post_fix,
        output_root=args.output_root,
    )

    try:
        if args.prepare_syn_dnms:
            logger.info("Combining DNM cohorts into the synonymous DNM Table...")
            res = resources.prepare_syn_dnms
            res.check_resource_existence()

            cohort_hts = {
                cohort: load_dnm_cohort(cohort, resource)
                for cohort, resource in get_dnm_cohort_resources().items()
            }
            context_ht = res.annotated_context_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht, use_gene_list=test_gene_list)

            ht, cohort_counts = combine_syn_dnms(cohort_hts, context_ht)
            ht = ht.checkpoint(res.syn_dnm_ht.path, overwrite=overwrite)

            # The synonymous DNM Table is small and checkpointed, so counting is
            # cheap.
            n_dnms = ht.count()
            n_before_dedup = sum(cohort_counts.values())
            logger.info("Synonymous DNMs per cohort (before deduplication):")
            for cohort, n in sorted(cohort_counts.items()):
                logger.info("    %s: %d", cohort, n)
            logger.info(
                "%d synonymous DNMs after deduplication (%d duplicates removed).",
                n_dnms,
                n_before_dedup - n_dnms,
            )
            logger.info(
                "Synonymous DNMs per contig: %s",
                dict(ht.aggregate(hl.agg.counter(ht.locus.contig))),
            )

            ht.export(res.syn_dnm_tsv)
            logger.info("Done creating the synonymous DNM Table.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path("prepare_syn_dnms", version=version)
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
        help=(
            "Whether to filter the context Table to only chromosome 20, chromosome X,"
            " and chromosome Y for testing."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-gene-list",
        help=(
            "Whether to filter the context Table to only a list of genes on"
            " chromosome 20, chromosome X, and chromosome Y for testing."
        ),
        action="store_true",
    )

    prepare_syn_dnms_args = parser.add_argument_group(
        "Prepare synonymous DNMs args",
        "Arguments used for creating the synonymous DNM Table.",
    )
    prepare_syn_dnms_args.add_argument(
        "--prepare-syn-dnms",
        help=(
            "Combine the DNM cohorts, filter to synonymous DNMs, and deduplicate by"
            " locus and alleles."
        ),
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
