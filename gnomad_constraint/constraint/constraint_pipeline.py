# noqa: D100

import argparse
import logging

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad.resources.grch38.reference_data import vep_context
from gnomad_constraint.resources.resource_utils import (
    get_processed_ht_path,
    get_logging_path,
    context_ht_path,
)
from gnomad_constraint.utils.constraint_basics import (
    add_vep_context_annotations,
    prepare_ht_for_constraint_calculations,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def main(args):
    """Execute the constraint pipeline."""
    trimers = args.trimers

    try:
        if args.pre_process_data:
            logger.info("Adding VEP context annotations...")
            add_vep_context_annotations(
                public_release("genomes").ht(), context_ht_path
            ).write(get_processed_ht_path("genomes"), overwrite=args.overwrite)
            add_vep_context_annotations(
                public_release("exomes").ht(), context_ht_path
            ).write(get_processed_ht_path("exomes"), overwrite=args.overwrite)
            logger.info("Done with preprocessing genome and exome Table.")

        full_context_ht = prepare_ht_for_constraint_calculations(
            hl.read_table(context_ht_path), trimers=trimers
        )
        full_genome_ht = prepare_ht_for_constraint_calculations(
            hl.read_table(get_processed_ht_path("genomes")), trimers=trimers
        )
        full_exome_ht = prepare_ht_for_constraint_calculations(
            hl.read_table(get_processed_ht_path("exomes")), trimers=trimers
        )
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("constraint_pipeline"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files", action="store_true"
    )
    parser.add_argument(
        "--trimers",
        help="Whether to use trimers instead of heptamers to define the context",
        action="store_true",
    )
    parser.add_argument(
        "--pre-process-data",
        help="Whether to add annotations from VEP context Table to genome and exome Table.",
        action="store_true",
    )
