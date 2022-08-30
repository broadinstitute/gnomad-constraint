# noqa: D100

import argparse
import logging

import hail as hl

from gnomad.resources.grch38.reference_data import vep_context
from gnomad.resources.grch37.gnomad import public_release
from gnomad_constraint.utils.constraint_basics import (
    add_vep_context_annotations,
    prepare_ht_for_constraint_calculations,
)

from gnomad_constraint.resources.resource_utils import get_processed_ht_path

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def main(args):
    """Execute the constraint pipeline."""
    trimers = args.trimers
    if args.pre_process_data:
        add_vep_context_annotations(
            public_release("genomes").ht(), vep_context.versions["101"].path
        ).write(get_processed_ht_path("genomes"), overwrite=args.overwrite)
        add_vep_context_annotations(
            public_release("exomes").ht(), vep_context.versions["101"].path
        ).write(get_processed_ht_path("exomes"), overwrite=args.overwrite)
        logger.info("Done with preprocessing genome and exome Table.")

    full_context_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(vep_context.versions["101"].path), trimers=trimers
    )
    full_genome_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(get_processed_ht_path("genomes")), trimers=trimers
    )
    full_exome_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(get_processed_ht_path("exomes")), trimers=trimers
    )


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
