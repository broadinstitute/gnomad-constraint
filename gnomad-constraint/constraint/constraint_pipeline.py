# noqa: D100

import argparse
import logging

import hail as hl

from gnomad.resources.grch38.reference_data import vep_context
from gnomad_qc.v2.resources.basics import get_gnomad_public_data
from utils.constraint_basics import (
    add_vep_context_annotations,
    prepare_ht_for_constraint_calculations,
)

from resources.resource_utils import processed_genomes_ht_path, processed_exomes_ht_path

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
            get_gnomad_public_data("genomes"), vep_context.versions["101"].path
        ).write(processed_genomes_ht_path, overwrite=args.overwrite)
        add_vep_context_annotations(
            get_gnomad_public_data("exomes"), vep_context.versions["101"].path
        ).write(processed_exomes_ht_path, overwrite=args.overwrite)
        logger.info("Done with preprocessing genome and exome Table.")

    full_context_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(vep_context.versions["101"].path), trimers=trimers
    )
    full_genome_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(processed_genomes_ht_path), trimers=trimers
    )
    full_exome_ht = prepare_ht_for_constraint_calculations(
        hl.read_table(processed_exomes_ht_path), trimers=trimers
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
