# noqa: D100

import argparse
import logging

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad.utils.filtering import (
    filter_x_nonpar,
    filter_y_nonpar,
)
from gnomad_constraint.resources.resource_utils import (
    get_preprocessed_ht,
    get_logging_path,
    annotated_context_ht,
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
    try:
        if args.preprocess_data:
            logger.info("Adding VEP context annotations...")
            # Add annotations from VEP context Table to gnomAD data.
            # TODO: Need to add function that annotates methylation, coverage, and gerp in the vep context table.
            preprocessed_genome_ht = add_vep_context_annotations(
                public_release("genomes").ht(), annotated_context_ht.ht()
            )
            preprocessed_exome_ht = add_vep_context_annotations(
                public_release("exomes").ht(), annotated_context_ht.ht()
            )
            # Filter input Table and add annotations used in constraint calculations.
            full_context_ht = prepare_ht_for_constraint_calculations(
                annotated_context_ht.ht()
            )
            full_genome_ht = prepare_ht_for_constraint_calculations(
                preprocessed_genome_ht
            )
            full_exome_ht = prepare_ht_for_constraint_calculations(
                preprocessed_exome_ht
            )
            # Filter to locus that is on an autosome or a pseudoautosomal region in context Table, exome Table, and genome Table.
            context_ht = full_context_ht.filter(
                full_context_ht.locus.in_autosome_or_par()
            ).write(get_preprocessed_ht("context").path, overwrite=args.overwrite)
            genome_ht = full_genome_ht.filter(
                full_genome_ht.locus.in_autosome_or_par()
            ).write(get_preprocessed_ht("genome").path, overwrite=args.overwrite)
            exome_ht = full_exome_ht.filter(
                full_exome_ht.locus.in_autosome_or_par()
            ).write(get_preprocessed_ht("exome").path, overwrite=args.overwrite)
            # Sex chromosomes are analyzed separately, since they are biologically different from the autosomes.
            context_x_ht = filter_x_nonpar(full_context_ht).write(
                get_preprocessed_ht("context", "chrx").path, overwrite=args.overwrite
            )
            context_y_ht = filter_y_nonpar(full_context_ht).write(
                get_preprocessed_ht("context", "chry").path, overwrite=args.overwrite
            )
            exome_x_ht = filter_x_nonpar(full_exome_ht).write(
                get_preprocessed_ht("exome", "chrx").path, overwrite=args.overwrite
            )
            exome_y_ht = filter_y_nonpar(full_exome_ht).write(
                get_preprocessed_ht("exome", "chrx").path, overwrite=args.overwrite
            )
            logger.info("Done with preprocessing genome and exome Table.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("constraint_pipeline"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files", action="store_true"
    )
    parser.add_argument(
        "--preprocess-data",
        help="Whether to prepare the exome, genome, and context Table for constraint calculations by adding necessary coverage, methylation level, and VEP annotations.",
        action="store_true",
    )
