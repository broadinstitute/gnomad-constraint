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
    CURRENT_VERSION,
    VERSIONS,
)
from gnomad_constraint.utils.constraint import (
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
    # The v4 Tables will be updated in the future. We will use v2 as for now.
    version = str(float(args.version))
    if version not in VERSIONS:
        version = CURRENT_VERSION
        logger.warning(
            "The requested version of resource Tables are not exist, will use v2 as default."
        )

    try:
        if args.preprocess_data:
            logger.info("Adding VEP context annotations...")
            # Add annotations from VEP context Table to gnomAD data.
            # TODO: Need to add function that annotates methylation, coverage, and gerp in the vep context table.
            preprocessed_genome_ht = add_vep_context_annotations(
                public_release("genomes").ht(),
                annotated_context_ht.versions[version].ht(),
            )
            preprocessed_exome_ht = add_vep_context_annotations(
                public_release("exomes").ht(),
                annotated_context_ht.versions[version].ht(),
            )
            # Filter input Table and add annotations used in constraint calculations.
            full_context_ht = prepare_ht_for_constraint_calculations(
                annotated_context_ht.versions[version].ht()
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
            ).write(
                get_preprocessed_ht("context", version).path, overwrite=args.overwrite
            )
            genome_ht = full_genome_ht.filter(
                full_genome_ht.locus.in_autosome_or_par()
            ).write(
                get_preprocessed_ht("genome", version).path, overwrite=args.overwrite
            )
            exome_ht = full_exome_ht.filter(
                full_exome_ht.locus.in_autosome_or_par()
            ).write(
                get_preprocessed_ht("exome", version).path, overwrite=args.overwrite
            )
            # Sex chromosomes are analyzed separately, since they are biologically different from the autosomes.
            context_x_ht = filter_x_nonpar(full_context_ht).write(
                get_preprocessed_ht("context", "chrx", version).path,
                overwrite=args.overwrite,
            )
            context_y_ht = filter_y_nonpar(full_context_ht).write(
                get_preprocessed_ht("context", "chry", version).path,
                overwrite=args.overwrite,
            )
            exome_x_ht = filter_x_nonpar(full_exome_ht).write(
                get_preprocessed_ht("exome", "chrx", version).path,
                overwrite=args.overwrite,
            )
            exome_y_ht = filter_y_nonpar(full_exome_ht).write(
                get_preprocessed_ht("exome", "chrx", version).path,
                overwrite=args.overwrite,
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
        "--version",
        help="Which version of the resource Tables will be used. Default is v2.",
        nargs="?",
        const=2,
        type=int,
        default=2,
    )
    parser.add_argument(
        "--preprocess-data",
        help="Whether to prepare the exome, genome, and context Table for constraint calculations by adding necessary coverage, methylation level, and VEP annotations.",
        action="store_true",
    )
