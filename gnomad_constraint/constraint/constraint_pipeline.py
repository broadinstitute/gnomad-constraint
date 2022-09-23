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
    preprocessed_ht,
    get_logging_path,
    training_dataset,
    annotated_context_ht,
    mutation_rate_ht,
    POPS,
)
from gnomad_constraint.utils.constraint_basics import (
    add_vep_context_annotations,
    prepare_ht_for_constraint_calculations,
    create_constraint_training_dataset,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def main(args):
    """Execute the constraint pipeline."""
    max_af = args.max_af
    partitions = args.partitions
    use_pops = args.use_pop

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
            # Filter to locus that is on an autosome or a pseudoautosomal region in context Table, exome Table, and genome Table
            full_context_ht.filter(full_context_ht.locus.in_autosome_or_par()).write(
                preprocessed_ht("context").path, overwrite=args.overwrite
            )
            full_genome_ht.filter(full_genome_ht.locus.in_autosome_or_par()).write(
                preprocessed_ht("genome").path, overwrite=args.overwrite
            )
            full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par()).write(
                preprocessed_ht("exome").path, overwrite=args.overwrite
            )
            logger.info("Done with preprocessing genome and exome Table.")

        if partitions:
            context_ht = (
                preprocessed_ht("context").ht()._filter_partitions(range(partitions))
            )
            exome_ht = (
                preprocessed_ht("genome").ht()._filter_partitions(range(partitions))
            )
        else:
            context_ht = preprocessed_ht("context").ht()
            exome_ht = preprocessed_ht("genome").ht()

        mutation_ht = mutation_rate_ht.ht().select("mu_snp")

        context_x_ht = filter_x_nonpar(full_context_ht)
        context_y_ht = filter_y_nonpar(full_context_ht)

        exome_x_ht = filter_x_nonpar(full_exome_ht)
        exome_y_ht = filter_y_nonpar(full_exome_ht)

        if args.create_training_set:
            create_constraint_training_dataset(
                exome_ht,
                context_ht,
                mutation_ht,
                max_af=max_af,
                pops=POPS if use_pops else None,
            ).write(training_dataset().path, overwrite=args.overwrite)
            create_constraint_training_dataset(
                exome_x_ht,
                context_x_ht,
                mutation_ht,
                max_af=max_af,
                pops=POPS if use_pops else None,
            ).write(training_dataset("chrx").path, overwrite=args.overwrite)
            create_constraint_training_dataset(
                exome_y_ht,
                context_y_ht,
                mutation_ht,
                max_af=max_af,
                pops=POPS if use_pops else None,
            ).write(training_dataset("chry").path, overwrite=args.overwrite)
            logger.info("Done with creating training dataset.")

        training_ht = training_dataset().ht()
        training_x_ht = training_dataset("chrx").ht()
        training_y_ht = training_dataset("chry").ht()

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("constraint_pipeline"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files", action="store_true"
    )
    parser.add_argument(
        "--use-pop",
        help="Whether to apply models on each population.",
        action="store_true",
    )
    parser.add_argument(
        "--max-af",
        help="Maximum variant allele frequency to keep.",
        nargs="?",
        const=0.001,
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--partitions",
        help="Whether to filters the exome Table and the context Table to only a few partitions.",
        nargs="?",
        const=10,
        type=int,
        default=10,
    )
    parser.add_argument(
        "--preprocess-data",
        help="Whether to add necessary coverage, methylation level, and VEP annotations to genome and exome Tables.",
        action="store_true",
    )
    parser.add_argument(
        "--create-training-set",
        help="Count the observed variants and possible variants by exome coverage at synonymous sites.",
        action="store_true",
    )
