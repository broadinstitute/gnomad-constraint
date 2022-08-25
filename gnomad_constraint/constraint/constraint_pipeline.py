# noqa: D100
# cSpell: disable
import argparse
import logging

import hail as hl

from gnomad.resources.grch37.gnomad import public_release
from gnomad_constraint.resources.resource_utils import (
    preprocessed_ht,
    get_logging_path,
    annotated_context_ht,
    mutation_rate_ht_path,
    training_ht_path,
)
from gnomad_constraint.utils.constraint_basics import (
    add_vep_context_annotations,
    prepare_ht_for_constraint_calculations,
    get_proportion_observed_by_coverage,
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

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("constraint_pipeline"))

    context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
    genome_ht = full_genome_ht.filter(full_genome_ht.locus.in_autosome_or_par())
    exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
    mutation_ht = hl.read_table(mutation_rate_ht_path).select("mu_snp")

    context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval("X")])
    context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
    context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval("Y")])
    context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())

    exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval("X")])
    exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
    exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval("Y")])
    exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

    if args.create_training_set:
        get_proportion_observed_by_coverage(
            exome_ht,
            context_ht,
            mutation_ht,
            True,
            args.dataset,
            not args.skip_af_filter_upfront,
        ).write(training_ht_path, overwrite=args.overwrite)
        hl.read_table(training_ht_path).export(
            training_ht_path.replace(".ht", ".txt.bgz")
        )
        get_proportion_observed_by_coverage(
            exome_x_ht,
            context_x_ht,
            mutation_ht,
            True,
            args.dataset,
            not args.skip_af_filter_upfront,
        ).write(training_ht_path.replace(".ht", "_x.ht"), overwrite=args.overwrite)
        get_proportion_observed_by_coverage(
            exome_y_ht,
            context_y_ht,
            mutation_ht,
            True,
            args.dataset,
            not args.skip_af_filter_upfront,
        ).write(training_ht_path.replace(".ht", "_y.ht"), overwrite=args.overwrite)
        logger.info("Done with creating training dataset.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files", action="store_true"
    )
    parser.add_argument(
        "--skip_af_filter_upfront",
        help="Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended",
        action="store_true",
    )
    parser.add_argument(
        "--dataset",
        help="Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)",
        default="gnomad",
    )
    parser.add_argument(
        "--preprocess-data",
        help="Whether to add necessary coverage, methylation level, and VEP annotations to genome and exome Tables.",
        action="store_true",
    )
    parser.add_argument(
        "--create-training-set",
        help="Count the observed variants and possible variants by exome coverage at synonymous site.",
        action="store_true",
    )
