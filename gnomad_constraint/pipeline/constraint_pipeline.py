"""
This script builds a constraint pipeline that calculates constraint metrics.

The constraint metrics pipeline will compute for LoF variants, missense variants, and
synonymous variants:
    - The number of observed variants
    - The number of expected variants
    - The observed: expected ratio
    - The confidence interval around the observed: expected ratio
    - pLI score (Probability of loss-of-function intolerance; probability that
      transcript falls into distribution of haploinsufficient genes)
    - pNull (Probability that transcript falls into distribution of unconstrained genes)
    - pRec (Probability that transcript falls into distribution of recessive genes)
    - z-scores (Measure how constrained or intolerant a gene or transcript is to
      missense variants and synonymous variants)

The constraint pipeline consists of the following parts:
    - preprocess data
    - create training set
    - build models
    - apply models
    - calculate metrics
"""

import argparse
import logging

import hail as hl
from gnomad.utils.filtering import filter_x_nonpar, filter_y_nonpar
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_constraint.resources.resource_utils import (
    CURRENT_VERSION,
    DATA_TYPES,
    VERSIONS,
    annotated_context_ht,
    get_logging_path,
    get_preprocessed_ht,
    get_sites_resource,
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
    hl.init(log="/constraint_pipeline.log")

    test = args.test
    overwrite = args.overwrite
    # TODO: gnomAD v4 is still in production, for now this will only use 2.1.1.
    version = args.version
    if version not in VERSIONS:
        version = CURRENT_VERSION
        logger.warning(
            "The requested version of resource Tables are not exist, will use gnomAD"
            " v2.1.1 as default."
        )
    preprocess_resources = {}
    training_resources = {}
    for region in GENOMIC_REGIONS:
        for data_type in DATA_TYPES:
            if (region == "autosome_par") | (data_type != "genomes"):
                preprocess_resources[(region, data_type)] = get_preprocessed_ht(
                    data_type, version, region, test
                )
        training_resources[region] = get_training_dataset(version, region, test)
    try:
        if args.preprocess_data:
            logger.info("Adding VEP context annotations...")
            # TODO: Need to add function that annotates methylation, coverage, and
            #  gerp in the vep context table.
            context_ht = annotated_context_ht.versions[version].ht()
            # Raise error if any of the output resources exist and --overwrite is not used.
            check_resource_existence(
                output_pipeline_step="--preprocess-data",
                output_resources=preprocess_resources.values(),
                overwrite=overwrite,
            )
            # Add annotations used in constraint calculations.
            for data_type in DATA_TYPES:
                if data_type != "context":
                    ht = get_sites_resource(data_type, version).ht()
                else:
                    ht = context_ht

                # Filtering the Table to chr20, chrX, and chrY for testing if applicable.
                if test:
                    rg = get_reference_genome(context_ht.locus)
                    contigs_keep = [
                        hl.parse_locus_interval(c, reference_genome=rg)
                        for c in [rg.contigs[19], rg.x_contigs[0], rg.y_contigs[0]]
                    ]
                    logger.info(
                        "Filtering the %s HT to chr20, chrX, and chrY for testing...",
                        data_type,
                    )
                    ht = hl.filter_intervals(ht, contigs_keep)

                # Add annotations from VEP context Table to genome and exome Tables.
                if data_type != "context":
                    ht = add_vep_context_annotations(ht, context_ht)

                # Filter input Table and add annotations used in constraint calculations.
                ht = prepare_ht_for_constraint_calculations(ht)
                # Filter to locus that is on an autosome or in a pseudoautosomal region.
                ht.filter(ht.locus.in_autosome_or_par()).write(
                    preprocess_resources[("autosome_par", data_type)].path,
                    overwrite=overwrite,
                )

                # Sex chromosomes are analyzed separately, since they are biologically
                # different from the autosomes.```
                if data_type != "genomes":
                    filter_x_nonpar(ht).write(
                        preprocess_resources[("chrx_nonpar", data_type)].path,
                        overwrite=overwrite,
                    )
                    filter_y_nonpar(ht).write(
                        preprocess_resources[("chry_nonpar", data_type)].path,
                        overwrite=overwrite,
                    )
            logger.info("Done with preprocessing genome and exome Table.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("constraint_pipeline"))


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
        help=(
            "Whether to filter the exome Table, genome Table and the context Table to"
            " only chromosome 20, chromosome X, and chromosome Y for testing."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--preprocess-data",
        help=(
            "Whether to prepare the exome, genome, and context Table for constraint"
            " calculations by adding necessary coverage, methylation level, and VEP"
            " annotations."
        ),
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
