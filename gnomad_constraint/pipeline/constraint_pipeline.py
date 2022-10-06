# noqa: D100
# cSpell: disable
import argparse
import logging

import hail as hl

from gnomad.utils.filtering import (
    filter_x_nonpar,
    filter_y_nonpar,
)
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_constraint.resources.resource_utils import (
    get_sites_resource,
    get_preprocessed_ht,
    get_logging_path,
    annotated_context_ht,
    CURRENT_VERSION,
    VERSIONS,
    DATA_TYPES,
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
            "The requested version of resource Tables are not exist, will use v2 as default."
        )

    try:
        if args.preprocess_data:
            logger.info("Adding VEP context annotations...")
            # TODO: Need to add function that annotates methylation, coverage, and gerp in the vep context table.
            context_ht = annotated_context_ht.versions[version].ht()
            # Add annotations used in constraint calculations.
            for data_type in DATA_TYPES:
                if data_type != "context":
                    ht = get_sites_resource(data_type, version).ht()
                else:
                    ht = context_ht

                # Filtering the Table to chr20, chrX, and chrY for testing if applicable
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
                preprocessed_ht = (
                    add_vep_context_annotations(ht, context_ht)
                    if data_type != "context"
                    else ht
                )

                # Filter input Table and add annotations used in constraint calculations.
                full_ht = prepare_ht_for_constraint_calculations(preprocessed_ht)
                # Filter to locus that is on an autosome or in a pseudoautosomal region.
                ht = full_ht.filter(full_ht.locus.in_autosome_or_par())
                ht.write(
                    get_preprocessed_ht(data_type, version, "autosome_par", test).path,
                    overwrite=overwrite,
                )

                # Sex chromosomes are analyzed separately, since they are biologically different from the autosomes.
                if data_type != "genomes":
                    filter_x_nonpar(full_ht).write(
                        get_preprocessed_ht(
                            data_type, version, "chrx_nonpar", test
                        ).path,
                        overwrite=overwrite,
                    )
                    filter_y_nonpar(full_ht).write(
                        get_preprocessed_ht(
                            data_type, version, "chry_nonpar", test
                        ).path,
                        overwrite=overwrite,
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
        help=f"Which version of the resource Tables will be used. Default is {CURRENT_VERSION}.",
        type=str,
        default=CURRENT_VERSION,
    )
    parser.add_argument(
        "--test",
        help="Whether to filter the exome Table, genome Table and the context Table to only chromosome 20, chromosome X, and chromosome Y for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--preprocess-data",
        help="Whether to prepare the exome, genome, and context Table for constraint calculations by adding necessary coverage, methylation level, and VEP annotations.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
