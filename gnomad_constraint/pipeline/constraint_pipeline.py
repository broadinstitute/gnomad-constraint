"""This script builds a constraint pipeline that calculates genic constraint metrics.

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
from gnomad.utils.constraint import build_models
from gnomad.utils.filtering import filter_x_nonpar, filter_y_nonpar
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_constraint.resources.resource_utils import (
    CURRENT_VERSION,
    CUSTOM_VEP_ANNOTATIONS,
    DATA_TYPES,
    GENOMIC_REGIONS,
    MODEL_TYPES,
    POPS,
    VERSIONS,
    annotated_context_ht,
    check_resource_existence,
    constraint_tmp_prefix,
    get_logging_path,
    get_models,
    get_predicted_proportion_observed_dataset,
    get_preprocessed_ht,
    get_sites_resource,
    get_training_dataset,
    mutation_rate_ht,
)
from gnomad_constraint.utils.constraint import (
    add_vep_context_annotations,
    apply_models,
    create_observed_and_possible_ht,
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
    hl.init(
        log="/constraint_pipeline.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    test = args.test
    overwrite = args.overwrite
    max_af = args.max_af
    training_set_partition_hint = args.training_set_partition_hint
    apply_obs_pos_count_partition_hint = args.apply_obs_pos_count_partition_hint
    apply_expected_variant_partition_hint = args.apply_expected_variant_partition_hint
    use_pops = args.use_pops
    use_weights = args.use_weights
    custom_vep_annotation = args.custom_vep_annotation

    # TODO: gnomAD v4 is still in production, for now this will only use 2.1.1.
    version = args.version
    if version not in VERSIONS:
        version = CURRENT_VERSION
        logger.warning(
            "The requested version of resource Tables are not exist, will use gnomAD"
            " v2.1.1 as default."
        )
    # Construct resources with paths for intermediate Tables generated in the pipeline.
    preprocess_resources = {}
    training_resources = {}
    models = {}
    applying_resources = {}

    for region in GENOMIC_REGIONS:
        for data_type in DATA_TYPES:
            if (region == "autosome_par") | (data_type != "genomes"):
                # Save a TableResource with a path to `preprocess_resources`
                preprocess_resources[(region, data_type)] = get_preprocessed_ht(
                    data_type, version, region, test
                )
        # Save a TableResource with a path to `training_resources`
        training_resources[region] = get_training_dataset(version, region, test)

        for model_type in MODEL_TYPES:
            # Save a path to `models`
            models[(region, model_type)] = get_models(model_type, version, region, test)

        # Save a TableResource with a path to `applying_resources`
        applying_resources[region] = get_predicted_proportion_observed_dataset(
            custom_vep_annotation, version, region, test
        )

    try:
        if args.preprocess_data:
            logger.info(
                "Adding VEP context annotations and preparing tables for constraint"
                " calculations..."
            )
            # TODO: Need to add function that annotates methylation, coverage, and
            #  gerp in the vep context table.
            context_ht = annotated_context_ht.versions[version].ht()
            # Raise error if any of the output resources exist and --overwrite is not
            # used.
            check_resource_existence(
                output_step_resources={
                    "--preprocess-data": preprocess_resources.values(),
                },
                overwrite=overwrite,
            )
            # Add annotations used in constraint calculations.
            for data_type in DATA_TYPES:
                if data_type != "context":
                    ht = get_sites_resource(data_type, version).ht()
                else:
                    ht = context_ht

                # Filtering the Table to chr20, chrX, and chrY for testing if
                # applicable.
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

                # Filter input Table and add annotations used in constraint
                # calculations.
                ht = prepare_ht_for_constraint_calculations(ht)
                # Filter to locus that is on an autosome or in a pseudoautosomal region.
                ht.filter(ht.locus.in_autosome_or_par()).write(
                    preprocess_resources[("autosome_par", data_type)].path,
                    overwrite=overwrite,
                )

                # Sex chromosomes are analyzed separately, since they are biologically
                # different from the autosomes.
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

        mutation_ht = mutation_rate_ht.versions[version].ht().select("mu_snp")

        # Create training datasets that include possible and observed variant counts
        # for building models.
        if args.create_training_set:
            logger.info("Counting possible and observed variant counts...")

            # Check if the input/output resources exist.
            check_resource_existence(
                {"--preprocess-data": preprocess_resources.values()},
                {"--create-training-set": training_resources.values()},
                overwrite,
            )
            # Create training datasets for sites on autosomes/pseudoautosomal regions,
            # chromosome X, and chromosome Y.
            for region in GENOMIC_REGIONS:
                create_observed_and_possible_ht(
                    preprocess_resources[(region, "exomes")].ht(),
                    preprocess_resources[(region, "context")].ht(),
                    mutation_ht,
                    max_af=max_af,
                    pops=POPS if use_pops else (),
                    partition_hint=training_set_partition_hint,
                    filter_to_canonical_synonymous=True,
                    global_annotation="training_dataset_params",
                ).write(training_resources[region].path, overwrite=overwrite)
            logger.info("Done with creating training dataset.")

        # Build plateau and coverage models for autosomes/pseudoautosomal regions,
        # chromosome X, and chromosome Y
        if args.build_models:
            # Check if the training datasets exist.
            check_resource_existence(
                {"--create-training-set": training_resources.values()},
                {"--build_models": models.values()},
                overwrite,
            )
            # Build plateau and coverage models.
            for region in GENOMIC_REGIONS:
                logger.info("Building %s plateau and coverage models...", region)

                coverage_model, plateau_models = build_models(
                    training_resources[region].ht(),
                    weighted=use_weights,
                    pops=POPS if use_pops else (),
                )
                hl.experimental.write_expression(
                    plateau_models,
                    models[(region, "plateau")].path,
                    overwrite=overwrite,
                )
                hl.experimental.write_expression(
                    coverage_model,
                    models[(region, "coverage")].path,
                    overwrite=overwrite,
                )
                logger.info("Done building %s plateau and coverage models.", region)

        # Apply coverage and plateau models to compute expected variant counts and
        # observed:expected ratio
        if args.apply_models:
            # Check if the input/output resources exist.
            check_resource_existence(
                {
                    "--preprocess-data": preprocess_resources.values(),
                    "--build_models": models.values(),
                },
                {"--apply_models": applying_resources.values()},
                overwrite,
            )
            # Apply coverage and plateau models for sites on autosomes/pseudoautosomal
            # regions, chromosome X, and chromosome Y.
            for region in GENOMIC_REGIONS:
                logger.info(
                    "Applying %s plateau and coverage models and computing expected"
                    " variant count and observed:expected ratio...",
                    region,
                )
                apply_models(
                    preprocess_resources[(region, "exomes")].ht(),
                    preprocess_resources[(region, "context")].ht(),
                    mutation_ht,
                    models[(region, "plateau")].he(),
                    models[(region, "coverage")].he(),
                    max_af=max_af,
                    pops=POPS if use_pops else (),
                    obs_pos_count_partition_hint=apply_obs_pos_count_partition_hint,
                    expected_variant_partition_hint=apply_expected_variant_partition_hint,
                    custom_vep_annotation=custom_vep_annotation,
                ).write(applying_resources[region].path, overwrite=overwrite)
            logger.info(
                "Done computing expected variant count and observed:expected ratio."
            )

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
        "--use-pops",
        help="Whether to apply models on each population.",
        action="store_true",
    )
    parser.add_argument(
        "--max-af",
        help="Maximum variant allele frequency to keep.",
        type=float,
        default=0.001,
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
    parser.add_argument(
        "--create-training-set",
        help=(
            "Count the observed variants and possible variants by exome coverage at"
            " synonymous sites."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--training-set-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting variants for"
            " training datasets."
        ),
        type=int,
        default=100,
    )
    parser.add_argument(
        "--build-models",
        help="Build plateau and coverage models.",
        action="store_true",
    )
    parser.add_argument(
        "--use-weights",
        help=(
            "Whether to generalize the models to weighted least squares using"
            "'possible_variants'."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-models",
        help=(
            "Apply plateau and coverage models to variants in exome sites Table and"
            " context Table to compute expected variant counts."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-obs-pos-count-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting observed and"
            " expected variants for model application"
        ),
        type=int,
        default=2000,
    )
    parser.add_argument(
        "--apply-expected-variant-partition-hint",
        help=(
            "Target number of partitions for sum aggregators after applying models to"
            " get expected variant counts."
        ),
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--custom-vep-annotation",
        help=(
            "Custom VEP annotation to be used to annotate transcript when"
            ' applying models (one of "transcript_consequences" or "worst_csq_by_gene")'
        ),
        type=str,
        default="transcript_consequences",
        choices=CUSTOM_VEP_ANNOTATIONS,
    )

    args = parser.parse_args()
    main(args)
