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
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.constraint import build_models
from gnomad.utils.filtering import filter_x_nonpar, filter_y_nonpar
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_qc.resource_utils import check_resource_existence  # TODO: fix throughout code
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import (
    CURRENT_VERSION,
    CUSTOM_VEP_ANNOTATIONS,
    DATA_TYPES,
    GENOMIC_REGIONS,
    MODEL_TYPES,
    POPS,
    VERSIONS,
    constraint_tmp_prefix,
    gerp_ht,
    get_annotated_context_ht,
    get_constraint_metrics_dataset,
    get_coverage_ht,
    get_logging_path,
    get_models,
    get_mutation_ht,
    get_predicted_proportion_observed_dataset,
    get_preprocessed_ht,
    get_sites_resource,
    get_training_dataset,
    methylation_ht,
)
from gnomad_constraint.utils.constraint import (
    add_vep_context_annotations,
    apply_models,
    calculate_mu_by_downsampling,
    compute_constraint_metrics,
    create_observed_and_possible_ht,
    prepare_ht_for_constraint_calculations,
    split_context_ht,
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
    prepare_context_ht = args.prepare_context_ht
    training_set_partition_hint = args.training_set_partition_hint
    apply_obs_pos_count_partition_hint = args.apply_obs_pos_count_partition_hint
    apply_expected_variant_partition_hint = args.apply_expected_variant_partition_hint
    use_pops = args.use_pops
    use_weights = args.use_weights
    use_v2_release_mutation_ht = args.use_v2_release_mutation_ht
    custom_vep_annotation = args.custom_vep_annotation

    pops = POPS if use_pops else ()  # TODO: fix throughout code

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

    # For genomes need a preprocessed ht for autosome_par.
    # For exomes need a preprocessed ht for autosome_par, chrX, and chrY.
    # For context need a preprocessed ht for autosome_par, chrX, chrY, and
    # full (autosome_par + chrX + chrY).
    for data_type in DATA_TYPES:
        if data_type != "genomes":
            for region in GENOMIC_REGIONS:
                preprocess_resources[(region, data_type)] = get_preprocessed_ht(
                    data_type, version, region, test
                )
        else:
            preprocess_resources[(region, data_type)] = get_preprocessed_ht(
                data_type, autosome_par, "autosome_par", test
            )

        if data_type == "context":
            preprocess_resources[(region, data_type)] = get_preprocessed_ht(
                data_type, version, "full", test
            )

    for region in GENOMIC_REGIONS:
        if region != "full":
            # Save a TableResource with a path to `training_resources`
            training_resources[region] = get_training_dataset(version, region, test)

            # Save a TableResource with a path to `applying_resources`
            applying_resources[region] = get_predicted_proportion_observed_dataset(
                custom_vep_annotation, version, region, test
            )
            # Save a path to `models`
            for model_type in MODEL_TYPES:
                models[(region, model_type)] = get_models(
                    model_type, version, region, test
                )

    # Save a TableResource with a path to mutation rate Table.
    mutation_rate_resource = get_mutation_ht(version, test, use_v2_release_mutation_ht)
    # Save a TableResource with a path to `constraint_metrics_ht`.
    constraint_metrics_ht = get_constraint_metrics_dataset(version, test)

    try:
        if args.preprocess_data:
            logger.info(
                "Adding VEP context annotations and preparing tables for constraint"
                " calculations..."
            )
            # Annotates methylation, coverage, and gerp in the vep context Table.
            if prepare_context_ht:
                split_context_ht(
                    vep_context.versions["101"].ht(),
                    {
                        "exomes": get_coverage_ht("exomes").ht(),
                        "genomes": get_coverage_ht("genomes").ht(),
                    },
                    methylation_ht.versions[version].ht(),
                    gerp_ht.versions[version].ht(),
                ).write(get_annotated_context_ht(version).path, overwrite)
                context_ht = get_annotated_context_ht(version).ht()
            else:
                context_ht = get_annotated_context_ht(use_old_data=True).ht()
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

                # Checkpoint the full preprocessed context ht containing autsome_par,
                # chrX, and chrY.
                if data_type == "context":
                    ht = ht.checkpoint(
                        preprocess_resources[("full", data_type)].path,
                        overwrite=overwrite,
                    )

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

        # Calculate mutation rate Table.
        if args.calculate_mutation_rate:
            # Save a TableResource with a path to mutation rate Table.
            mutation_rate_resource = get_mutation_ht(version, test)
            logger.info("Calculating mutation rate...")
            # Check if the input/output resources exist.
            check_resource_existence(
                input_step_resources={
                    "--preprocess-data": preprocess_resources.values()
                },
                output_step_resources={
                    "--calculate-mutation-rate": [mutation_rate_resource]
                },
                overwrite=overwrite,
            )

            # Calculate mutation rate using the downsampling with size 1000 genomes in
            # genome site Table.
            calculate_mu_by_downsampling(
                preprocess_resources[("autosome_par", "genomes")].ht(),
                preprocess_resources[
                    ("autosome_par", "context")
                ].ht(),  # TODO: check why this was 'full' before (doesn't exist and immediately filtered to autosomes)
                recalculate_all_possible_summary=True,
                recalculate_all_possible_summary_unfiltered=False,
                pops=pops,
                min_cov=args.min_cov,
                max_cov=args.max_cov,
                gerp_lower_cutoff=args.gerp_lower_cutoff,
                gerp_upper_cutoff=args.gerp_upper_cutoff,
            ).write(mutation_rate_resource.path, overwrite=overwrite)

        # Create training datasets that include possible and observed variant counts
        # for building models.
        if args.create_training_set:
            logger.info("Counting possible and observed variant counts...")

            # Check if the input/output resources exist.
            check_resource_existence(
                input_step_resources={
                    "--preprocess-data": preprocess_resources.values(),
                    "--calculate-mutation-rate": [mutation_rate_resource],
                },
                output_step_resources={
                    "--create-training-set": training_resources.values()
                },
                overwrite=overwrite,
            )
            # Create training datasets for sites on autosomes/pseudoautosomal regions,
            # chromosome X, and chromosome Y.
            for region in GENOMIC_REGIONS:
                create_observed_and_possible_ht(
                    preprocess_resources[(region, "exomes")].ht(),
                    preprocess_resources[(region, "context")].ht(),
                    mutation_rate_resource.ht().select("mu_snp"),
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
                input_step_resources={
                    "--create-training-set": training_resources.values()
                },
                output_step_resources={"--build-models": models.values()},
                overwrite=overwrite,
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
                input_step_resources={
                    "--preprocess-data": preprocess_resources.values(),
                    "--calculate-mutation-rate": [mutation_rate_resource],
                    "--build-models": models.values(),
                },
                output_step_resources={"--apply-models": applying_resources.values()},
                overwrite=overwrite,
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
                    mutation_rate_resource.ht().select("mu_snp"),
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

        if args.compute_constraint_metrics:
            logger.info(
                "Computing constraint metrics, including pLI scores, z scores, oe"
                " ratio, and confidence interval around oe ratio..."
            )

            # Check if the input/output resources exist.
            check_resource_existence(
                input_step_resources={"--apply-models": applying_resources.values()},
                output_step_resources={
                    "--compute-constraint-metrics": [constraint_metrics_ht]
                },
                overwrite=overwrite,
            )
            # Combine Tables of expected variant counts at autosomes/pseudoautosomal
            # regions, chromosome X, and chromosome Y sites.
            hts = [applying_resources[region].ht() for region in GENOMIC_REGIONS]
            union_ht = hts[0].union(*hts[1:])
            union_ht = union_ht.repartition(
                args.compute_constraint_metrics_partitions
            ).checkpoint(new_temp_file(prefix="constraint_apply_union", extension="ht"))

            # Compute constraint metrics
            compute_constraint_metrics(
                union_ht,
                pops=POPS if use_pops else (),
                expected_values={
                    "Null": args.expectation_null,
                    "Rec": args.expectation_rec,
                    "LI": args.expectation_li,
                },
                min_diff_convergence=args.min_diff_convergence,
                raw_z_outlier_threshold=args.raw_z_outlier_threshold,
            ).write(constraint_metrics_ht.path, overwrite=overwrite)
            logger.info("Done with computing constraint metrics.")

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
    preprocess_data_args = parser.add_argument_group("Preprocess data args", "Arguments used for preprocessing the data.")
    preprocess_data_args.add_argument(
        "--prepare-context-ht",
        help=(
            "Whether to split multiallelic sites and add 'methylation', 'coverage', and"
            " 'gerp' annotation to context Table with VEP annotation."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--calculate-mutation-rate",
        help="Calculate baseline mutation rate for each substitution and context using downsampling data.",
        action="store_true",
    )
    mutation_rate_args = parser.add_argument_group("Calcualte mutation rate args", "Arguments used for calculating the muataion rate.")
    mutation_rate_args.add_argument(
        "--min-cov",
        help=(
            "Minimum coverage required to keep a site when calculating the mutation"
            " rate. Default is 15."
        ),
        type=int,
        default=15,
    )
    mutation_rate_args.add_argument(
        "--max-cov",
        help=(
            "Maximum coverage required to keep a site when calculating the mutation"
            " rate. Default is 60."
        ),
        type=int,
        default=60,
    )
    mutation_rate_args.add_argument(
        "--gerp-lower-cutoff",
        help=(
            "Minimum GERP score for variant to be included when calculating the"
            " mutation rate. Default is -3.9885."
        ),
        type=float,
        default=-3.9885,
    )
    mutation_rate_args.add_argument(
        "--gerp-upper-cutoff",
        help=(
            "Maximum GERP score for variant to be included when calculating the"
            " mutation rate. Default is 2.6607."
        ),
        type=float,
        default=2.6607,
    )
    parser.add_argument(
        "--create-training-set",
        help=(
            "Count the observed variants and possible variants by exome coverage at"
            " synonymous sites."
        ),
        action="store_true",
    )
    training_set_args = parser.add_argument_group("Training set args", "Arguments used for creating the training set.")
    training_set_args.add_argument(
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
    build_models_args = parser.add_argument_group("Build models args", "Arguments used for building models.")
    build_models_args.add_argument(
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
    apply_models_args = parser.add_argument_group("Apply models args", "Arguments used for applying the plateau and coverage models.")
    apply_models_args.add_argument(
        "--apply-obs-pos-count-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting observed and"
            " expected variants for model application"
        ),
        type=int,
        default=2000,
    )
    apply_models_args.add_argument(
        "--apply-expected-variant-partition-hint",
        help=(
            "Target number of partitions for sum aggregators after applying models to"
            " get expected variant counts."
        ),
        type=int,
        default=1000,
    )
    apply_models_args.add_argument(
        "--custom-vep-annotation",
        help=(
            "Custom VEP annotation to be used to annotate transcript when"
            ' applying models (one of "transcript_consequences" or "worst_csq_by_gene")'
        ),
        type=str,
        default="transcript_consequences",
        choices=CUSTOM_VEP_ANNOTATIONS,
    )
    compute_constraint_args = parser.add_argument_group("Computate constraint metrics args", "Arguments used for computing constraint metrics.")
    compute_constraint_args.add_argument(
        "--compute-constraint-metrics",
        help=(
            "Compute constraint metrics including pLI scores, z scores, oe ratio, and"
            " confidence interval around oe ratio."
        ),
        action="store_true",
    )
    compute_constraint_args.add_argument(
        "--compute-constraint-metrics-partitions",
        help=(
            "Number of partitions to which the unioned Table of expected variant counts"
            " for autosomes/pseudoautosomal regions, chromosome X, and chromosome Y "
            " should be reaprtitioned."
        ),
        type=int,
        default=1000,
    )

    compute_constraint_args.add_argument(
        "--min-diff-convergence",
        help=(
            "Minimum iteration change in pLI (Probability of loss-of-function"
            " intolerance for haploinsufficient genes) to consider the EM"
            " (expectation-maximization) model convergence criteria as met when"
            " calculating pLI scores."
        ),
        type=float,
        default=0.001,
    )

    compute_constraint_args.add_argument(
        "--expectation-null",
        help=(
            "Expected observed/expected rate of truncating variation for genes where"
            " protein truncating variation is completely tolerated by natural"
            " selection."
        ),
        type=float,
        default=1.0,
    )
    compute_constraint_args.add_argument(
        "--expectation-rec",
        help=(
            "Expected observed/expected rate of truncating variation for recessive"
            " disease genes."
        ),
        type=float,
        default=0.463,
    )
    compute_constraint_args.add_argument(
        "--expectation-li",
        help=(
            "Expected observed/expected rate of truncating variation for severe"
            " haploinsufficient genes."
        ),
        type=float,
        default=0.089,
    )
    compute_constraint_args.add_argument(
        "--raw-z-outlier-threshold",
        help=(
            "Value at which the raw z-score is considered an outlier. Values below the"
            " negative of '--raw-z-outlier-threshold' will be considered outliers for"
            " lof and missense varaint counts (indicating too many variants), whereas"
            " values either above '--raw-z-outlier-threshold' or below the negative of"
            " '--raw-z-outlier-threshold' will be considered outliers for synonymous"
            " varaint counts (indicating too few or too many variants)."
        ),
        type=int,
        default=5,
    )

    args = parser.parse_args()
    main(args)
