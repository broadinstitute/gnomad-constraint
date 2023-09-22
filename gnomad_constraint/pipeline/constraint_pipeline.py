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
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.utils.constraint import (
    add_vep_context_annotations,
    annotate_context_ht,
    apply_models,
    calculate_gerp_cutoffs,
    calculate_mu_by_downsampling,
    compute_constraint_metrics,
    create_observed_and_possible_ht,
    prepare_ht_for_constraint_calculations,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def filter_for_test(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Filter `ht` to chr20, chrX, and chrY for testing.

    :param ht: Table to filter.
    :param data_type: Data type of `ht`.
    :return: Filtered Table for testing.
    """
    rg = get_reference_genome(ht.locus)
    contigs_keep = [
        hl.parse_locus_interval(c, reference_genome=rg)
        for c in [rg.contigs[19], rg.x_contigs[0], rg.y_contigs[0]]
    ]
    logger.info(
        "Filtering the %s HT to chr20, chrX, and chrY for testing...",
        data_type,
    )
    ht = hl.filter_intervals(ht, contigs_keep)

    return ht


def get_constraint_resources(
    version: str,
    use_v2_release_mutation_ht: bool,
    use_v2_release_context_ht: bool,
    custom_vep_annotation: str,
    overwrite: bool,
    test: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the constraint pipeline.

    :param version: Version of constraint resources to use.
    :param use_v2_release_mutation_ht: Whether to use the v2 release mutation ht.
    :param use_v2_release_context_ht: Whether to use the v2 release context ht.
    :param custom_vep_annotation: Custom VEP annotation to use for applying models
        resources.
    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :return: PipelineResourceCollection containing resources for all steps of the
        constraint pipeline.
    """
    data_types = constraint_res.DATA_TYPES
    regions = constraint_res.GENOMIC_REGIONS
    # Initialize constraint pipeline resource collection.
    constraint_pipeline = PipelineResourceCollection(
        pipeline_name="constraint",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the constraint pipeline.
    context_res = constraint_res.vep_context_ht.versions[version]
    context_build = get_reference_genome(context_res.ht().locus).name
    prepare_context = PipelineStepResourceCollection(
        "--prepare-context-ht",
        output_resources={
            "annotated_context_ht": constraint_res.get_annotated_context_ht(
                version, use_v2_release_context_ht
            )
        },
        input_resources={
            "gnomAD resources": {
                "context_ht": context_res,
                "exomes_coverage_ht": constraint_res.get_coverage_ht("exomes"),
                "genomes_coverage_ht": constraint_res.get_coverage_ht("genomes"),
                "methylation_ht": constraint_res.get_methylation_ht(context_build),
            },
        },
    )
    # For genomes need a preprocessed ht for autosome_par.
    # For exomes and context need a preprocessed ht for autosome_par, chrX,
    # and chrY.
    preprocess_data = PipelineStepResourceCollection(
        "--preprocess-data",
        output_resources={
            f"preprocessed_{r}_{d}_ht": constraint_res.get_preprocessed_ht(
                d, version, r, test
            )
            for r in regions
            for d in data_types
            if (r == "autosome_par") | (d != "genomes")
        },
        pipeline_input_steps=[prepare_context],
        add_input_resources={
            "gnomAD sites resources": {
                f"{d}_sites_ht": constraint_res.get_sites_resource(d, version)
                for d in data_types
                if d != "context"
            }
        },
    )
    calculate_gerp_cutoffs = PipelineStepResourceCollection(
        "--calculate-gerp-cutoffs",
        output_resources={},
        pipeline_input_steps=[preprocess_data],
    )
    calculate_mutation_rate = PipelineStepResourceCollection(
        "--calculate-mutation-rate",
        output_resources={
            "mutation_ht": constraint_res.get_mutation_ht(
                version, test, use_v2_release_mutation_ht
            )
        },
        pipeline_input_steps=[preprocess_data],
    )
    create_training_set = PipelineStepResourceCollection(
        "--create-training-set",
        output_resources={
            f"train_{r}_ht": constraint_res.get_training_dataset(version, r, test)
            for r in regions
        },
        pipeline_input_steps=[preprocess_data, calculate_mutation_rate],
    )
    build_models = PipelineStepResourceCollection(
        "--build-models",
        output_resources={
            f"model_{r}_{m}_ht": constraint_res.get_models(m, version, r, test)
            for m in constraint_res.MODEL_TYPES
            for r in regions
        },
        pipeline_input_steps=[create_training_set],
    )
    apply_models = PipelineStepResourceCollection(
        "--apply-models",
        output_resources={
            f"apply_{r}_ht": constraint_res.get_predicted_proportion_observed_dataset(
                custom_vep_annotation, version, r, test
            )
            for r in regions
        },
        pipeline_input_steps=[preprocess_data, calculate_mutation_rate, build_models],
    )
    compute_constraint_metrics = PipelineStepResourceCollection(
        "--compute-constraint-metrics",
        output_resources={
            "constraint_metrics_ht": constraint_res.get_constraint_metrics_dataset(
                version, test
            )
        },
        pipeline_input_steps=[apply_models],
    )

    # Add all steps to the constraint pipeline resource collection.
    constraint_pipeline.add_steps(
        {
            "prepare_context": prepare_context,
            "preprocess_data": preprocess_data,
            "calculate_gerp_cutoffs": calculate_gerp_cutoffs,
            "calculate_mutation_rate": calculate_mutation_rate,
            "create_training_set": create_training_set,
            "build_models": build_models,
            "apply_models": apply_models,
            "compute_constraint_metrics": compute_constraint_metrics,
        }
    )

    return constraint_pipeline


def main(args):
    """Execute the constraint pipeline."""
    hl.init(
        log="/constraint_pipeline.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    regions = constraint_res.GENOMIC_REGIONS
    version = args.version
    test = args.test
    overwrite = args.overwrite

    max_af = args.max_af
    use_pops = args.use_pops
    use_v2_release_mutation_ht = args.use_v2_release_mutation_ht
    custom_vep_annotation = args.custom_vep_annotation
    gerp_lower_cutoff = args.gerp_lower_cutoff
    gerp_upper_cutoff = args.gerp_upper_cutoff

    pops = constraint_res.POPS if use_pops else ()

    if version not in constraint_res.VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = get_constraint_resources(
        version,
        use_v2_release_mutation_ht,
        args.use_v2_release_context_ht,
        custom_vep_annotation,
        overwrite,
        test,
    )

    try:
        if args.prepare_context_ht:
            logger.info(
                "Annotating methylation, coverage, and GERP on the VEP context Table..."
            )
            res = resources.prepare_context
            res.check_resource_existence()
            context_ht = res.context_ht.ht()
            coverage_hts = {
                "exomes": res.exomes_coverage_ht.ht(),
                "genomes": res.genomes_coverage_ht.ht(),
            }
            annotate_context_ht(
                context_ht,
                coverage_hts,
                res.methylation_ht.ht(),
                constraint_res.get_gerp_ht(get_reference_genome(context_ht.locus).name),
            ).write(res.annotated_context_ht.path, overwrite)

        if args.preprocess_data:
            logger.info(
                "Adding VEP context annotations and preparing tables for constraint"
                " calculations..."
            )
            res = resources.preprocess_data
            res.check_resource_existence()
            context_ht = res.annotated_context_ht.ht()

            # Add annotations used in constraint calculations.
            for data_type in constraint_res.DATA_TYPES:
                if data_type != "context":
                    ht = getattr(res, f"{data_type}_sites_ht").ht()
                else:
                    ht = context_ht

                if test:
                    ht = filter_for_test(ht, data_type)

                # Add annotations from VEP context Table to genome and exome Tables.
                if data_type != "context":
                    ht = add_vep_context_annotations(ht, context_ht)

                # Filter input Table and add annotations used in constraint
                # calculations.
                ht = prepare_ht_for_constraint_calculations(
                    ht, require_exome_coverage=(data_type == "exomes")
                )
                # Filter to locus that is on an autosome or in a pseudoautosomal region.
                ht.filter(ht.locus.in_autosome_or_par()).write(
                    getattr(res, f"preprocessed_autosome_par_{data_type}_ht").path,
                    overwrite=overwrite,
                )
                # Sex chromosomes are analyzed separately, since they are biologically
                # different from the autosomes.
                if data_type != "genomes":
                    filter_x_nonpar(ht).write(
                        getattr(res, f"preprocessed_chrx_nonpar_{data_type}_ht").path,
                        overwrite=overwrite,
                    )
                    filter_y_nonpar(ht).write(
                        getattr(res, f"preprocessed_chry_nonpar_{data_type}_ht").path,
                        overwrite=overwrite,
                    )
            logger.info("Done with preprocessing genome and exome Table.")

        if args.calculate_gerp_cutoffs:
            logger.warning(
                "Calculating new GERP cutoffs to be used instead of"
                " '--gerp-lower-cutoff' and '--gerp-upper-cutoff' defaults."
            )
            res = resources.calculate_gerp_cutoffs
            res.check_resource_existence()
            gerp_lower_cutoff, gerp_upper_cutoff = calculate_gerp_cutoffs(
                res.preprocessed_autosome_par_context_ht.ht()
            )
            logger.info(
                "Calculated new GERP cutoffs: using a lower GERP cutoff of %f and an"
                " upper GERP cutoff of %f.",
                gerp_lower_cutoff,
                gerp_upper_cutoff,
            )

        if args.calculate_mutation_rate:
            logger.info("Calculating mutation rate...")
            res = resources.calculate_mutation_rate
            res.check_resource_existence()

            # Calculate mutation rate using the downsampling with size 1000 genomes in
            # genome site Table.
            calculate_mu_by_downsampling(
                res.preprocessed_autosome_par_genomes_ht.ht(),
                res.preprocessed_autosome_par_context_ht.ht(),
                recalculate_all_possible_summary=True,
                pops=pops,
                min_cov=args.min_cov,
                max_cov=args.max_cov,
                gerp_lower_cutoff=gerp_lower_cutoff,
                gerp_upper_cutoff=gerp_upper_cutoff,
            ).repartition(args.mutation_rate_partitions).write(
                res.mutation_ht.path, overwrite=overwrite
            )

        # Create training datasets that include possible and observed variant counts
        # for building models.
        if args.create_training_set:
            logger.info("Counting possible and observed variant counts...")
            res = resources.create_training_set
            res.check_resource_existence()

            # Create training datasets for sites on autosomes/pseudoautosomal regions,
            # chromosome X, and chromosome Y.
            for r in regions:
                op_ht = create_observed_and_possible_ht(
                    getattr(res, f"preprocessed_{r}_exomes_ht").ht(),
                    getattr(res, f"preprocessed_{r}_context_ht").ht(),
                    res.mutation_ht.ht().select("mu_snp"),
                    max_af=max_af,
                    pops=pops,
                    partition_hint=args.training_set_partition_hint,
                    transcript_for_synonymous_filter=(
                        "mane_select" if int(version[0]) >= 4 else "canonical"
                    ),  # Switch to using MANE Select transcripts rather than canonical for gnomAD v4 and later versions.
                    global_annotation="training_dataset_params",
                )
                if use_v2_release_mutation_ht:
                    op_ht = op_ht.annotate_globals(use_v2_release_mutation_ht=True)
                op_ht.write(getattr(res, f"train_{r}_ht").path, overwrite=overwrite)
            logger.info("Done with creating training dataset.")

        if args.build_models:
            res = resources.build_models
            res.check_resource_existence()

            # Build plateau and coverage models for autosomes/pseudoautosomal regions,
            # chromosome X, and chromosome Y
            for r in regions:
                logger.info("Building %s plateau and coverage models...", r)
                coverage_model, plateau_models = build_models(
                    getattr(res, f"train_{r}_ht").ht(),
                    weighted=args.use_weights,
                    pops=pops,
                )
                hl.experimental.write_expression(
                    plateau_models,
                    getattr(res, f"model_{r}_plateau").path,
                    overwrite=overwrite,
                )
                hl.experimental.write_expression(
                    coverage_model,
                    getattr(res, f"model_{r}_coverage").path,
                    overwrite=overwrite,
                )
                logger.info("Done building %s plateau and coverage models.", r)

        if args.apply_models:
            res = resources.apply_models
            res.check_resource_existence()

            # Apply separate plateau models for sites on autosomes/pseudoautosomal
            # regions, chromosome X, and chromosome Y. Use autosomes/pseudoautosomal
            # coverage models for all contigs (Note: should test separate coverage models
            # for XX/XY in the future).
            for r in regions:
                logger.info(
                    "Applying %s plateau and autosome coverage models and computing"
                    " expected variant count and observed:expected ratio...",
                    r,
                )
                oe_ht = apply_models(
                    getattr(res, f"preprocessed_{r}_exomes_ht").ht(),
                    getattr(res, f"preprocessed_{r}_context_ht").ht(),
                    # mutation_ht,
                    res.mutation_ht.ht().select("mu_snp"),
                    getattr(res, f"model_{r}_plateau").he(),
                    getattr(res, "model_autosome_par_coverage").he(),
                    max_af=max_af,
                    pops=pops,
                    obs_pos_count_partition_hint=args.apply_obs_pos_count_partition_hint,
                    expected_variant_partition_hint=args.apply_expected_variant_partition_hint,
                    custom_vep_annotation=custom_vep_annotation,
                    use_mane_select_instead_of_canonical=(
                        True if int(version[0]) >= 4 else False
                    ),  # Group by MANE Select transcripts rather than canonical for gnomAD v4 and later versions.
                )
                if use_v2_release_mutation_ht:
                    oe_ht = oe_ht.annotate_globals(use_v2_release_mutation_ht=True)
                oe_ht.write(getattr(res, f"apply_{r}_ht").path, overwrite=overwrite)
            logger.info(
                "Done computing expected variant count and observed:expected ratio."
            )

        if args.compute_constraint_metrics:
            logger.info(
                "Computing constraint metrics, including pLI scores, z scores, oe"
                " ratio, and confidence interval around oe ratio..."
            )
            res = resources.compute_constraint_metrics
            res.check_resource_existence()

            # Combine Tables of expected variant counts at autosomes/pseudoautosomal
            # regions, chromosome X, and chromosome Y sites.
            hts = [getattr(res, f"apply_{r}_ht").ht() for r in regions]
            union_ht = hts[0].union(*hts[1:])
            union_ht = union_ht.repartition(args.compute_constraint_metrics_partitions)
            union_ht = union_ht.checkpoint(
                new_temp_file(prefix="constraint_apply_union", extension="ht")
            )

            # Compute constraint metrics
            compute_constraint_metrics(
                union_ht,
                pops=pops,
                keys=tuple(
                    [
                        i
                        for i in list(union_ht.key)
                        if i in ["gene", "transcript", "canonical", "mane_select"]
                    ]
                ),
                expected_values={
                    "Null": args.expectation_null,
                    "Rec": args.expectation_rec,
                    "LI": args.expectation_li,
                },
                min_diff_convergence=args.min_diff_convergence,
                raw_z_outlier_threshold=args.raw_z_outlier_threshold,
            ).write(res.constraint_metrics_ht.path, overwrite=overwrite)
            logger.info("Done with computing constraint metrics.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(constraint_res.get_logging_path("constraint_pipeline", version))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--version",
        help=(
            "Which version of the resource Tables will be used. Default is"
            f" {constraint_res.CURRENT_VERSION}."
        ),
        type=str,
        default=constraint_res.CURRENT_VERSION,
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to filter the exome Table, genome Table and the context Table to"
            " only chromosome 20, chromosome X, and chromosome Y for testing."
        ),
        action="store_true",
    )

    prepare_context_args = parser.add_argument_group(
        "Prepare context Table args", "Arguments used for preparing the context Table."
    )
    prepare_context_args.add_argument(
        "--prepare-context-ht",
        help=(
            "Prepare the context Table by splitting multiallelic sites and adding "
            "'methylation', 'coverage', and 'gerp' annotations to context Table with "
            "VEP annotation."
        ),
        action="store_true",
    )

    preprocess_data_args = parser.add_argument_group(
        "Preprocess data args", "Arguments used for preprocessing the data."
    )
    preprocess_data_args.add_argument(
        "--use-v2-release-context-ht",
        help="Whether to use the annotated context Table for the v2 release.",
        action="store_true",
    )
    preprocess_data_args.add_argument(
        "--preprocess-data",
        help=(
            "Whether to prepare the exome, genome, and context Table for constraint"
            " calculations by adding necessary coverage, methylation level, and VEP"
            " annotations."
        ),
        action="store_true",
    )

    parser.add_argument(
        "--calculate-gerp-cutoffs",
        help=(
            "Calculate GERP lower and upper cutoffs based on 5th and 95th percentiles"
            " of GERP scores in the context Table. Note that if"
            " '--calculate-gerp-cutoffs' is specified, the computed values will be used"
            " in all downstream steps of the constraint pipeline instead of the values"
            " defined by --gerp-lower-cutoff and --gerp-upper-cutoff."
        ),
        action="store_true",
    )

    mutation_rate_args = parser.add_argument_group(
        "Calculate mutation rate args",
        "Arguments used for calculating the muataion rate.",
    )

    recalculate_mutation_rate = mutation_rate_args.add_argument(
        "--calculate-mutation-rate",
        help=(
            "Calculate baseline mutation rate for each substitution and context using"
            " downsampling data."
        ),
        action="store_true",
    )

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
            " mutation rate. Default is -3.9885 (precalculated on the GRCh37 context"
            " Table to define the 95th percentile)."
        ),
        type=float,
        default=-3.9885,
    )
    mutation_rate_args.add_argument(
        "--gerp-upper-cutoff",
        help=(
            "Maximum GERP score for variant to be included when calculating the"
            " mutation rate. Default is 2.6607 (precalculated on the GRCh37 context"
            " Table to define the 5th percentile)."
        ),
        type=float,
        default=2.6607,
    )
    mutation_rate_args.add_argument(
        "--mutation-rate-partitions",
        help=(
            "Number of partitions to which the mutation rate Table should be"
            " repartitioned."
        ),
        type=int,
        default=1,
    )

    training_set_args = parser.add_argument_group(
        "Training set args", "Arguments used for creating the training set."
    )

    training_set_args.add_argument(
        "--create-training-set",
        help=(
            "Count the observed variants and possible variants by exome coverage at"
            " synonymous sites."
        ),
        action="store_true",
    )

    training_set_args.add_argument(
        "--training-set-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting variants for"
            " training datasets."
        ),
        type=int,
        default=100,
    )

    # `max-af` is an arg for both `--create-training-set` and `--apply-models`
    maximum_af = training_set_args.add_argument(
        "--max-af",
        help="Maximum variant allele frequency to keep.",
        type=float,
        default=0.001,
    )

    # `use_populations` is an arg for `--create-training-set`, `--apply-models`, `--build-models`, and `compute_constraint_args`
    use_populations = training_set_args.add_argument(
        "--use-pops",
        help=(
            "Whether to apply models on each population. If not specified, will use"
            " 'global'."
        ),
        action="store_true",
    )

    use_v2_release_mutation_rate = training_set_args.add_argument(
        "--use-v2-release-mutation-ht",
        help="Whether to use the mutatation rate computed for the v2 release.",
        action="store_true",
    )

    mutation_rate_parser = parser.add_mutually_exclusive_group(required=False)
    mutation_rate_parser._group_actions.append(use_v2_release_mutation_rate)
    mutation_rate_parser._group_actions.append(recalculate_mutation_rate)

    build_models_args = parser.add_argument_group(
        "Build models args", "Arguments used for building models."
    )

    build_models_args.add_argument(
        "--build-models",
        help="Build plateau and coverage models.",
        action="store_true",
    )

    build_models_args.add_argument(
        "--use-weights",
        help=(
            "Whether to generalize the models to weighted least squares using"
            "'possible_variants'."
        ),
        action="store_true",
    )

    build_models_args._group_actions.append(use_populations)

    apply_models_args = parser.add_argument_group(
        "Apply models args",
        "Arguments used for applying the plateau and coverage models.",
    )

    apply_models_args.add_argument(
        "--apply-models",
        help=(
            "Apply plateau and coverage models to variants in exome sites Table and"
            " context Table to compute expected variant counts."
        ),
        action="store_true",
    )

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
        choices=constraint_res.CUSTOM_VEP_ANNOTATIONS,
    )
    apply_models_args._group_actions.append(maximum_af)
    apply_models_args._group_actions.append(use_populations)
    apply_models_args._group_actions.append(use_v2_release_mutation_rate)

    compute_constraint_args = parser.add_argument_group(
        "Computate constraint metrics args",
        "Arguments used for computing constraint metrics.",
    )
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

    compute_constraint_args._group_actions.append(use_populations)

    args = parser.parse_args()
    main(args)
