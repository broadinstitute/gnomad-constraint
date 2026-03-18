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
from gnomad.utils.constraint import (
    build_models,
    calculate_gerp_cutoffs,
    explode_downsamplings_oe,
)
from gnomad.utils.file_utils import print_global_struct
from gnomad.utils.reference_genome import get_reference_genome

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    CUSTOM_VEP_ANNOTATIONS,
    RELEASE_KEY_ORDER,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import (
    filter_for_test,
    get_adj_r_ht,
    get_syn_adj_r_ht,
)
from gnomad_constraint.utils.constraint import (
    aggregate_by_constraint_groups,
    aggregate_per_variant_expected_ht,
    calculate_mu_by_downsampling,
    compute_constraint_metrics,
    compute_gene_quality_metrics,
    create_aggregated_expected_ht,
    create_per_variant_expected_ht,
    create_training_set,
    flatten_release_ht,
    prepare_context_ht,
    prepare_ht_for_constraint_calculations,
    prepare_release_ht,
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
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    directory_post_fix = args.directory_post_fix
    path_post_fix = args.path_post_fix
    overwrite = args.overwrite
    custom_vep_annotation = args.custom_vep_annotation
    skip_coverage_model = args.skip_coverage_model
    log10_coverage = args.use_logarithmic_coverage_model

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    if version == "2.1.1":
        raise ValueError(
            "Version 2.1.1 is no longer supported by this constraint pipeline script."
            "Please refer to Commit 39928d1 for the last version of the script that"
            "supports v2.1.1."
        )

    # Generate both "plateau" and "coverage" models unless specified to skip the
    # coverage model.
    models = ["plateau", "coverage"] if not skip_coverage_model else ["plateau"]

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = constraint_res.get_constraint_resources(
        version,
        custom_vep_annotation,
        overwrite,
        test,
        models,
        directory_post_fix,
        path_post_fix,
    )

    try:
        if args.prepare_context_ht:
            logger.info(
                "Annotating methylation, coverage, and GERP on the VEP context Table..."
            )
            res = resources.prepare_context
            res.check_resource_existence()

            # We use naive_coalesce on the context Table because it has a large
            # number of partitions which caused issues with Hail 0.2.133.
            ht = res.context_ht.ht().naive_coalesce(5000)
            if test:
                ht = filter_for_test(ht, use_gene_list=test_gene_list)

            dts = ["exomes", "genomes"]
            ht = prepare_context_ht(
                ht,
                coverage_hts={d: getattr(res, f"{d}_coverage_ht").ht() for d in dts},
                an_hts={d: getattr(res, f"{d}_an_ht").ht() for d in dts},
                freq_hts={
                    d: getattr(res, f"{d}_sites_ht").ht().select("freq") for d in dts
                },
                filter_hts={
                    d: getattr(res, f"{d}_sites_ht").ht().select("filters") for d in dts
                },
                methylation_ht=res.methylation_ht.ht(),
                gerp_ht=constraint_res.get_gerp_ht(get_reference_genome(ht.locus).name),
                adj_r_ht=get_adj_r_ht(),
                syn_adj_r_ht=get_syn_adj_r_ht(),
            )
            ht.write(res.annotated_context_ht.path, overwrite)

            logger.info("Done annotating the VEP context Table.")

        if args.compute_gene_quality_metrics:
            logger.info("Computing per-transcript gene quality metrics...")
            res = resources.compute_gene_quality_metrics
            res.check_resource_existence()

            gencode_cds_ht = constraint_res.get_gencode_cds_ht(version).ht()
            exomes_sites_ht = res.exomes_sites_ht.ht()
            if test:
                gencode_cds_ht = filter_for_test(
                    gencode_cds_ht, use_gene_list=test_gene_list
                )
                exomes_sites_ht = filter_for_test(
                    exomes_sites_ht, use_gene_list=test_gene_list
                )
            gene_quality_ht = compute_gene_quality_metrics(
                res.temp_preprocess_data_ht.ht(),
                exomes_sites_ht,
                gencode_cds_ht,
            )
            gene_quality_ht.write(res.gene_quality_metrics_ht.path, overwrite=overwrite)
            logger.info("Done computing gene quality metrics.")

        if args.calculate_gerp_cutoffs:
            logger.warning(
                "Calculating new GERP cutoffs to be used instead of"
                " '--gerp-lower-cutoff' and '--gerp-upper-cutoff' defaults."
            )
            res = resources.calculate_gerp_cutoffs
            res.check_resource_existence()
            ht = res.annotated_context_ht.ht()
            gerp_lower_cutoff, gerp_upper_cutoff = calculate_gerp_cutoffs(
                ht.filter(ht.genomic_region == "autosome_par")
            )
            logger.info(
                "Calculated new GERP cutoffs: using a lower GERP cutoff of %f "
                "and an upper GERP cutoff of %f.",
                gerp_lower_cutoff,
                gerp_upper_cutoff,
            )

        if args.preprocess_data:
            logger.info(
                "Preprocessing the context Table for all downstream constraint steps..."
            )
            res = resources.preprocess_data
            res.check_resource_existence()
            ht = res.annotated_context_ht.ht()
            ht = filter_for_test(ht, use_gene_list=test_gene_list) if test else ht
            ht = prepare_ht_for_constraint_calculations(
                ht,
                exome_coverage_metric=args.exome_coverage_metric,
                gen_ancs=args.genetic_ancestry_groups,
                include_downsamplings=args.include_downsamplings,
                calculate_mutation_rate_min_cov=args.calculate_mutation_rate_min_cov,
                calculate_mutation_rate_max_cov=args.calculate_mutation_rate_max_cov,
                calculate_mutation_rate_gerp_lower_cutoff=args.calculate_mutation_rate_gerp_lower_cutoff,
                calculate_mutation_rate_gerp_upper_cutoff=args.calculate_mutation_rate_gerp_upper_cutoff,
                max_af=args.max_af,
                build_model_low_cov_cutoff=args.pipeline_low_coverage_filter,
                build_model_high_cov_cutoff=args.build_model_high_cov_definition,
                build_model_upper_cov_cutoff=args.build_model_upper_cov_cutoff,
                apply_model_low_cov_cutoff=args.pipeline_low_coverage_filter,
                apply_model_high_cov_cutoff=args.apply_model_high_cov_definition,
                skip_coverage_model=skip_coverage_model,
            )
            ht.write(res.temp_preprocess_data_ht.path, overwrite=overwrite)

            logger.info("Done preprocessing the context Table.")

        if args.calculate_mutation_rate:
            logger.info("Calculating mutation rate...")
            res = resources.calculate_mutation_rate
            res.check_resource_existence()

            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            ht = calculate_mu_by_downsampling(res.temp_preprocess_data_ht.ht())
            ht = ht.repartition(args.mutation_rate_partitions)
            ht.write(res.mutation_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info("Done calculating mutation rate.")

        if args.create_training_set:
            logger.info(
                "Computing the observed and possible counts of synonymous variants to"
                "use as a training set for the plateau and coverage models..."
            )
            res = resources.create_training_set
            res.check_resource_existence()

            ht = create_training_set(
                res.temp_preprocess_data_ht.ht(),
                res.mutation_ht.ht(),
                partition_hint=args.training_set_partition_hint,
            )

            # TODO: Remove repartition once partition_hint bugs are resolved.
            ht = ht.repartition(args.training_set_partition_hint)
            ht = ht.checkpoint(res.train_ht.path, overwrite=overwrite)
            ht.export(res.train_tsv)

            logger.info("Done with creating training dataset.")

        if args.build_models:
            logger.info("Building plateau and coverage models...")
            res = resources.build_models
            res.check_resource_existence()
            ht = res.train_ht.ht()
            print_global_struct(ht.build_models_globals)
            coverage_model, plateau_models = build_models(
                ht,
                ht.exomes_coverage,
                model_group_expr=ht.build_model,
                skip_coverage_model=skip_coverage_model,
                log10_coverage=log10_coverage,
            )
            hl.experimental.write_expression(
                plateau_models, res.model_plateau.path, overwrite=overwrite
            )
            if not args.skip_coverage_model:
                hl.experimental.write_expression(
                    coverage_model, res.model_coverage.path, overwrite=overwrite
                )

            logger.info("Done building models.")

        if args.apply_models_per_variant:
            logger.info(
                "Applying plateau and coverage models (if specified) per variant to "
                "compute the per-variant expected variant count..."
            )
            res = resources.apply_models_per_variant
            res.check_resource_existence()

            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            ht = res.temp_preprocess_data_ht.ht()
            print_global_struct(ht.apply_models_globals)
            ht = create_per_variant_expected_ht(
                ht,
                res.mutation_ht.ht().select("mu_snp"),
                res.model_plateau.he(),
                coverage_model=None if skip_coverage_model else res.model_coverage.he(),
                log10_coverage=log10_coverage,
                custom_vep_annotation=custom_vep_annotation,
                use_mane_select=True,
            )
            ht.write(res.per_variant_apply_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info("Done computing per-variant expected variant count.")

        if args.aggregate_per_variant_expected:
            logger.info(
                "Aggregating per-variant expected variant count by transcript, "
                "consequence annotations, and consequence modifier annotations..."
            )
            res = resources.aggregate_per_variant_expected
            res.check_resource_existence()

            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            ht = res.per_variant_apply_ht.ht()
            ht = aggregate_per_variant_expected_ht(ht)
            ht.write(res.apply_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info(
                "Done aggregating per-variant expected variant count by transcript, "
                "consequence annotations, and consequence modifier annotations."
            )

        if args.apply_models_aggregated:
            logger.info("Aggregating counts and applying models on aggregated data...")
            res = resources.apply_models_aggregated
            res.check_resource_existence()

            hl._set_flags(use_new_shuffle="1")

            ht = res.temp_preprocess_data_ht.ht()
            print_global_struct(ht.apply_models_globals)
            ht = create_aggregated_expected_ht(
                ht,
                res.mutation_ht.ht().select("mu_snp"),
                res.model_plateau.he(),
                coverage_model=(
                    None if skip_coverage_model else res.model_coverage.he()
                ),
                log10_coverage=log10_coverage,
                custom_vep_annotation=custom_vep_annotation,
                use_mane_select=True,
            )
            ht.write(res.aggregated_expected_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info("Done with aggregated model application.")

        if args.aggregate_by_constraint_groups:
            logger.info(
                "Aggregating observed and expected variant counts by constraint groups..."
            )
            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            if args.use_aggregated_expected:
                res = resources.apply_models_aggregated
                ht = res.aggregated_expected_ht.ht()
            else:
                res = resources.aggregate_by_constraint_groups
                res.check_resource_existence()
                ht = res.apply_ht.ht()

            out_res = resources.aggregate_by_constraint_groups
            aggregate_by_constraint_groups(
                ht,
                keys=tuple(k for k in ht.key if k in RELEASE_KEY_ORDER),
            ).write(out_res.constraint_group_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)
            logger.info("Done with aggregating by constraint groups.")

        if args.compute_constraint_metrics:
            logger.info(
                "Computing constraint metrics, including pLI scores, z scores, oe"
                " ratio, and confidence interval around oe ratio..."
            )
            res = resources.compute_constraint_metrics
            res.check_resource_existence()

            # Compute constraint metrics.
            ht = res.constraint_group_ht.ht(read_args={"_n_partitions": 10000})
            compute_constraint_metrics(
                ht=ht,
                gencode_ht=constraint_res.get_gencode_ht(version),
                gene_quality_metrics_ht=res.gene_quality_metrics_ht.ht(),
                expected_values={
                    "Null": args.expectation_null,
                    "Rec": args.expectation_rec,
                    "LI": args.expectation_li,
                },
                min_diff_convergence=args.min_diff_convergence,
                raw_z_outlier_threshold_lower_lof=args.raw_z_outlier_threshold_lower_lof,
                raw_z_outlier_threshold_lower_missense=args.raw_z_outlier_threshold_lower_missense,
                raw_z_outlier_threshold_lower_syn=args.raw_z_outlier_threshold_lower_syn,
                raw_z_outlier_threshold_upper_syn=args.raw_z_outlier_threshold_upper_syn,
            ).write(res.constraint_metrics_ht.path, overwrite=overwrite)
            logger.info("Done with computing constraint metrics.")

        if args.prepare_release:
            logger.info("Preparing constraint metrics Table for release...")
            res = resources.prepare_release
            res.check_resource_existence()

            constraint_ht = res.constraint_metrics_ht.ht()

            release_ht = prepare_release_ht(
                constraint_ht,
                release_version=args.release_version,
            )
            release_ht.write(res.release_ht.path, overwrite=overwrite)
            logger.info("Done preparing release Table.")

        if args.export_release_tsv or args.export_release_downsampling_tsv:
            res = resources.export_release_tsv
            res.check_resource_existence()
            release_ht = hl.read_table(res.release_ht.path)

            if args.export_release_tsv:
                logger.info("Exporting release TSV...")
                flatten_release_ht(release_ht).export(res.release_tsv)
                logger.info("Done exporting release TSV.")

            if args.export_release_downsampling_tsv:
                logger.info("Exporting release downsampling TSV...")
                downsampling_ht = explode_downsamplings_oe(
                    release_ht,
                    downsampling_meta=hl.eval(release_ht.downsamplings),
                    metrics=["syn", "mis", "lof_hc_lc", "lof"],
                )
                downsampling_ht.export(res.release_downsampling_tsv)
                logger.info("Done exporting release downsampling TSV.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path("constraint_pipeline", version=version)
        )


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
        "--directory-post-fix",
        help="Post-fix to append to the output directory path.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--path-post-fix",
        help="Post-fix to append to the output file path.",
        type=str,
        default=None,
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
        "--test-gene-list",
        help=(
            "Whether to filter the exome Table, genome Table and the context Table to"
            " only a list of genes on chromosome 20, chromosome X, and chromosome Y for"
            " testing."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--prepare-context-ht",
        help=(
            "Prepare the context Table by splitting multiallelic sites and adding "
            "'methylation', 'coverage', and 'gerp' annotations to context Table with "
            "VEP annotation."
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

    preprocess_args = parser.add_argument_group(
        "Preprocess data args",
        "All arguments used for preprocessing data for downstream steps.",
    )
    preprocess_args.add_argument(
        "--preprocess-data",
        help="Preprocess the context Table for downstream constraint steps.",
        action="store_true",
    )
    preprocess_args.add_argument(
        "--calculate-mutation-rate-min-cov",
        help=(
            "Minimum coverage required to keep a site when calculating the mutation"
            " rate. Default is 15."
        ),
        type=int,
        default=15,
    )
    preprocess_args.add_argument(
        "--calculate-mutation-rate-max-cov",
        help=(
            "Maximum coverage required to keep a site when calculating the mutation"
            " rate. Default is 60."
        ),
        type=int,
        default=60,
    )
    preprocess_args.add_argument(
        "--calculate-mutation-rate-gerp-lower-cutoff",
        help=(
            "Minimum GERP score for variant to be included when calculating the"
            " mutation rate. Default is -3.9885 (precalculated on the GRCh37 context"
            " Table to define the 95th percentile)."
        ),
        type=float,
        default=-3.9885,
    )
    preprocess_args.add_argument(
        "--calculate-mutation-rate-gerp-upper-cutoff",
        help=(
            "Maximum GERP score for variant to be included when calculating the"
            " mutation rate. Default is 2.6607 (precalculated on the GRCh37 context"
            " Table to define the 5th percentile)."
        ),
        type=float,
        default=2.6607,
    )
    preprocess_args.add_argument(
        "--exome-coverage-metric",
        help=(
            "Name of metric to use to assess exome coverage, such as 'median', 'AN', or"
            "'AN_percent'. Default is 'AN_percent'."
        ),
        type=str,
        default="AN_percent",
    )
    preprocess_args.add_argument(
        "--pipeline-low-coverage-filter",
        help=(
            "Lower exome coverage cutoff to use throughout the pipeline. Sites with"
            " coverage below this cutoff will be excluded when creating the training"
            " set, building and applying models, and computing constraint metrics."
            "  Default is None."
        ),
        type=int,
        default=None,
    )
    preprocess_args.add_argument(
        "--max-af",
        help=(
            "Maximum variant allele frequency to use when filtering variants for "
            "training and applying models."
        ),
        type=float,
        default=0.001,
    )
    preprocess_args.add_argument(
        "--genetic-ancestry-groups",
        nargs="+",
        help=(
            "Populations on which to build models, apply models, and or compute metrics "
            "on. Default is None."
        ),
        choices=["afr", "amr", "eas", "nfe", "sas"],
        default=None,
    )
    preprocess_args.add_argument(
        "--include-downsamplings",
        help="Include downsamplings in the constraint pipeline.",
        action="store_true",
    )
    preprocess_args.add_argument(
        "--skip-coverage-model",
        help="Omit computing and applying the coverage model.",
        action="store_true",
    )
    preprocess_args.add_argument(
        "--build-model-upper-cov-cutoff",
        help=(
            "Upper exome coverage cutoff. Sites with coverage above this cutoff are"
            " excluded from the high coverage Table when building the models. Default"
            " is None."
        ),
        type=int,
        default=None,
    )
    preprocess_args.add_argument(
        "--build-model-high-cov-definition",
        help=(
            "Lower exome coverage cutoff to use to define high coverage sites when "
            "building models. Sites with coverage below this cutoff are excluded from "
            "the high coverage Table when building models. Default is 90."
        ),
        type=int,
        default=90,
    )
    preprocess_args.add_argument(
        "--apply-model-high-cov-definition",
        help=(
            "Lower exome coverage cutoff to use to define high coverage sites when "
            "applying models. Sites with coverage below this cutoff are excluded from "
            "the high coverage Table when applying models. Default is 90."
        ),
        type=int,
        default=90,
    )

    mutation_rate_args = parser.add_argument_group(
        "Calculate mutation rate args",
        "Arguments used for calculating the mutation rate.",
    )
    mutation_rate_args.add_argument(
        "--calculate-mutation-rate",
        help=(
            "Calculate baseline mutation rate for each substitution and context using"
            " downsampling data."
        ),
        action="store_true",
    )
    mutation_rate_args.add_argument(
        "--mutation-rate-partitions",
        help=(
            "Number of partitions to which the mutation rate Table should be "
            "repartitioned."
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
            "Count the observed variants and possible variants by exome coverage at "
            "synonymous sites."
        ),
        action="store_true",
    )
    training_set_args.add_argument(
        "--training-set-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting variants for "
            "training datasets."
        ),
        type=int,
        default=100,
    )

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
    cov_model_type = build_models_args.add_argument(
        "--use-logarithmic-coverage-model",
        help=(
            "Use a logarithmic model for low coverage sites when building and applying "
            "the coverage model."
        ),
        action="store_true",
    )

    parser.add_argument(
        "--apply-models-per-variant",
        help=(
            "Apply plateau and coverage models to variants in exome sites Table and"
            " context Table to compute expected variant counts per variant."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-models-aggregated",
        help=(
            "Apply plateau and coverage models and aggregate to constraint groups"
            " in a single step, without writing per-variant intermediates. This is"
            " an alternative to running --apply-models-per-variant,"
            " --aggregate-per-variant-expected, and"
            " --aggregate-by-constraint-groups separately."
        ),
        action="store_true",
    )

    aggregate_per_variant_expected_args = parser.add_argument_group(
        "Aggregate per variant expected args",
        "Arguments used for applying aggregating the per variant expected values.",
    )
    aggregate_per_variant_expected_args.add_argument(
        "--aggregate-per-variant-expected",
        help=(
            "Aggregate the per-variant expected variant counts to get the expected "
            "variant counts for each transcript by consequence annotation and "
            "modifier."
        ),
        action="store_true",
    )

    aggregate_per_variant_expected_args.add_argument(
        "--apply-obs-pos-count-partition-hint",
        help=(
            "Target number of partitions for aggregation when counting observed and"
            " expected variants for model application"
        ),
        type=int,
        default=2000,
    )
    aggregate_per_variant_expected_args.add_argument(
        "--apply-expected-variant-partition-hint",
        help=(
            "Target number of partitions for sum aggregators after applying models to"
            " get expected variant counts."
        ),
        type=int,
        default=1000,
    )
    aggregate_per_variant_expected_args.add_argument(
        "--custom-vep-annotation",
        help=(
            "Custom VEP annotation to be used to annotate transcript when"
            ' applying models (one of "transcript_consequences" or "worst_csq_by_gene")'
        ),
        type=str,
        default="transcript_consequences",
        choices=CUSTOM_VEP_ANNOTATIONS,
    )
    aggregate_per_variant_expected_args._group_actions.append(cov_model_type)

    aggregate_by_constraint_groups_args = parser.add_argument_group(
        "Aggregate by constraint groups args",
        "Arguments used for aggregating by constraint groups.",
    )
    aggregate_by_constraint_groups_args.add_argument(
        "--aggregate-by-constraint-groups",
        help=(
            "Aggregate the observed and expected variant counts by constraint groups"
            " to get the constraint metrics."
        ),
        action="store_true",
    )
    aggregate_by_constraint_groups_args.add_argument(
        "--use-aggregated-expected",
        help=(
            "Read from the --apply-models-aggregated output instead of the"
            " --aggregate-per-variant-expected output as input for"
            " --aggregate-by-constraint-groups."
        ),
        action="store_true",
    )

    gene_quality_args = parser.add_argument_group(
        "Compute gene quality metrics args",
        "Arguments used for computing per-transcript gene quality metrics.",
    )
    gene_quality_args.add_argument(
        "--compute-gene-quality-metrics",
        help=(
            "Compute per-transcript gene quality metrics (coverage, mapping quality,"
            " segdup, LCR) from the preprocessed context Table and gnomAD exomes"
            " sites Table."
        ),
        action="store_true",
    )
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
            "Number of partitions to which the Table of expected variant counts should "
            "be reaprtitioned."
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
            " selection. Default is 1.0"
        ),
        type=float,
        default=1.0,
    )
    compute_constraint_args.add_argument(
        "--expectation-rec",
        help=(
            "Expected observed/expected rate of truncating variation for recessive"
            " disease genes. Default is 0.706. Note that v2 used a value of 0.463."
        ),
        type=float,
        default=0.706,
    )
    compute_constraint_args.add_argument(
        "--expectation-li",
        help=(
            "Expected observed/expected rate of truncating variation for severe"
            " haploinsufficient genes. Default is 0.207. Note that v2 used a value of"
            " 0.089."
        ),
        type=float,
        default=0.207,
    )
    compute_constraint_args.add_argument(
        "--raw-z-outlier-threshold-lower-lof",
        help=(
            "Value at which the raw z-score is considered an outlier for lof variants."
            " Values below this threshold will be considered outliers. Default is -8.0."
        ),
        type=float,
        default=-8.0,
    )
    compute_constraint_args.add_argument(
        "--raw-z-outlier-threshold-lower-missense",
        help=(
            "Value at which the raw z-score is considered an outlier for missense"
            " variants. Values below this threshold will be considered outliers."
            " Default is -8.0."
        ),
        type=float,
        default=-8.0,
    )
    compute_constraint_args.add_argument(
        "--raw-z-outlier-threshold-lower-syn",
        help=(
            "Lower value at which the raw z-score is considered an outlier for"
            " synonymous variants. Values below this threshold will be considered"
            " outliers. Default is -8.0."
        ),
        type=float,
        default=-8.0,
    )
    compute_constraint_args.add_argument(
        "--raw-z-outlier-threshold-upper-syn",
        help=(
            "Upper value at which the raw z-score is considered an outlier for"
            " synonymous variants. Values above this threshold will be considered"
            " outliers. Default is 8.0."
        ),
        type=float,
        default=8.0,
    )
    prepare_release_args = parser.add_argument_group(
        "Prepare release args",
        "Arguments used for preparing the constraint metrics Table for release.",
    )
    prepare_release_args.add_argument(
        "--prepare-release",
        help=(
            "Prepare the constraint metrics Table for public release by restructuring "
            "constraint groups into named top-level fields and consolidating globals."
        ),
        action="store_true",
    )
    prepare_release_args.add_argument(
        "--release-version",
        help=(
            "Version string to set in the release Table globals. If not specified, "
            "the existing version global is retained."
        ),
        type=str,
        default=None,
    )
    prepare_release_args.add_argument(
        "--export-release-tsv",
        help=(
            "Flatten the release Hail Table and export it as a TSV. Output paths are"
            " determined by the release resource functions."
        ),
        action="store_true",
    )
    prepare_release_args.add_argument(
        "--export-release-downsampling-tsv",
        help=(
            "Export per-genetic-ancestry downsampling observed and expected counts from"
            " the release Hail Table as a TSV. Reads from the release HT path."
        ),
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
