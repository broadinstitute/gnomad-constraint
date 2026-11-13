"""
Fit the models that expected variant counts come from.

The steps, in order, after prepare_data.py:
    - --calculate-mutation-rate: the mutation rate for each substitution and
      context.
    - --create-training-set: count observed and possible synonymous variants for
      model training.
    - --build-models: build the plateau and coverage models.

And the regional mutation rate (rr) model steps, which need only the
prepare_data.py outputs:
    - --create-rr-site-features: create the synonymous de novo mutation (DNM) sites
      and a sampled background of non-mutated synonymous sites, with their window's
      genomic features, leaving out genes with a QC gene flag.
    - --create-rr-training-set: prepare the sites for the models: each site's
      label, context, chromosome fold, and its context's features.
    - --tune-rr-penalty: choose the ridge penalty by the held-out log loss of
      chromosome-split models.
    - --train-rr-models: fit the per-context rr models.

Then run expected_counts.py.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.constraint import (
    build_models,
    choose_regional_mutation_rate_penalty,
    fit_regional_mutation_rate_models,
    prepare_regional_mutation_rate_training_ht,
    tune_regional_mutation_rate_penalty,
)
from gnomad.utils.file_utils import print_global_struct
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_CONTEXT_FEATURES,
    RR_RIDGE_PENALTIES,
    RR_RIDGE_PENALTY,
    RR_TRAINING_N_PARTITIONS,
    VERSIONS,
)
from gnomad_constraint.utils.constraint import (
    calculate_mu_by_downsampling,
    create_training_set,
)
from gnomad_constraint.utils.regional_rate import (
    create_rr_site_features_ht,
    get_rr_fold_expr,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("expected_models")
logger.setLevel(logging.INFO)


def check_skip_coverage_model(
    params_expr: hl.expr.StructExpression, skip_coverage_model: bool
) -> None:
    """
    Check that ``--skip-coverage-model`` matches the setting used in preprocessing.

    Preprocessing assigns low coverage sites to the coverage model unless it is
    skipped, so building or applying the models with a different setting gives wrong
    expected counts.

    :param params_expr: Parameter global of the input Table (``build_models_globals``
        or ``apply_models_globals``).
    :param skip_coverage_model: Whether ``--skip-coverage-model`` was passed.
    :return: None.
    """
    preprocess_skip = hl.eval(params_expr.skip_coverage_model)
    if preprocess_skip != skip_coverage_model:
        raise ValueError(
            f"--skip-coverage-model is {skip_coverage_model}, but the input Table was"
            f" preprocessed with skip_coverage_model={preprocess_skip}. Use the same"
            " setting as --preprocess-data."
        )


def main(args: argparse.Namespace) -> None:
    """Execute the expected models steps of the constraint pipeline."""
    hl.init(
        log="/expected_models.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    directory_post_fix = args.directory_post_fix
    path_post_fix = args.path_post_fix
    overwrite = args.overwrite
    skip_coverage_model = args.skip_coverage_model
    log10_coverage = args.use_logarithmic_coverage_model

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    # Generate both "plateau" and "coverage" models unless specified to skip the
    # coverage model.
    models = ["plateau", "coverage"] if not skip_coverage_model else ["plateau"]

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = constraint_res.get_constraint_resources(
        version,
        overwrite,
        test,
        models,
        directory_post_fix,
        path_post_fix,
    )

    try:
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
                " use as a training set for the plateau and coverage models..."
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
            check_skip_coverage_model(ht.build_models_globals, skip_coverage_model)
            coverage_model, plateau_models = build_models(
                ht,
                ht.exomes_coverage,
                weighted=args.use_weights,
                model_group_expr=ht.build_model,
                skip_coverage_model=skip_coverage_model,
                log10_coverage=log10_coverage,
            )
            hl.experimental.write_expression(
                plateau_models, res.model_plateau.path, overwrite=overwrite
            )
            if not args.skip_coverage_model:
                # Save the coverage scale with the model, so it's applied on the same
                # scale.
                hl.experimental.write_expression(
                    hl.struct(
                        coefficients=coverage_model, log10_coverage=log10_coverage
                    ),
                    res.model_coverage.path,
                    overwrite=overwrite,
                )

            logger.info("Done building models.")

        if args.create_rr_site_features:
            logger.info("Creating the rr site features...")
            res = resources.create_rr_site_features
            res.check_resource_existence()

            ht = create_rr_site_features_ht(
                res.annotated_context_ht.ht(),
                res.dnm_ht.ht(),
                res.genomic_features_ht.ht(),
                RR_BACKGROUND_FRACTION,
                RR_BACKGROUND_SEED,
                gene_flags_ht=(
                    None if args.skip_qc_flagged_gene_filter else res.gene_flags_ht.ht()
                ),
            )
            # The sites are a few hundred thousand of the whole context
            # Table. naive_coalesce merges upstream partitions, so coalescing before
            # the checkpoint would run that scan in only RR_TRAINING_N_PARTITIONS
            # tasks; coalescing after it only merges the small, mostly empty result.
            ht = ht.checkpoint(new_temp_file("rr_site_features", "ht"))
            ht = ht.naive_coalesce(RR_TRAINING_N_PARTITIONS)
            ht.write(res.rr_site_features_ht.path, overwrite=overwrite)
            logger.info("Done creating the rr site features.")

        if args.create_rr_training_set:
            logger.info("Creating the rr training set...")
            res = resources.create_rr_training_set
            res.check_resource_existence()

            ht = res.rr_site_features_ht.ht()
            ht = prepare_regional_mutation_rate_training_ht(
                ht, RR_CONTEXT_FEATURES, fold_expr=get_rr_fold_expr(ht.locus)
            )
            ht.write(res.rr_training_ht.path, overwrite=overwrite)
            logger.info("Done creating the rr training set.")

        if args.tune_rr_penalty:
            logger.info("Tuning the rr ridge penalty over %s...", args.rr_penalties)
            res = resources.tune_rr_penalty
            res.check_resource_existence()

            ht = tune_regional_mutation_rate_penalty(
                res.rr_training_ht.ht(), args.rr_penalties
            )
            penalty = choose_regional_mutation_rate_penalty(ht, RR_RIDGE_PENALTY)
            logger.info("Chosen ridge penalty: %s", penalty)
            ht = ht.annotate_globals(ridge_penalty=penalty)
            ht.write(res.rr_penalty_tuning_ht.path, overwrite=overwrite)
            logger.info("Done tuning the rr ridge penalty.")

        if args.train_rr_models:
            res = resources.train_rr_models
            res.check_resource_existence()
            penalty = hl.eval(res.rr_penalty_tuning_ht.ht().ridge_penalty)
            logger.info("Training the rr models (ridge penalty %s)...", penalty)

            ht = fit_regional_mutation_rate_models(res.rr_training_ht.ht(), penalty)
            # The models Table has one small row per context.
            non_converged = ht.aggregate(
                hl.agg.filter(~ht.converged, hl.agg.collect(ht.context))
            )
            ht = ht.annotate_globals(
                ridge_penalty=penalty,
                non_converged_contexts=hl.literal(non_converged, hl.tarray(hl.tstr)),
            )
            ht.write(res.rr_models_ht.path, overwrite=overwrite)
            logger.info("Done training the rr models.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path(
                "expected_models",
                version=version,
                directory_post_fix=directory_post_fix,
                path_post_fix=path_post_fix,
            )
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
        help="Post-fix to append to the output directory path, for every output.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--path-post-fix",
        help=(
            "Post-fix to append to output file names. Applies to the outputs from"
            " --create-training-set onward (training set, models, rr Tables, expected"
            " counts, constraint groups, constraint metrics, and the release"
            " constraint metrics Table and TSVs), but not to the annotated context,"
            " preprocessed, gene quality metrics, mutation rate, or release mutation"
            " rate outputs."
        ),
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
            " only six genes for testing: PCSK9 (chr1), ADNP and ANGPT4 (chr20), OFD1"
            " and FAAH2 (chrX), and ZFY (chrY)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-coverage-model",
        help="Omit computing and applying the coverage model.",
        action="store_true",
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
            " 'possible_variants'."
        ),
        action="store_true",
    )
    build_models_args.add_argument(
        "--use-logarithmic-coverage-model",
        help=(
            "Build the coverage model on log10 coverage instead of linear coverage."
            " The scale is saved with the model, and the model is applied on it. For"
            " a coverage model saved without the scale, pass this flag when applying"
            " to apply it on log10 coverage; an error is raised if its coefficients"
            " look like the other scale."
        ),
        action="store_true",
    )

    rr_model_args = parser.add_argument_group(
        "Regional mutation rate (rr) model args",
        "Arguments used for training the rr models.",
    )
    rr_model_args.add_argument(
        "--create-rr-site-features",
        help=(
            "Create the rr site features: synonymous DNM sites and a sampled"
            " background of non-mutated synonymous sites with their window's"
            " features, leaving out genes with a QC gene flag."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--skip-qc-flagged-gene-filter",
        help=(
            "With --create-rr-site-features, keep the sites in genes with a QC gene"
            " flag in the v4.1.2 constraint metrics, which are left out by default."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--create-rr-training-set",
        help=(
            "Prepare the rr site features for the models: each site's label, context,"
            " chromosome fold, and its context's features in RR_CONTEXT_FEATURES."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--tune-rr-penalty",
        help=(
            "Choose the ridge penalty of the rr models: the penalty in --rr-penalties"
            " with the lowest held-out log loss of models fit on odd chromosomes and"
            " scored on even ones, and vice versa, summed over contexts."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--rr-penalties",
        help=(
            "Per-site ridge penalties tried by --tune-rr-penalty. Default is"
            f" {RR_RIDGE_PENALTIES}."
        ),
        type=float,
        nargs="+",
        default=RR_RIDGE_PENALTIES,
    )
    rr_model_args.add_argument(
        "--train-rr-models",
        help=(
            "Fit the per-context rr models with the penalty chosen by"
            " --tune-rr-penalty."
        ),
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
