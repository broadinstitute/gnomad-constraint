"""
Fit the models that expected variant counts come from.

The steps, in order, after prepare_data.py:
    - --calculate-mutation-rate: the mutation rate for each substitution and
      context.
    - --create-training-set: count observed and possible synonymous variants for
      model training.
    - --build-models: build the plateau and coverage models.

Then run expected_counts.py.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.constraint import build_models
from gnomad.utils.file_utils import print_global_struct

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import CURRENT_VERSION, VERSIONS
from gnomad_constraint.utils.constraint import (
    calculate_mu_by_downsampling,
    create_training_set,
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
            " --create-training-set onward (training set, models, expected counts,"
            " constraint groups, constraint metrics, and the release constraint"
            " metrics Table and TSVs), but not to the annotated context, preprocessed,"
            " gene quality metrics, mutation rate, or release mutation rate outputs."
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

    args = parser.parse_args()
    main(args)
