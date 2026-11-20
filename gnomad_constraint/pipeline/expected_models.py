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
      and a sampled background of non-mutated synonymous sites, in the plateau
      models' high coverage range, with their window's genomic features, leaving out
      genes with a QC gene flag.
    - --screen-rr-features: optionally, test each feature of each context on its
      own, for Gnocchi's feature selection (--rr-feature-selection).
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
    screen_regional_mutation_rate_features,
    tune_regional_mutation_rate_penalty,
)
from gnomad.utils.file_utils import print_global_struct
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_BONFERRONI_CUTOFF,
    RR_CONTEXT_FEATURES,
    RR_FEATURE_SELECTION_ALL,
    RR_GENE_FLAGS_VERSION,
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
    select_rr_features,
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


def rr_feature_selection_arg(value: str) -> str:
    """
    Check an ``--rr-feature-selection`` value.

    :param value: ``RR_FEATURE_SELECTION_ALL``, ``RR_BONFERRONI_CUTOFF``, or a p-value
        cutoff above 0 and at most 1.
    :return: The value.
    """
    if value in (RR_FEATURE_SELECTION_ALL, RR_BONFERRONI_CUTOFF):
        return value
    try:
        cutoff = float(value)
    except ValueError:
        cutoff = None
    if cutoff is None or not 0 < cutoff <= 1:
        raise argparse.ArgumentTypeError(
            f"must be {RR_FEATURE_SELECTION_ALL!r}, {RR_BONFERRONI_CUTOFF!r}, or a"
            f" p-value cutoff in (0, 1], not {value!r}"
        )
    return value


def check_rr_input_record(
    ht: hl.Table, input_ht: hl.Table, name: str, step: str
) -> None:
    """
    Check that an rr step's output was made from its input Table as it is now.

    The step that wrote ``ht`` records its input's globals and row count as the
    ``{name}_globals`` and ``{name}_row_count`` globals. If the input has other
    globals or another row count now, it was rewritten after the step ran, so ``ht``
    describes other data. The globals record how the input was made (for the
    training Table, its features and feature selection, which don't change its row
    count).

    :param ht: Output of ``step``, with the record of its input.
    :param input_ht: The input Table, as read from disk, so counting its rows reads
        only metadata.
    :param name: Prefix of the record's globals.
    :param step: Flag of the step that wrote ``ht``, for the error message.
    :return: None.
    """
    if f"{name}_globals" not in ht.globals or f"{name}_row_count" not in ht.globals:
        raise ValueError(
            f"The {step} output doesn't record the Table it was made from. Rerun"
            f" {step}."
        )
    recorded = hl.eval(ht[f"{name}_row_count"])
    n_rows = input_ht.count()
    if recorded != n_rows:
        raise ValueError(
            f"{step} ran on a Table with {recorded} rows, but it has {n_rows} now."
            f" Rerun {step}."
        )
    recorded = dict(hl.eval(ht[f"{name}_globals"]))
    now = dict(hl.eval(input_ht.globals))
    changed = sorted(
        k for k in recorded.keys() | now.keys() if recorded.get(k) != now.get(k)
    )
    if changed:
        raise ValueError(
            f"{step} ran on a Table with other {', '.join(changed)} globals. Rerun"
            f" {step}."
        )


def get_rr_ridge_penalty(tuning_ht: hl.Table, training_ht: hl.Table) -> float:
    """
    Get the tuned ridge penalty, checking it was tuned on this rr training Table.

    ``--tune-rr-penalty`` records the training Table's globals and row count
    (:func:`check_rr_input_record`), so a training Table rebuilt since, with other
    sites or another ``--rr-feature-selection``, is caught.

    :param tuning_ht: Penalty tuning Table with ``ridge_penalty``,
        ``rr_training_globals``, and ``rr_training_row_count`` globals.
    :param training_ht: rr training Table, as read from disk.
    :return: The tuned ridge penalty.
    """
    check_rr_input_record(tuning_ht, training_ht, "rr_training", "--tune-rr-penalty")

    return hl.eval(tuning_ht.ridge_penalty)


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
        rr_feature_screen=args.rr_feature_selection != RR_FEATURE_SELECTION_ALL,
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
                res.temp_preprocess_data_ht.ht(),
                res.dnm_ht.ht(),
                res.genomic_features_ht.ht(),
                RR_BACKGROUND_FRACTION,
                RR_BACKGROUND_SEED,
                gene_flags_ht=(
                    None if args.skip_qc_flagged_gene_filter else res.gene_flags_ht.ht()
                ),
                gene_flags_version=RR_GENE_FLAGS_VERSION,
            )
            # The sites are a few hundred thousand of the whole context
            # Table. naive_coalesce merges upstream partitions, so coalescing before
            # the checkpoint would run that scan in only RR_TRAINING_N_PARTITIONS
            # tasks; coalescing after it only merges the small, mostly empty result.
            ht = ht.checkpoint(new_temp_file("rr_site_features", "ht"))
            ht = ht.naive_coalesce(RR_TRAINING_N_PARTITIONS)
            ht.write(res.rr_site_features_ht.path, overwrite=overwrite)
            logger.info("Done creating the rr site features.")

        if args.screen_rr_features:
            logger.info("Screening the rr features one at a time...")
            res = resources.screen_rr_features
            res.check_resource_existence()

            site_ht = res.rr_site_features_ht.ht()
            ht = screen_regional_mutation_rate_features(
                prepare_regional_mutation_rate_training_ht(site_ht, RR_CONTEXT_FEATURES)
            )
            # --create-rr-training-set checks the screen is of its site features.
            ht = ht.annotate_globals(
                rr_site_features_globals=site_ht.index_globals(),
                rr_site_features_row_count=site_ht.count(),
            )
            ht.write(res.rr_feature_screen_ht.path, overwrite=overwrite)
            logger.info("Done screening the rr features.")

        if args.create_rr_training_set:
            logger.info("Creating the rr training set...")
            res = resources.create_rr_training_set
            res.check_resource_existence()

            site_ht = res.rr_site_features_ht.ht()
            features = RR_CONTEXT_FEATURES
            if args.rr_feature_selection != RR_FEATURE_SELECTION_ALL:
                screen_ht = res.rr_feature_screen_ht.ht()
                check_rr_input_record(
                    screen_ht, site_ht, "rr_site_features", "--screen-rr-features"
                )
                features = select_rr_features(screen_ht, args.rr_feature_selection)
                logger.info(
                    "Features selected by %s: %s", args.rr_feature_selection, features
                )
                no_features = sorted(set(RR_CONTEXT_FEATURES) - set(features))
                if no_features:
                    logger.warning(
                        "No feature was selected for %s, so they get no rr model.",
                        no_features,
                    )
            ht = prepare_regional_mutation_rate_training_ht(
                site_ht, features, fold_expr=get_rr_fold_expr(site_ht.locus)
            )
            # gnomad_methods keeps only the features global, so carry the record of
            # how the sites and features were chosen.
            ht = ht.annotate_globals(
                rr_training_sites=site_ht.index_globals().rr_training_sites,
                rr_feature_selection=args.rr_feature_selection,
            )
            ht.write(res.rr_training_ht.path, overwrite=overwrite)
            logger.info("Done creating the rr training set.")

        if args.tune_rr_penalty:
            logger.info("Tuning the rr ridge penalty over %s...", args.rr_penalties)
            res = resources.tune_rr_penalty
            res.check_resource_existence()

            training_ht = res.rr_training_ht.ht()
            ht = tune_regional_mutation_rate_penalty(training_ht, args.rr_penalties)
            penalty = choose_regional_mutation_rate_penalty(ht, RR_RIDGE_PENALTY)
            logger.info("Chosen ridge penalty: %s", penalty)
            # --train-rr-models checks the penalty was tuned on its training Table.
            ht = ht.annotate_globals(
                ridge_penalty=penalty,
                rr_training_globals=training_ht.index_globals(),
                rr_training_row_count=training_ht.count(),
            )
            ht.write(res.rr_penalty_tuning_ht.path, overwrite=overwrite)
            logger.info("Done tuning the rr ridge penalty.")

        if args.train_rr_models:
            res = resources.train_rr_models
            res.check_resource_existence()
            tuning_ht = res.rr_penalty_tuning_ht.ht()
            training_ht = res.rr_training_ht.ht()
            penalty = get_rr_ridge_penalty(tuning_ht, training_ht)
            logger.info("Training the rr models (ridge penalty %s)...", penalty)

            # The fit records the penalty in a ridge_penalty global.
            ht = fit_regional_mutation_rate_models(training_ht, penalty)
            training_globals = training_ht.index_globals()
            ht = ht.annotate_globals(
                rr_training_sites=training_globals.rr_training_sites,
                rr_feature_selection=training_globals.rr_feature_selection,
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
            " --create-training-set onward (training set, models, rr Table, expected"
            " counts, constraint groups, constraint metrics, and the release"
            " constraint metrics Table and TSVs), but not to the annotated context,"
            " preprocessed, gene quality metrics, mutation rate, rr model, or release"
            " mutation rate outputs."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to read and write the test Tables in the testing bucket, made"
            " from the subset of the data that prepare_data.py --test or"
            " --test-gene-list keeps."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-gene-list",
        help=(
            "Same as --test in this script; prepare_data.py --test-gene-list applies"
            " the six-gene filter."
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
        "--screen-rr-features",
        help=(
            "Test each rr feature of each context on its own (an unpenalized"
            " univariate logistic fit), for --rr-feature-selection."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--create-rr-training-set",
        help=(
            "Prepare the rr site features for the models: each site's label, context,"
            " chromosome fold, and its context's features in RR_CONTEXT_FEATURES, or"
            " those --rr-feature-selection selects."
        ),
        action="store_true",
    )
    rr_model_args.add_argument(
        "--rr-feature-selection",
        help=(
            "Features each context's rr model uses, set at --create-rr-training-set:"
            f" {RR_FEATURE_SELECTION_ALL!r} for every feature in RR_CONTEXT_FEATURES,"
            f" or Gnocchi's selection from the --screen-rr-features output, with"
            f" {RR_BONFERRONI_CUTOFF!r} for its Bonferroni cutoff or a p-value cutoff"
            " such as 1e-2; only the smallest selected scale of each feature is kept."
            " To fit the selected features without a penalty, as Gnocchi did, also"
            " pass --rr-penalties 0. Default is"
            f" {RR_FEATURE_SELECTION_ALL!r}."
        ),
        type=rr_feature_selection_arg,
        default=RR_FEATURE_SELECTION_ALL,
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
