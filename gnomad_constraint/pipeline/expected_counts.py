"""
Apply the models to compute expected variant counts.

The steps, after expected_models.py:
    - --apply-models-per-variant: the expected count of each variant.
    - --compute-rr: the regional mutation rate (rr) of each 1kb window and context,
      from the expected_models.py rr models, normalized on the
      --apply-models-per-variant output that
      constraint_metrics.py --aggregate-per-variant-expected --apply-regional-rr
      applies it to.
    - --apply-models-aggregated: an alternative to --apply-models-per-variant and
      constraint_metrics.py --aggregate-per-variant-expected: sum the observed and
      possible counts first, then apply the models.

Then run constraint_metrics.py.
"""

import argparse
import logging
import math
from typing import Optional, Tuple

import hail as hl
from gnomad.utils.file_utils import file_exists, print_global_struct
from gnomad_qc.resource_utils import PipelineStepResourceCollection

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.pipeline.expected_models import check_skip_coverage_model
from gnomad_constraint.resources.constants import CURRENT_VERSION, VERSIONS
from gnomad_constraint.utils.constraint import (
    create_aggregated_expected_ht,
    create_per_variant_expected_ht,
)
from gnomad_constraint.utils.regional_rate import compute_rr_ht

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("expected_counts")
logger.setLevel(logging.INFO)


def check_apply_models_resources(
    res: PipelineStepResourceCollection, skip_coverage_model: bool
) -> None:
    """
    Check an apply-models step's ``--skip-coverage-model`` setting, then its resources.

    The setting is checked first, against the preprocessed Table if it exists. With
    a mismatch the coverage model is missing, and the resource check would report it
    as a missing ``--build-models`` output instead of as the mismatch.

    :param res: Resources of ``--apply-models-per-variant`` or
        ``--apply-models-aggregated``.
    :param skip_coverage_model: Whether ``--skip-coverage-model`` was passed.
    :return: None.
    """
    preprocessed_ht = res.temp_preprocess_data_ht
    if file_exists(preprocessed_ht.path):
        check_skip_coverage_model(
            preprocessed_ht.ht().apply_models_globals, skip_coverage_model
        )
    res.check_resource_existence()


def get_coverage_model(
    res: PipelineStepResourceCollection,
    skip_coverage_model: bool,
    log10_coverage: bool,
    high_cov_cutoff: int,
) -> Tuple[Optional[hl.expr.ArrayNumericExpression], bool]:
    """
    Read the coverage model and the coverage scale to apply it on.

    ``--build-models`` saves the coverage model as a struct of its ``coefficients``
    (``[intercept, slope]``) and the ``log10_coverage`` setting it was built with,
    and the model is applied on that scale. A model saved before the setting was
    stored has only the coefficients, so ``--use-logarithmic-coverage-model`` sets
    the scale for it, with a warning.

    For such a model, the coefficients are checked against the flag. The model is
    fit to o/e scaled so high coverage sites are at 1, so on the scale it was built
    on, its correction at the high coverage cutoff is close to 1, and on the other
    scale it is far from 1. If the flag's scale gives the correction further from 1,
    an error is raised.

    :param res: Resources of ``--apply-models-per-variant`` or
        ``--apply-models-aggregated``.
    :param skip_coverage_model: Whether ``--skip-coverage-model`` was passed.
    :param log10_coverage: Whether ``--use-logarithmic-coverage-model`` was passed.
    :param high_cov_cutoff: Coverage at and above which sites are high coverage
        when applying the models.
    :return: Tuple of the coverage model coefficients, or None when the coverage
        model is skipped, and whether to apply them to log10 coverage.
    """
    if skip_coverage_model:
        return None, log10_coverage

    coverage_model = res.model_coverage.he()
    if not isinstance(coverage_model, hl.expr.StructExpression):
        intercept, slope = hl.eval(coverage_model)
        correction = {
            log10: intercept
            + slope * (math.log10(high_cov_cutoff) if log10 else high_cov_cutoff)
            for log10 in (False, True)
        }
        if abs(correction[log10_coverage] - 1) > abs(
            correction[not log10_coverage] - 1
        ):
            scale = {False: "linear", True: "log10"}
            raise ValueError(
                "The coverage model was saved without the coverage scale it was built"
                " on, and its coefficients look like it was built on"
                f" {scale[not log10_coverage]} coverage: at the high coverage cutoff"
                f" ({high_cov_cutoff}), its correction is"
                f" {correction[log10_coverage]:.3g} on {scale[log10_coverage]}"
                f" coverage and {correction[not log10_coverage]:.3g} on"
                f" {scale[not log10_coverage]} coverage, where it should be close to"
                f" 1. {'Remove' if log10_coverage else 'Add'}"
                " --use-logarithmic-coverage-model, or rebuild the models."
            )
        logger.warning(
            "The coverage model was saved without the coverage scale it was built on,"
            " so it is applied to %s coverage, as set by"
            " --use-logarithmic-coverage-model. Rebuild the models to save the scale"
            " with them.",
            "log10" if log10_coverage else "linear",
        )
        return coverage_model, log10_coverage

    saved_log10_coverage = hl.eval(coverage_model.log10_coverage)
    if log10_coverage and not saved_log10_coverage:
        logger.warning(
            "--use-logarithmic-coverage-model is ignored: the coverage model was built"
            " on linear coverage."
        )
    logger.info(
        "Applying the coverage model to %s coverage, the scale it was built on.",
        "log10" if saved_log10_coverage else "linear",
    )
    return coverage_model.coefficients, saved_log10_coverage


def main(args: argparse.Namespace) -> None:
    """Execute the expected counts steps of the constraint pipeline."""
    hl.init(
        log="/expected_counts.log",
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
        if args.apply_models_per_variant:
            logger.info(
                "Applying plateau and coverage models (if specified) per variant to "
                "compute the per-variant expected variant count..."
            )
            res = resources.apply_models_per_variant
            check_apply_models_resources(res, skip_coverage_model)

            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            ht = res.temp_preprocess_data_ht.ht()
            print_global_struct(ht.apply_models_globals)
            coverage_model, apply_log10_coverage = get_coverage_model(
                res,
                skip_coverage_model,
                log10_coverage,
                hl.eval(ht.apply_models_globals.high_cov_cutoff),
            )
            ht = create_per_variant_expected_ht(
                ht,
                res.mutation_ht.ht().select("mu_snp"),
                res.model_plateau.he(),
                coverage_model=coverage_model,
                log10_coverage=apply_log10_coverage,
                use_mane_select=True,
            )
            ht.write(res.per_variant_apply_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info("Done computing per-variant expected variant count.")

        if args.compute_rr:
            logger.info("Computing rr per window and context...")
            res = resources.compute_rr
            res.check_resource_existence()

            ht = compute_rr_ht(
                res.genomic_features_ht.ht(),
                res.rr_models_ht.ht(),
                res.per_variant_apply_ht.ht(),
            )
            ht.write(res.regional_rr_ht.path, overwrite=overwrite)
            logger.info("Done computing rr.")

        if args.apply_models_aggregated:
            logger.info("Aggregating counts and applying models on aggregated data...")
            res = resources.apply_models_aggregated
            check_apply_models_resources(res, skip_coverage_model)

            hl._set_flags(use_new_shuffle="1")

            ht = res.temp_preprocess_data_ht.ht()
            print_global_struct(ht.apply_models_globals)
            coverage_model, apply_log10_coverage = get_coverage_model(
                res,
                skip_coverage_model,
                log10_coverage,
                hl.eval(ht.apply_models_globals.high_cov_cutoff),
            )
            ht = create_aggregated_expected_ht(
                ht,
                res.mutation_ht.ht().select("mu_snp"),
                res.model_plateau.he(),
                coverage_model=coverage_model,
                log10_coverage=apply_log10_coverage,
                use_mane_select=True,
                partition_hint=args.apply_models_aggregated_partition_hint,
            )
            ht.write(res.aggregated_expected_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info("Done with aggregated model application.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path(
                "expected_counts",
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

    parser.add_argument(
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

    parser.add_argument(
        "--apply-models-per-variant",
        help=(
            "Apply plateau and coverage models to variants in exome sites Table and"
            " context Table to compute expected variant counts per variant."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-rr",
        help=(
            "Compute the regional mutation rate (rr) per 1kb window and context from"
            " the expected_models.py rr models, normalized by its per-context"
            " expected-weighted mean over exome synonymous MANE Select sites in the"
            " --apply-models-per-variant output."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-models-aggregated",
        help=(
            "Sum the observed and possible counts by transcript, consequence,"
            " context, and coverage first, then apply the plateau and coverage"
            " models, without writing per-variant expected counts. This replaces"
            " --apply-models-per-variant and --aggregate-per-variant-expected; run"
            " --aggregate-by-constraint-groups --use-aggregated-expected afterward."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-models-aggregated-partition-hint",
        help=(
            "Target number of partitions for the group_by that sums the observed and"
            " possible counts in --apply-models-aggregated. Default is 100."
        ),
        type=int,
        default=100,
    )

    args = parser.parse_args()
    main(args)
