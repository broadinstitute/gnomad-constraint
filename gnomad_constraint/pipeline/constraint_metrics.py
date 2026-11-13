"""
Compute the constraint metrics and prepare the release files.

The steps, in order, after expected_counts.py:
    - --aggregate-per-variant-expected: sum the per-variant expected counts by
      transcript and consequence.
    - --aggregate-by-constraint-groups: sum by constraint group (synonymous,
      missense, pLoF).
    - --compute-constraint-metrics: o/e, confidence intervals, z-scores, flags,
      pLI, ranks, and bins.
    - --prepare-release, --prepare-release-mutation-rate: format the constraint
      metrics and mutation rate Tables for release.
    - --export-release-tsv, --export-release-downsampling-tsv: export the release
      Table and its per-ancestry downsampling counts to TSV.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.constraint import explode_downsamplings_oe

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RELEASE_KEY_ORDER,
    VERSIONS,
)
from gnomad_constraint.utils.constraint import (
    aggregate_by_constraint_groups,
    aggregate_per_variant_expected_ht,
    compute_constraint_metrics,
    compute_constraint_percentile_bins,
    flatten_release_ht,
    lof_bin_thresholds_to_ht,
    prepare_release_ht,
    prepare_release_mutation_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_metrics")
logger.setLevel(logging.INFO)


def get_release_downsamplings(release_ht: hl.Table) -> hl.Struct:
    """
    Return the downsampling levels of the release Table.

    The release Table only has a ``downsamplings`` global when preprocessing used
    ``--include-downsamplings``.

    :param release_ht: Release constraint metrics Table.
    :return: Downsampling levels for each genetic ancestry group.
    """
    if "downsamplings" not in release_ht.globals:
        raise ValueError(
            "The release Table has no downsamplings global, so there is no"
            " downsampling TSV to export. Rerun --preprocess-data with"
            " --include-downsamplings, then the steps after it."
        )
    return hl.eval(release_ht.downsamplings)


def check_constraint_group_input(
    aggregate_per_variant_expected: bool,
    apply_models_aggregated: bool,
    use_aggregated_expected: bool,
) -> None:
    """
    Warn when ``--aggregate-by-constraint-groups`` reads an output this command didn't write.

    The step reads the ``--aggregate-per-variant-expected`` output, or with
    ``--use-aggregated-expected`` the ``--apply-models-aggregated`` output. If only
    the other one runs in the same command, the step groups a Table that may be from
    an earlier run instead of the one just written.

    :param aggregate_per_variant_expected: Whether
        ``--aggregate-per-variant-expected`` runs in this command.
    :param apply_models_aggregated: Whether ``--apply-models-aggregated`` runs in this
        command.
    :param use_aggregated_expected: Whether ``--use-aggregated-expected`` was passed.
    :return: None.
    """
    if apply_models_aggregated and not (
        use_aggregated_expected or aggregate_per_variant_expected
    ):
        logger.warning(
            "--apply-models-aggregated runs in this command, but without"
            " --use-aggregated-expected, --aggregate-by-constraint-groups reads the"
            " --aggregate-per-variant-expected output, which may be from an earlier"
            " run. Add --use-aggregated-expected to group the Table written here."
        )
    if (
        aggregate_per_variant_expected
        and use_aggregated_expected
        and not apply_models_aggregated
    ):
        logger.warning(
            "--aggregate-per-variant-expected runs in this command, but with"
            " --use-aggregated-expected, --aggregate-by-constraint-groups reads the"
            " --apply-models-aggregated output, which may be from an earlier run."
            " Drop --use-aggregated-expected to group the Table written here."
        )


def main(args: argparse.Namespace) -> None:
    """Execute the constraint metrics steps of the constraint pipeline."""
    hl.init(
        log="/constraint_metrics.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    directory_post_fix = args.directory_post_fix
    path_post_fix = args.path_post_fix
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = constraint_res.get_constraint_resources(
        version,
        overwrite,
        test,
        directory_post_fix=directory_post_fix,
        path_post_fix=path_post_fix,
        skip_pre_rank_metrics=args.skip_pre_rank_metrics,
        use_aggregated_expected=args.use_aggregated_expected,
    )

    try:
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

        if args.aggregate_by_constraint_groups:
            logger.info(
                "Aggregating observed and expected variant counts by constraint groups..."
            )
            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            check_constraint_group_input(
                aggregate_per_variant_expected=args.aggregate_per_variant_expected,
                # --apply-models-aggregated runs in expected_counts.py, never here.
                apply_models_aggregated=False,
                use_aggregated_expected=args.use_aggregated_expected,
            )
            res = resources.aggregate_by_constraint_groups
            res.check_resource_existence()
            if args.use_aggregated_expected:
                ht = res.aggregated_expected_ht.ht()
            else:
                ht = res.apply_ht.ht()

            aggregate_by_constraint_groups(
                ht,
                keys=tuple(k for k in ht.key if k in RELEASE_KEY_ORDER),
            ).write(res.constraint_group_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)
            logger.info("Done with aggregating by constraint groups.")

        if args.compute_constraint_metrics:
            logger.info(
                "Computing constraint metrics, including pLI scores, z scores, oe"
                " ratio, and confidence interval around oe ratio..."
            )
            res = resources.compute_constraint_metrics
            res.check_resource_existence()

            # Compute constraint metrics, excluding rank and bin annotations.
            if args.skip_pre_rank_metrics:
                logger.info(
                    "Skipping metrics computation, reusing %s.",
                    res.pre_rank_constraint_metrics_ht.path,
                )
            else:
                ht = res.constraint_group_ht.ht()
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
                ).write(res.pre_rank_constraint_metrics_ht.path, overwrite=overwrite)

            # Add rank and bin annotations as a separate phase so they can be
            # recomputed without rerunning the metrics above.
            logger.info("Adding rank and percentile bin annotations...")
            compute_constraint_percentile_bins(
                res.pre_rank_constraint_metrics_ht.ht(),
                use_mane_select_over_canonical=args.use_mane_select_over_canonical,
            ).write(res.constraint_metrics_ht.path, overwrite=overwrite)
            logger.info("Done with computing constraint metrics.")

        if args.prepare_release:
            logger.info("Preparing constraint metrics Table for release...")
            res = resources.prepare_release
            res.check_resource_existence()

            constraint_ht = res.constraint_metrics_ht.ht()

            release_ht = prepare_release_ht(
                constraint_ht,
                release_version=args.release_version or version,
            ).naive_coalesce(1000)
            release_ht.write(res.release_ht.path, overwrite=overwrite)
            logger.info("Done preparing release Table.")

        if args.prepare_release_mutation_rate:
            logger.info("Preparing mutation rate Table for release...")
            res = resources.prepare_release_mutation_rate
            res.check_resource_existence()

            mutation_ht = res.mutation_ht.ht()
            release_mutation_ht = prepare_release_mutation_ht(
                mutation_ht,
                release_version=args.release_version or version,
            )
            release_mutation_ht.write(res.release_mutation_ht.path, overwrite=overwrite)

            logger.info("Exporting release mutation rate TSV...")
            release_mutation_ht.export(res.release_mutation_tsv)
            logger.info("Done preparing and exporting release mutation rate Table.")

        if args.export_release_tsv:
            res = resources.export_release_tsv
            res.check_resource_existence()
            release_ht = res.release_ht.ht()

            logger.info("Exporting release TSV...")
            flatten_release_ht(release_ht).export(res.release_tsv)
            logger.info("Done exporting release TSV.")

            logger.info("Exporting LoF OE CI upper bin thresholds TSV...")
            lof_bin_thresholds_to_ht(release_ht).export(res.lof_threshold_tsv)
            logger.info("Done exporting LoF threshold TSV.")

        if args.export_release_downsampling_tsv:
            res = resources.export_release_downsampling_tsv
            res.check_resource_existence()
            release_ht = res.release_ht.ht()

            logger.info("Exporting release downsampling TSV...")
            downsampling_ht = explode_downsamplings_oe(
                release_ht,
                downsampling_meta=get_release_downsamplings(release_ht),
                metrics=["syn", "mis", "lof_hc_lc", "lof"],
            )
            downsampling_ht.export(res.release_downsampling_tsv)
            logger.info("Done exporting release downsampling TSV.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path(
                "constraint_metrics",
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
        "--skip-pre-rank-metrics",
        help=(
            "Skip computing the constraint metrics and reuse the existing pre-rank"
            " Table, recomputing only the rank and percentile bin annotations. Used to"
            " reissue a release with corrected ranks without rerunning the pipeline."
        ),
        action="store_true",
    )
    compute_constraint_args.add_argument(
        "--use-mane-select-over-canonical",
        help=(
            "Use MANE Select rather than canonical transcripts when determining which"
            " transcripts to rank, falling back to canonical for genes without a MANE"
            " Select transcript."
        ),
        action=argparse.BooleanOptionalAction,
        default=True,
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
            "Version string to set as the version global of the release Tables."
            " Default is --version."
        ),
        type=str,
        default=None,
    )
    prepare_release_args.add_argument(
        "--prepare-release-mutation-rate",
        help=(
            "Prepare the mutation rate Table for public release by selecting"
            " the scalar mutation rate (mu), trinucleotide-class flags, and"
            " restructuring globals. Also exports a TSV."
        ),
        action="store_true",
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
