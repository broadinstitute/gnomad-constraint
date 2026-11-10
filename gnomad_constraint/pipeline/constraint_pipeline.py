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
    - optionally, apply the regional mutation rate (rr) adjustment to the per-variant
      expected counts (--aggregate-per-variant-expected --apply-regional-rr; rr is
      computed by regional_rate_pipeline.py)
    - calculate metrics
"""

import argparse
import logging
import math
from typing import Optional, Tuple

import hail as hl
from gnomad.utils.constraint import (
    build_models,
    calculate_gerp_cutoffs,
    explode_downsamplings_oe,
)
from gnomad.utils.file_utils import file_exists, print_global_struct
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_qc.resource_utils import PipelineStepResourceCollection

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RELEASE_KEY_ORDER,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import filter_for_test
from gnomad_constraint.utils.constraint import (
    aggregate_by_constraint_groups,
    aggregate_per_variant_expected_ht,
    calculate_mu_by_downsampling,
    compute_constraint_metrics,
    compute_constraint_percentile_bins,
    compute_gene_quality_metrics,
    create_aggregated_expected_ht,
    create_per_variant_expected_ht,
    create_training_set,
    flatten_release_ht,
    lof_bin_thresholds_to_ht,
    prepare_context_ht,
    prepare_ht_for_constraint_calculations,
    prepare_release_ht,
    prepare_release_mutation_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
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


def check_gerp_cutoffs_step(
    calculate_gerp_cutoffs: bool, preprocess_data: bool
) -> None:
    """
    Check that ``--calculate-gerp-cutoffs`` runs with ``--preprocess-data``.

    The cutoffs are not saved, and only ``--preprocess-data`` in the same command
    uses them.

    :param calculate_gerp_cutoffs: Whether ``--calculate-gerp-cutoffs`` runs in this
        command.
    :param preprocess_data: Whether ``--preprocess-data`` runs in this command.
    :return: None.
    """
    if calculate_gerp_cutoffs and not preprocess_data:
        raise ValueError(
            "--calculate-gerp-cutoffs doesn't save its cutoffs, and only"
            " --preprocess-data in the same command uses them. Add --preprocess-data."
        )


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
    skip_coverage_model = args.skip_coverage_model
    gerp_lower_cutoff = args.calculate_mutation_rate_gerp_lower_cutoff
    gerp_upper_cutoff = args.calculate_mutation_rate_gerp_upper_cutoff
    log10_coverage = args.use_logarithmic_coverage_model

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")
    check_gerp_cutoffs_step(args.calculate_gerp_cutoffs, args.preprocess_data)

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
        skip_pre_rank_metrics=args.skip_pre_rank_metrics,
        use_aggregated_expected=args.use_aggregated_expected,
        apply_regional_rr=args.apply_regional_rr,
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
                adj_r_ht=res.adj_r_ht.ht(),
                syn_adj_r_ht=res.syn_adj_r_ht.ht(),
            )
            ht.write(res.annotated_context_ht.path, overwrite)

            logger.info("Done annotating the VEP context Table.")

        if args.calculate_gerp_cutoffs:
            logger.warning(
                "Calculating new GERP cutoffs to be used instead of"
                " '--calculate-mutation-rate-gerp-lower-cutoff' and"
                " '--calculate-mutation-rate-gerp-upper-cutoff'."
            )
            res = resources.calculate_gerp_cutoffs
            res.check_resource_existence()
            ht = res.annotated_context_ht.ht()
            gerp_lower_cutoff, gerp_upper_cutoff = calculate_gerp_cutoffs(
                ht.filter(ht.genomic_region == "autosome_or_par")
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
                calculate_mutation_rate_gerp_lower_cutoff=gerp_lower_cutoff,
                calculate_mutation_rate_gerp_upper_cutoff=gerp_upper_cutoff,
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
                res.seg_dup_intervals.ht(),
                res.lcr_intervals.ht(),
            )
            gene_quality_ht.write(res.gene_quality_metrics_ht.path, overwrite=overwrite)
            logger.info("Done computing gene quality metrics.")

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
            ht = aggregate_per_variant_expected_ht(
                ht,
                regional_rr_ht=(
                    res.regional_rr_ht.ht() if args.apply_regional_rr else None
                ),
            )
            ht.write(res.apply_ht.path, overwrite=overwrite)
            hl._set_flags(use_new_shuffle=None)

            logger.info(
                "Done aggregating per-variant expected variant count by transcript, "
                "consequence annotations, and consequence modifier annotations."
            )

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

        if args.aggregate_by_constraint_groups:
            logger.info(
                "Aggregating observed and expected variant counts by constraint groups..."
            )
            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")

            check_constraint_group_input(
                aggregate_per_variant_expected=args.aggregate_per_variant_expected,
                apply_models_aggregated=args.apply_models_aggregated,
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
                "constraint_pipeline",
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
            " of GERP scores at autosome/PAR sites in the context Table. The cutoffs"
            " are not saved, so this must run in the same command as"
            " --preprocess-data, which uses them instead of"
            " --calculate-mutation-rate-gerp-lower-cutoff and"
            " --calculate-mutation-rate-gerp-upper-cutoff."
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
            " Table to define the 5th percentile)."
        ),
        type=float,
        default=-3.9885,
    )
    preprocess_args.add_argument(
        "--calculate-mutation-rate-gerp-upper-cutoff",
        help=(
            "Maximum GERP score for variant to be included when calculating the"
            " mutation rate. Default is 2.6607 (precalculated on the GRCh37 context"
            " Table to define the 95th percentile)."
        ),
        type=float,
        default=2.6607,
    )
    preprocess_args.add_argument(
        "--exome-coverage-metric",
        help=(
            "Metric to use for exome coverage: 'median' (median depth), 'AN' (allele"
            " number), or 'AN_percent' (percent of samples with a genotype). Default is"
            " 'AN_percent'."
        ),
        type=str,
        default="AN_percent",
        choices=["median", "AN", "AN_percent"],
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
        "--apply-regional-rr",
        help=(
            "Multiply the per-variant expected counts by the regional mutation rate"
            " (rr) of the variant's 1kb window and context before aggregating, after"
            " model calibration. Reads the rr Table written by"
            " regional_rate_pipeline.py --compute-rr for the same version, --test, and"
            " --directory-post-fix."
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
