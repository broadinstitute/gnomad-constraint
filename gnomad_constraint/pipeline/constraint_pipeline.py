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
from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import all_sites_an
from gnomad.utils.constraint import (
    annotate_with_mu,
    assemble_constraint_context_ht,
    build_models,
    explode_downsamplings_oe,
)
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import update_loftee_end_trunc_filter
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.utils.constraint import (
    aggregate_by_constraint_groups,
    aggregate_per_variant_expected_ht,
    calculate_gerp_cutoffs,
    calculate_mu_by_downsampling,
    compute_constraint_metrics,
    create_per_variant_expected_ht,
    create_training_set,
    prepare_ht_for_constraint_calculations,
    print_global_struct,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def filter_for_test(
    ht: hl.Table,
    use_gene_list: bool = False,
) -> hl.Table:
    """
    Filter `ht` to chr20, chrX, and chrY or a gene list for testing.

    :param ht: Table to filter.
    :param use_gene_list: Whether to use a gene list for testing instead of all of
        chr20, chrX, and chrY for testing.
    :return: Filtered Table for testing.
    """
    rg = get_reference_genome(ht.locus)
    if use_gene_list:
        if rg == "GRCh37":
            keep_regions = [
                "20:49505585-49547958",  # ADNP
                "20:853296-896977",  # ANGPT4
                "X:13752832-13787480",  # OFD1
                "X:57313139-57515629",  # FAAH2
                "Y:2803112-2850547",  # ZFY
            ]
        else:
            keep_regions = [
                "chr20:50888916-50931437",  # ADNP
                "chr20:869900-916334",  # ANGPT4
                "chrX:13734743-13777955",  # OFD1
                "chrX:57286706-57489193",  # FAAH2
                "chrY:2935281-2982506",  # ZFY
            ]
        keep = [hl.parse_locus_interval(c, reference_genome=rg) for c in keep_regions]
    else:
        keep = [
            hl.parse_locus_interval(c, reference_genome=rg)
            for c in [rg.contigs[19], rg.x_contigs[0], rg.y_contigs[0]]
        ]
        logger.info("Filtering the context HT to chr20, chrX, and chrY for testing...")

    ht = hl.filter_intervals(ht, keep)

    return ht


def run_prepare_context(
    resources: PipelineResourceCollection,
    test: bool = False,
    test_gene_list: bool = False,
) -> hl.Table:
    """
    Annotate the context Table with coverage, AN, and frequency annotations.

    Uses `assemble_constraint_context_ht` to annotate the context Table with annotations
    that are used in downstream steps of the constraint pipeline.

    :param resources: PipelineResourceCollection containing resources for the constraint
        pipeline.
    :param test: Whether to filter the context Table to only chr20, chrX, and chrY for
        testing.
    :param test_gene_list: Whether to filter the context Table to a gene list for
        testing.
    :return: Annotated context Table.
    """
    # We use naive_coalesce on the context Table because it has a large number of
    # partitions which caused some issues with Hail 0.2.133. 5000 partitions was a
    # number that worked well for the context Table in the past.
    ht = resources.context_ht.ht().naive_coalesce(5000)

    if test:
        ht = filter_for_test(ht, use_gene_list=test_gene_list)

    def _build_ht_dict(ht_name: str, keep: List[str] = None):
        dts = ["exomes", "genomes"]
        hts = {d: getattr(resources, f"{d}_{ht_name}_ht").ht() for d in dts}
        return {d: t.select(*keep) for d, t in hts.items()} if keep else hts

    # There was a bug in the GERP cutoffs used to filter transcripts with the
    # "END_TRUNC" filter in the LOFTEE VEP plugin resulting in some transcripts
    # being considered "HC" when they should have been "LC". We use the
    # `update_loftee_end_trunc_filter` function to correct this issue.
    ht = ht.annotate(
        vep=ht.vep.annotate(
            transcript_consequences=update_loftee_end_trunc_filter(
                ht.vep.transcript_consequences
            )
        )
    )
    ht = assemble_constraint_context_ht(
        ht,
        coverage_hts=_build_ht_dict("coverage"),
        an_hts=_build_ht_dict("an"),
        freq_hts=_build_ht_dict("sites", ["freq"]),
        filter_hts=_build_ht_dict("sites", ["filters"]),
        methylation_ht=resources.methylation_ht.ht(),
        gerp_ht=constraint_res.get_gerp_ht(get_reference_genome(ht.locus).name),
        transformation_funcs=None,
    )

    # Add annotation for exome coverage and genomic region (autosome/PAR, X non-PAR,
    # Y non-PAR).
    genomic_region_expr = (
        hl.case()
        .when(ht.locus.in_autosome_or_par(), "autosome_or_par")
        .when(ht.locus.in_x_nonpar(), "chrx_nonpar")
        .when(ht.locus.in_y_nonpar(), "chry_nonpar")
        .or_missing()
    )

    # Add annotation for SFS bin.
    sfs_bin_cutoffs = [0, 1e-6, 2e-6, 4e-6, 2e-5, 5e-5, 5e-4, 5e-3, 0.5]
    af_expr = ht.freq.exomes[0].AF
    sfs_bin_expr = hl.case().when(hl.is_missing(af_expr), 0)
    for i, af in enumerate(sfs_bin_cutoffs):
        sfs_bin_expr = sfs_bin_expr.when(af_expr <= af, i)

    sfs_bin_expr = sfs_bin_expr.or_missing()

    adj_r_ht = hl.read_table(
        "gs://gnomad/v4.1/constraint/resources/annotations/ht/adj_r_per_context_methyl_genome_1kb_autosome.agg.ht"
    )

    ht = ht.annotate(
        coverage=hl.struct(
            exomes=ht.coverage.exomes.select("mean", "median_approx"),
            genomes=ht.coverage.genomes.select("mean", "median_approx"),
        ),
        AN=hl.struct(
            exomes=ht.AN.exomes[0],
            genomes=ht.AN.genomes[0],
        ),
        genomic_region=genomic_region_expr,
        adj_r=adj_r_ht[ht.locus].adj_r[ht.context],
        sfs_bin=sfs_bin_expr,
    )

    return ht


def get_constraint_resources(
    version: str,
    custom_vep_annotation: str,
    overwrite: bool,
    test: bool,
    models: List[str] = ["plateau", "coverage"],
    directory_post_fix: Optional[str] = None,
    path_post_fix: Optional[str] = None,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the constraint pipeline.

    :param version: Version of constraint resources to use.
    :param custom_vep_annotation: Custom VEP annotation to use for applying models
        resources.
    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :param models: List of models to use. Default is ["plateau", "coverage"].
    :param directory_post_fix: Post-fix to add to the directory path of the resources.
    :param path_post_fix: Post-fix to add to the path of the resources.
    :return: PipelineResourceCollection containing resources for all steps of the
        constraint pipeline.
    """
    # Initialize constraint pipeline resource collection.
    constraint_pipeline = PipelineResourceCollection(
        pipeline_name="constraint",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the constraint pipeline.
    context_res = constraint_res.get_vep_context_ht(version)
    context_build = get_reference_genome(context_res.ht().locus).name

    # Make dictionary for prepare_context input Tables.
    input_hts = {
        "context_ht": context_res,
        "methylation_ht": constraint_res.get_methylation_ht(context_build),
    }
    for d in ["exomes", "genomes"]:
        input_hts[f"{d}_coverage_ht"] = constraint_res.get_coverage_ht(d, version)
        input_hts[f"{d}_sites_ht"] = constraint_res.get_sites_resource(d, version)
        input_hts[f"{d}_an_ht"] = all_sites_an(d)

    common_params = {
        "version": version,
        "test": test,
        "directory_post_fix": directory_post_fix,
    }

    prepare_context = PipelineStepResourceCollection(
        "--prepare-context-ht",
        output_resources={
            "annotated_context_ht": constraint_res.get_annotated_context_ht(
                **common_params
            )
        },
        input_resources={"gnomAD resources": input_hts},
    )
    preprocess_data = PipelineStepResourceCollection(
        "preprocess data for downstream steps",
        output_resources={
            "temp_preprocess_data_ht": constraint_res.get_preprocessed_ht(
                **common_params
            ),
        },
        pipeline_input_steps=[prepare_context],
    )
    calculate_gerp_cutoffs = PipelineStepResourceCollection(
        "--calculate-gerp-cutoffs",
        output_resources={},
        pipeline_input_steps=[prepare_context],
    )
    calculate_mutation_rate = PipelineStepResourceCollection(
        "--calculate-mutation-rate",
        output_resources={
            "mutation_ht": constraint_res.get_mutation_ht(**common_params)
        },
        pipeline_input_steps=[preprocess_data],
    )
    create_training_set = PipelineStepResourceCollection(
        "--create-training-set",
        output_resources={
            f"train_ht": constraint_res.get_training_dataset(
                **common_params, path_post_fix=path_post_fix
            ),
            f"train_tsv": constraint_res.get_training_tsv_path(
                **common_params, path_post_fix=path_post_fix
            ),
        },
        pipeline_input_steps=[preprocess_data, calculate_mutation_rate],
    )
    build_models = PipelineStepResourceCollection(
        "--build-models",
        output_resources={
            f"model_{m}": constraint_res.get_models(
                m, **common_params, path_post_fix=path_post_fix
            )
            for m in models
        },
        pipeline_input_steps=[create_training_set],
    )
    apply_models_per_variant = PipelineStepResourceCollection(
        "--apply-models-per-variant",
        output_resources={
            "per_variant_apply_ht": constraint_res.get_per_variant_expected_dataset(
                custom_vep_annotation, **common_params, path_post_fix=path_post_fix
            )
        },
        pipeline_input_steps=[preprocess_data, calculate_mutation_rate, build_models],
    )
    aggregate_per_variant_expected = PipelineStepResourceCollection(
        "--aggregate-per-variant-expected",
        output_resources={
            f"apply_ht": constraint_res.get_aggregated_per_variant_expected(
                custom_vep_annotation, **common_params, path_post_fix=path_post_fix
            )
        },
        pipeline_input_steps=[
            apply_models_per_variant,
            calculate_mutation_rate,
            build_models,
        ],
    )
    aggregate_by_constraint_groups = PipelineStepResourceCollection(
        "--aggregate-by-constraint-groups",
        output_resources={
            f"constraint_group_ht": constraint_res.get_constraint_group_ht(
                custom_vep_annotation, **common_params, path_post_fix=path_post_fix
            )
        },
        pipeline_input_steps=[aggregate_per_variant_expected],
    )
    compute_constraint_metrics = PipelineStepResourceCollection(
        "--compute-constraint-metrics",
        output_resources={
            "constraint_metrics_ht": constraint_res.get_constraint_metrics_dataset(
                custom_vep_annotation, **common_params, path_post_fix=path_post_fix
            )
        },
        pipeline_input_steps=[aggregate_by_constraint_groups],
    )
    export_tsv = PipelineStepResourceCollection(
        "--export-tsv",
        output_resources={
            "constraint_metrics_tsv": constraint_res.get_constraint_tsv_path(
                **common_params, path_post_fix=path_post_fix
            ),
            "downsampling_constraint_metrics_tsv": (
                constraint_res.get_downsampling_constraint_tsv_path(**common_params)
            ),
        },
        pipeline_input_steps=[compute_constraint_metrics],
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
            "apply_models_per_variant": apply_models_per_variant,
            "aggregate_per_variant_expected": aggregate_per_variant_expected,
            "aggregate_by_constraint_groups": aggregate_by_constraint_groups,
            "compute_constraint_metrics": compute_constraint_metrics,
            "export_tsv": export_tsv,
        }
    )

    return constraint_pipeline


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

    if version not in constraint_res.VERSIONS:
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
    resources = get_constraint_resources(
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
            ht = run_prepare_context(res, test=test, test_gene_list=test_gene_list)
            ht.write(res.annotated_context_ht.path, overwrite)

            logger.info("Done annotating the VEP context Table.")

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
            # hl._set_flags(use_new_shuffle="1")

            ht = res.per_variant_apply_ht.ht()
            # ht = res.per_variant_apply_ht.ht(read_args={"_n_partitions": 8000})
            ht = aggregate_per_variant_expected_ht(
                ht, include_mu_annotations_in_grouping=True
            )
            ht = ht.checkpoint(
                "gs://gnomad/v4.1/constraint_coverage_corrected/apply_models/transcript_consequences/gnomad.v4.1.per_variant_expected.aggregated_with_mu_annotations.coverage_corrected.with_downsamplings.ht",
                overwrite=overwrite,
            )
            # hl._set_flags(use_new_shuffle=None)

            ht = aggregate_per_variant_expected_ht(ht)
            ht.write(res.apply_ht.path, overwrite=overwrite)

            logger.info(
                "Done aggregating per-variant expected variant count by transcript, "
                "consequence annotations, and consequence modifier annotations."
            )

        if args.aggregate_by_constraint_groups:
            logger.info(
                "Aggregating observed and expected variant counts by constraint groups..."
            )
            res = resources.aggregate_by_constraint_groups
            res.check_resource_existence()

            # Use new shuffle method to prevent shuffle errors.
            hl._set_flags(use_new_shuffle="1")
            ht = hl.read_table(
                "gs://gnomad/v4.1/constraint_coverage_corrected/apply_models/transcript_consequences/gnomad.v4.1.per_variant_expected.aggregated_with_mu_annotations.coverage_corrected.with_downsamplings.ht"
            )
            aggregate_by_constraint_groups(
                ht,
                keys=tuple(
                    [
                        i
                        for i in list(ht.key)
                        if i
                        in [
                            "gene",
                            "transcript",
                            "canonical",
                            "mane_select",
                            "gene_id",
                            "context",
                            "ref",
                            "alt",
                            "methylation_level",
                        ]
                    ]
                ),
            ).write(
                "gs://gnomad/v4.1/constraint_coverage_corrected/apply_models/transcript_consequences/gnomad.v4.1.constraint_group_with_mu_annotations.coverage_corrected.with_downsamplings.ht",
                overwrite=overwrite,
            )

            ht = res.apply_ht.ht()
            aggregate_by_constraint_groups(
                ht,
                keys=tuple(
                    [
                        i
                        for i in list(ht.key)
                        if i
                        in ["gene", "transcript", "canonical", "mane_select", "gene_id"]
                    ]
                ),
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

            # Compute constraint metrics.
            ht = res.constraint_group_ht.ht(read_args={"_n_partitions": 10000})
            compute_constraint_metrics(
                ht=ht,
                gencode_ht=constraint_res.get_gencode_ht(version),
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
                # ).select_globals(
                #   "version", "apply_model_params", "constraint_meta", "sd_raw_z"
            ).write(res.constraint_metrics_ht.path, overwrite=overwrite)
            logger.info("Done with computing constraint metrics.")

        if args.export_tsv:
            res = resources.export_tsv
            res.check_resource_existence()
            logger.info("Exporting constraint tsv...")

            ht = res.constraint_metrics_ht.ht()
            # If downsamplings per genetic ancestry group are present, export
            # downsamplings to a separate tsv and drop from the main metrics tsv.
            if args.genetic_ancestry_groups:
                downsampling_ht = explode_downsamplings_oe(
                    ht,
                    downsampling_meta=hl.eval(ht.apply_model_params.downsampling_meta),
                )

                # Drop downsampling annotations from the main metrics Table.
                ht = ht.annotate(
                    **{
                        i: ht[i].drop(*["gen_anc_exp", "gen_anc_obs"])
                        for i in ["lof_hc_lc", "lof", "syn", "mis"]
                    }
                )
                # Export separate downsampling Table.
                downsampling_ht.export(res.downsampling_constraint_metrics_tsv)
            ht = ht.flatten()
            ht.export(res.constraint_metrics_tsv)

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
            f" {constraint_res.CURRENT_VERSION}."
        ),
        type=str,
        default=constraint_res.CURRENT_VERSION,
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
        choices=constraint_res.CUSTOM_VEP_ANNOTATIONS,
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
    compute_constraint_args.add_argument(
        "--export-tsv",
        help="Export constraint metrics to tsv file.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
