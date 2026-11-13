"""
Import the external data and prepare the per-site data every constraint step reads.

The steps, in order:
    - --import-dnms: import the published de novo mutations (DNMs) the regional
      mutation rate (rr) is modeled from, as a Table of DNM sites.
    - --import-genomic-features: import the Gnocchi 1kb regional genomic feature
      matrix rr is modeled from.
    - --prepare-context-ht: annotate the VEP context Table with methylation,
      coverage, allele number, frequency, filters, GERP, and genomic region.
    - --calculate-gerp-cutoffs: compute GERP percentile cutoffs for
      --preprocess-data in the same command.
    - --preprocess-data: add per-site observed and possible counts, mutation rate
      inputs, and model groups.
    - --compute-gene-quality-metrics: per-transcript coverage, mapping quality,
      segmental duplication, and low-complexity region metrics.

Then run expected_models.py.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.constraint import calculate_gerp_cutoffs
from gnomad.utils.reference_genome import get_reference_genome
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_MISSING,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import filter_for_test
from gnomad_constraint.utils.constraint import (
    compute_gene_quality_metrics,
    prepare_context_ht,
    prepare_ht_for_constraint_calculations,
)
from gnomad_constraint.utils.import_data import format_dnm_ht, format_genomic_features

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("prepare_data")
logger.setLevel(logging.INFO)


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


def main(args: argparse.Namespace) -> None:
    """Execute the prepare data steps of the constraint pipeline."""
    hl.init(
        log="/prepare_data.log",
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

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")
    check_gerp_cutoffs_step(args.calculate_gerp_cutoffs, args.preprocess_data)

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = constraint_res.get_constraint_resources(
        version,
        overwrite,
        test,
        directory_post_fix=directory_post_fix,
        path_post_fix=path_post_fix,
    )

    try:
        if args.import_dnms:
            logger.info("Importing the DNMs...")
            res = resources.import_dnms
            res.check_resource_existence()

            # The files are gzip (not bgzip) compressed and quote fields with commas.
            # TODO: The calls aren't QC-filtered. The ASD file has Confidence and Role
            #  columns and the DDD file has LCR and gnomad_non_neuro_AF columns that a
            #  filter could use.
            ht = hl.Table.union(
                *[
                    hl.import_table(path, delimiter=",", quote='"', force=True).select(
                        "Variant"
                    )
                    for path in [res.asd_dnm_path, res.ddd_dnm_path]
                ]
            )
            format_dnm_ht(ht).write(res.dnm_ht.path, overwrite=overwrite)
            logger.info("Done importing the DNMs.")

        if args.import_genomic_features:
            logger.info("Importing the regional genomic feature matrix...")
            res = resources.import_genomic_features
            res.check_resource_existence()

            ht = hl.import_table(
                res.published_features_path,
                types={c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
                missing=GENOMIC_FEATURE_MISSING,
                min_partitions=100,
            )
            ht = format_genomic_features(ht)
            ht = ht.checkpoint(new_temp_file("genomic_features", "ht"))
            # Check before writing, so a bad matrix never replaces a good Table.
            n_duplicated = ht.count() - ht.distinct().count()
            if n_duplicated:
                raise ValueError(f"{n_duplicated} windows are duplicated!")
            ht.write(res.genomic_features_ht.path, overwrite=overwrite)
            logger.info("Done importing the regional genomic feature matrix.")

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

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path(
                "prepare_data",
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
        "--import-dnms",
        help="Import the published DNMs as a Table of DNM sites.",
        action="store_true",
    )
    parser.add_argument(
        "--import-genomic-features",
        help="Import the Gnocchi 1kb regional genomic feature matrix.",
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

    args = parser.parse_args()
    main(args)
