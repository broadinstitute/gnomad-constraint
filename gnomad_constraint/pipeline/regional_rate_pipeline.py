"""Script to compute the regional mutation rate (rr) adjustment of expected counts.

Port of the Gnocchi regional adjustment (https://github.com/atgu/gnomad_nc_constraint)
retrained on synonymous de novo mutations (DNMs). For each strand-collapsed
trinucleotide context, a ridge-penalized logistic model separates synonymous DNM sites
from non-mutated synonymous sites using regional genomic features, and the rr of a 1kb
window is the model's predicted probability for the window relative to the average
window, normalized so that each context's expected-weighted mean rr over exome
synonymous sites is 1. The constraint pipeline multiplies per-variant expected counts by
rr after model calibration (constraint_pipeline.py --aggregate-per-variant-expected
--apply-regional-rr).

The genomic features cover chr1-22 only, so chrX and chrY get no rr (an rr of 1).

The constraint pipeline's annotated context Table and per-variant expected Table are
read-only inputs; the constraint models are not re-fit or re-applied.

The rr pipeline consists of the following steps:
    - --prepare-syn-dnms: combine the DNM cohorts into the synonymous DNM Table.
    - --import-genomic-features: import the 1kb regional genomic feature matrix.
    - --create-rr-training-set: create the DNM / non-mutated background training set,
      leaving out genes with a QC gene flag.
    - --train-rr-models: fit the per-context ridge-penalized rr models.
    - --compute-rr: compute the normalized rr per window and context.
"""

import argparse
import logging
from typing import Union

import hail as hl
from gnomad.resources.resource_utils import TableResource
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    GENOMIC_FEATURE_COLUMNS,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_RIDGE_PENALTY,
    RR_TRAINING_N_PARTITIONS,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import (
    filter_for_test,
    get_dnm_cohort_resources,
)
from gnomad_constraint.utils.regional_rate import (
    combine_syn_dnms,
    compute_rr_ht,
    compute_rr_normalization_factors,
    create_rr_training_ht,
    fit_rr_models,
    format_dnm_cohort,
    format_genomic_features,
    get_syn_mane_select_ht,
    get_syn_sites_ht,
    normalize_rr_ht,
    remove_qc_flagged_genes,
    rr_models_to_ht,
    rr_training_ht_to_pandas,
    summarize_genomic_features,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_pipeline")
logger.setLevel(logging.INFO)


def import_dnm_cohort(cohort: str, resource: Union[str, TableResource]) -> hl.Table:
    """
    Import one DNM cohort (see ``get_dnm_cohort_resources``).

    :param cohort: Cohort name.
    :param resource: Path (text cohorts) or TableResource (Hail Table cohorts).
    :return: Imported cohort Table.
    """
    if isinstance(resource, TableResource):
        return resource.ht()
    if cohort.startswith("satterstrom_2026"):
        # The files are gzip (not bgzip) compressed and quote fields with commas.
        return hl.import_table(resource, delimiter=",", quote='"', force=True)

    return hl.import_table(resource, impute=True)


def main(args: argparse.Namespace) -> None:
    """Compute the regional mutation rate (rr) adjustment."""
    hl.init(
        log="/regional_rate_pipeline.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    overwrite = args.overwrite

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_regional_rate_resources(
        version,
        overwrite,
        test,
        directory_post_fix=args.directory_post_fix,
        path_post_fix=args.path_post_fix,
    )

    try:
        if args.prepare_syn_dnms:
            logger.info("Combining DNM cohorts into the synonymous DNM Table...")
            res = resources.prepare_syn_dnms
            res.check_resource_existence()

            cohort_hts = {
                cohort: format_dnm_cohort(cohort, import_dnm_cohort(cohort, resource))
                for cohort, resource in get_dnm_cohort_resources().items()
            }
            context_ht = res.annotated_context_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht, use_gene_list=test_gene_list)

            ht, cohort_counts = combine_syn_dnms(cohort_hts, context_ht)
            ht = ht.checkpoint(res.syn_dnm_ht.path, overwrite=overwrite)

            # The synonymous DNM Table is small and checkpointed, so counting is
            # cheap.
            n_dnms = ht.count()
            logger.info(
                "Synonymous DNMs per cohort (before deduplication): %s", cohort_counts
            )
            logger.info(
                "%d synonymous DNMs after deduplication (%d duplicates removed).",
                n_dnms,
                sum(cohort_counts.values()) - n_dnms,
            )
            logger.info("Done creating the synonymous DNM Table.")

        if args.import_genomic_features:
            logger.info("Importing the regional genomic feature matrix...")
            res = resources.import_genomic_features
            res.check_resource_existence()

            path = res.published_features
            ht = hl.import_table(
                path,
                types={c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
                missing=["NA", "nan", "NaN", ""],
                min_partitions=100,
                force=path.endswith(".gz"),
            )
            ht = format_genomic_features(ht)
            ht = ht.checkpoint(res.genomic_features_ht.path, overwrite=overwrite)

            summary = summarize_genomic_features(ht)
            if summary["n_duplicated_windows"]:
                raise ValueError(
                    f"{summary['n_duplicated_windows']} windows are duplicated in the"
                    " imported feature matrix!"
                )
            for contig, s in sorted(summary["contigs"].items()):
                logger.info(
                    "    %s: %d windows; features missing in all windows: %s",
                    contig,
                    s["n_windows"],
                    s["missing_features"],
                )
            logger.info("Done importing the regional genomic feature matrix.")

        if args.create_rr_training_set:
            logger.info("Creating the rr training set...")
            res = resources.create_rr_training_set
            res.check_resource_existence()

            context_ht = res.annotated_context_ht.ht()
            exomes_ht = res.exomes_sites_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht, use_gene_list=test_gene_list)
                exomes_ht = filter_for_test(exomes_ht, use_gene_list=test_gene_list)
            syn_sites_ht = get_syn_sites_ht(context_ht).checkpoint(
                new_temp_file("rr_syn_sites", "ht")
            )
            ht = create_rr_training_ht(
                syn_sites_ht,
                res.syn_dnm_ht.ht(),
                exomes_ht,
                res.genomic_features_ht.ht(),
                RR_BACKGROUND_FRACTION,
                RR_BACKGROUND_SEED,
            )
            ht = remove_qc_flagged_genes(ht, res.gene_flags_ht.ht())
            ht = ht.checkpoint(new_temp_file("rr_training", "ht"))
            # Sampled sites keep the context Table partitioning, so coalesce the
            # small training Table before writing it.
            ht = ht.naive_coalesce(RR_TRAINING_N_PARTITIONS).checkpoint(
                res.rr_training_ht.path, overwrite=overwrite
            )
            logger.info(
                "rr training set sites (1 = DNM, 0 = background): %s",
                dict(ht.aggregate(hl.agg.counter(ht.y))),
            )
            logger.info("Done creating the rr training set.")

        if args.train_rr_models:
            logger.info(
                "Training the rr models (ridge penalty %s)...", RR_RIDGE_PENALTY
            )
            res = resources.train_rr_models
            res.check_resource_existence()

            df = rr_training_ht_to_pandas(res.rr_training_ht.ht())
            models = fit_rr_models(df)
            rr_models_to_ht(models).write(res.rr_models_ht.path, overwrite=overwrite)
            logger.info("Done training the rr models.")

        if args.compute_rr:
            logger.info("Computing rr per window and context...")
            res = resources.compute_rr
            res.check_resource_existence()

            models_ht = res.rr_models_ht.ht()
            rr_ht = compute_rr_ht(res.genomic_features_ht.ht(), models_ht.collect())
            rr_ht = rr_ht.checkpoint(new_temp_file("rr_unnormalized", "ht"))

            per_variant_ht = res.per_variant_expected_ht.ht()
            if test:
                per_variant_ht = filter_for_test(
                    per_variant_ht, use_gene_list=test_gene_list
                )
            syn_ht = get_syn_mane_select_ht(per_variant_ht).checkpoint(
                new_temp_file("rr_syn_mane_select", "ht")
            )
            factors = compute_rr_normalization_factors(syn_ht, rr_ht)
            logger.info("rr normalization factors: %s", factors)

            rr_ht = normalize_rr_ht(rr_ht, factors)
            rr_ht = rr_ht.annotate_globals(
                rr_globals=hl.struct(
                    ridge_penalty=RR_RIDGE_PENALTY,
                    normalization_factors=hl.literal(
                        factors, hl.tdict(hl.tstr, hl.tfloat64)
                    ),
                    models_ht_path=res.rr_models_ht.path,
                    features_ht_path=res.genomic_features_ht.path,
                )
            )
            rr_ht = rr_ht.checkpoint(res.regional_rr_ht.path, overwrite=overwrite)

            # Log rr quantiles per context from the written Table.
            stats = rr_ht.aggregate(
                hl.agg.explode(
                    lambda kv: hl.agg.group_by(
                        kv[0],
                        hl.agg.approx_quantiles(kv[1], [0.0, 0.05, 0.5, 0.95, 1.0]),
                    ),
                    rr_ht.adj_r.items(),
                )
            )
            logger.info(
                "rr quantiles (min, 5%%, median, 95%%, max) per context:\n%s",
                "\n".join(f"    {c}: {q}" for c, q in sorted(stats.items())),
            )
            logger.info("Done computing rr.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path(
                "regional_rate_pipeline",
                version=version,
                directory_post_fix=args.directory_post_fix,
                path_post_fix=args.path_post_fix,
            )
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
        help="Post-fix to append to output file names.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to filter the context, exomes sites, and per-variant expected"
            " Tables to only chromosome 20, chromosome X, and chromosome Y for"
            " testing."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-gene-list",
        help=(
            "Whether to filter the context, exomes sites, and per-variant expected"
            " Tables to only six genes for testing (see --test-gene-list in"
            " constraint_pipeline.py)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--prepare-syn-dnms",
        help="Combine the DNM cohorts into the synonymous DNM Table.",
        action="store_true",
    )
    parser.add_argument(
        "--import-genomic-features",
        help="Import the 1kb regional genomic feature matrix.",
        action="store_true",
    )
    parser.add_argument(
        "--create-rr-training-set",
        help=(
            "Create the rr training set of synonymous DNM sites and a sampled"
            " background of non-mutated synonymous sites, leaving out genes with a QC"
            " gene flag."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--train-rr-models",
        help="Fit the per-context ridge-penalized rr models.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-rr",
        help=(
            "Compute rr per 1kb window and context, normalized by its per-context"
            " expected-weighted mean over exome synonymous MANE Select sites."
        ),
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
