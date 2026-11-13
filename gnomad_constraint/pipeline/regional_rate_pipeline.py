"""
Script to compute the regional mutation rate (rr) of each 1kb window and context.

For each strand-collapsed trinucleotide context, a ridge-penalized logistic model
separates synonymous de novo mutation (DNM) sites from non-mutated synonymous sites
using regional genomic features. A window's rr is the model's prediction for the
window relative to the average window, normalized so that each context's mean rr,
weighted by the adj expected count, over exome synonymous Ensembl MANE Select sites
is 1.
``constraint_metrics.py --aggregate-per-variant-expected --apply-regional-rr``
applies it.

The DNM and genomic feature Tables are written by ``prepare_data.py``. The
annotated context Table is read with the same --version, --test, and
--directory-post-fix as the rr outputs, and the per-variant expected Table also with
the same --path-post-fix: the values ``--apply-regional-rr`` must be run with, so rr
is normalized on the Table it is applied to.

The rr pipeline consists of the following steps:
    - --create-rr-site-features: create the DNM and non-mutated background sites
      with their window's features, leaving out genes with a QC gene flag.
    - --create-rr-training-set: prepare the sites for the models: each site's
      label, context, chromosome fold, and its context's features.
    - --tune-rr-penalty: choose the ridge penalty by the held-out log loss of
      chromosome-split models.
    - --train-rr-models: fit the per-context rr models.
    - --compute-rr: compute the normalized rr per window and context.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.constraint import (
    choose_regional_mutation_rate_penalty,
    fit_regional_mutation_rate_models,
    prepare_regional_mutation_rate_training_ht,
    tune_regional_mutation_rate_penalty,
)
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
from gnomad_constraint.utils.regional_rate import (
    compute_rr_ht,
    create_rr_site_features_ht,
    get_rr_fold_expr,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_pipeline")
logger.setLevel(logging.INFO)


def main(args: argparse.Namespace) -> None:
    """Compute the regional mutation rate (rr)."""
    hl.init(
        log="/regional_rate_pipeline.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test = args.test
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

        if args.compute_rr:
            logger.info("Computing rr per window and context...")
            res = resources.compute_rr
            res.check_resource_existence()

            ht = compute_rr_ht(
                res.genomic_features_ht.ht(),
                res.rr_models_ht.ht(),
                res.per_variant_expected_ht.ht(),
            )
            ht.write(res.regional_rr_ht.path, overwrite=overwrite)
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
        help=(
            "Post-fix to append to the output directory path. Also applies to the"
            " per-variant expected input, as in the constraint pipeline scripts."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--path-post-fix",
        help=(
            "Post-fix to append to output file names. Also applies to the"
            " per-variant expected input, as in the constraint pipeline scripts."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to run on the test data: reads the --test outputs of"
            " prepare_data.py and expected_counts.py, and writes the outputs to the"
            " testing bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--create-rr-site-features",
        help=(
            "Create the rr site features: synonymous DNM sites and a sampled"
            " background of non-mutated synonymous sites with their window's"
            " features, leaving out genes with a QC gene flag."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--create-rr-training-set",
        help=(
            "Prepare the rr site features for the models: each site's label, context,"
            " chromosome fold, and its context's features in RR_CONTEXT_FEATURES."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-qc-flagged-gene-filter",
        help=(
            "With --create-rr-site-features, keep the sites in genes with a QC gene"
            " flag in the v4.1.2 constraint metrics, which are left out by default."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--tune-rr-penalty",
        help=(
            "Choose the ridge penalty of the rr models: the penalty in --rr-penalties"
            " with the lowest held-out log loss of models fit on odd chromosomes and"
            " scored on even ones, and vice versa, summed over contexts."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--rr-penalties",
        help=(
            "Per-site ridge penalties tried by --tune-rr-penalty. Default is"
            f" {RR_RIDGE_PENALTIES}."
        ),
        type=float,
        nargs="+",
        default=RR_RIDGE_PENALTIES,
    )
    parser.add_argument(
        "--train-rr-models",
        help=(
            "Fit the per-context rr models with the penalty chosen by"
            " --tune-rr-penalty."
        ),
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

    args = parser.parse_args()
    main(args)
