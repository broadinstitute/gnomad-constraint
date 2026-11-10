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
read-only inputs; the constraint models are not re-fit or re-applied. The per-variant
expected Table is read with the same --version, --test, --directory-post-fix, and
--path-post-fix as the rr outputs, which are the values constraint_pipeline.py
--apply-regional-rr must be run with, so rr is normalized on the Table it is applied
to. The DNM Table and genomic feature Table are written by ``import_data.py
--import-dnms`` and ``--import-genomic-features``.

The rr pipeline consists of the following steps:
    - --create-rr-training-set: create the DNM / non-mutated background training set,
      leaving out genes with a QC gene flag.
    - --tune-rr-penalty: choose the ridge penalty by the held-out log loss of
      chromosome-split models.
    - --train-rr-models: fit the per-context ridge-penalized rr models.
    - --compute-rr: compute the normalized rr per window and context.
"""

import argparse
import logging

import hail as hl
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_RIDGE_PENALTIES,
    RR_TRAINING_N_PARTITIONS,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import filter_for_test
from gnomad_constraint.utils.regional_rate import (
    choose_rr_penalty,
    compute_rr_ht,
    compute_rr_normalization_factors,
    create_rr_training_ht,
    fit_rr_models,
    get_syn_mane_select_ht,
    get_syn_sites_ht,
    normalize_rr_ht,
    remove_qc_flagged_genes,
    rr_models_to_ht,
    rr_training_ht_to_pandas,
    tune_rr_penalty,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_pipeline")
logger.setLevel(logging.INFO)


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
        if args.create_rr_training_set:
            logger.info("Creating the rr training set...")
            res = resources.create_rr_training_set
            res.check_resource_existence()

            context_ht = res.annotated_context_ht.ht()
            exomes_ht = res.exomes_sites_ht.ht()
            dnm_ht = res.dnm_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht, use_gene_list=test_gene_list)
                exomes_ht = filter_for_test(exomes_ht, use_gene_list=test_gene_list)
                dnm_ht = filter_for_test(dnm_ht, use_gene_list=test_gene_list)
            syn_sites_ht = get_syn_sites_ht(context_ht).checkpoint(
                new_temp_file("rr_syn_sites", "ht")
            )
            ht = create_rr_training_ht(
                syn_sites_ht,
                dnm_ht,
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

        if args.tune_rr_penalty:
            logger.info("Tuning the rr ridge penalty over %s...", args.rr_penalties)
            res = resources.tune_rr_penalty
            res.check_resource_existence()

            df = rr_training_ht_to_pandas(res.rr_training_ht.ht())
            tuning_df = tune_rr_penalty(df, args.rr_penalties)
            penalty = choose_rr_penalty(tuning_df)
            if not tuning_df.empty:
                logger.info(
                    "Held-out log loss per penalty, summed over contexts and"
                    " folds:\n%s",
                    tuning_df.groupby("penalty")["log_loss"].sum().to_string(),
                )
            logger.info("Chosen ridge penalty: %s", penalty)
            # The tuning Table is empty in a --test run, so give the schema
            # explicitly instead of inferring it from the DataFrame.
            ht = hl.Table.parallelize(
                tuning_df.to_dict("records"),
                hl.tstruct(
                    penalty=hl.tfloat64,
                    context=hl.tstr,
                    fold=hl.tstr,
                    n_test=hl.tint64,
                    log_loss=hl.tfloat64,
                ),
                key=["penalty", "context", "fold"],
            )
            ht = ht.annotate_globals(ridge_penalty=penalty)
            ht.write(res.rr_penalty_tuning_ht.path, overwrite=overwrite)
            logger.info("Done tuning the rr ridge penalty.")

        if args.train_rr_models:
            res = resources.train_rr_models
            res.check_resource_existence()
            penalty = hl.eval(res.rr_penalty_tuning_ht.ht().ridge_penalty)
            logger.info("Training the rr models (ridge penalty %s)...", penalty)

            df = rr_training_ht_to_pandas(res.rr_training_ht.ht())
            models = fit_rr_models(df, penalty)
            non_converged = [m["context"] for m in models if not m["converged"]]
            rr_models_to_ht(models).annotate_globals(
                ridge_penalty=penalty,
                non_converged_contexts=hl.literal(non_converged, hl.tarray(hl.tstr)),
            ).write(res.rr_models_ht.path, overwrite=overwrite)
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
            # Coalesce after the checkpoint so the filter of the per-variant Table
            # runs over all of its partitions.
            syn_ht = syn_ht.naive_coalesce(200)
            factors = compute_rr_normalization_factors(syn_ht, rr_ht)
            logger.info("rr normalization factors: %s", factors)

            rr_ht = normalize_rr_ht(rr_ht, factors)
            rr_ht = rr_ht.annotate_globals(
                rr_globals=hl.struct(
                    ridge_penalty=models_ht.index_globals().ridge_penalty,
                    non_converged_contexts=(
                        models_ht.index_globals().non_converged_contexts
                    ),
                    normalization_factors=hl.literal(
                        factors, hl.tdict(hl.tstr, hl.tfloat64)
                    ),
                    models_ht_path=res.rr_models_ht.path,
                    features_ht_path=res.genomic_features_ht.path,
                    per_variant_expected_ht_path=res.per_variant_expected_ht.path,
                    # The fitted models of the per-variant Table rr is normalized
                    # on; --apply-regional-rr checks they still match.
                    per_variant_apply_models_globals=(
                        per_variant_ht.index_globals().apply_models_globals
                    ),
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
        help=(
            "Post-fix to append to the output directory path. Also applies to the"
            " per-variant expected input, as in constraint_pipeline.py."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--path-post-fix",
        help=(
            "Post-fix to append to output file names. Also applies to the"
            " per-variant expected input, as in constraint_pipeline.py."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to filter the context, exomes sites, DNM, and per-variant"
            " expected Tables to only chromosome 20, chromosome X, and chromosome Y"
            " for testing. Reads the import_data.py --test outputs and the"
            " constraint_pipeline.py --test per-variant expected Table, and writes the"
            " outputs to the testing bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-gene-list",
        help=(
            "Whether to filter the context, exomes sites, DNM, and per-variant"
            " expected Tables to only six genes for testing (see --test-gene-list in"
            " constraint_pipeline.py). Implies --test."
        ),
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
            "Fit the per-context ridge-penalized rr models with the penalty chosen by"
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

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
