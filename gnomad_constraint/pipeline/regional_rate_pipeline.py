"""
Script to compute the regional mutation rate (rr) of each 1kb window and context.

For each strand-collapsed trinucleotide context, a ridge-penalized logistic model
separates synonymous de novo mutation (DNM) sites from non-mutated synonymous sites
using regional genomic features. A window's rr is the model's prediction for the
window relative to the average window, normalized so that each context's
expected-weighted mean rr over exome synonymous sites is 1.
``constraint_pipeline.py --aggregate-per-variant-expected --apply-regional-rr``
applies it.

The DNM and genomic feature Tables are written by ``import_data.py``. The
per-variant expected Table is read with the same --version, --test,
--directory-post-fix, and --path-post-fix as the rr outputs, the values
``--apply-regional-rr`` must be run with, so rr is normalized on the Table it is
applied to.

The rr pipeline consists of the following steps:
    - --create-rr-training-set: create the DNM and non-mutated background training
      set, leaving out genes with a QC gene flag.
    - --tune-rr-penalty: choose the ridge penalty by the held-out log loss of
      chromosome-split models.
    - --train-rr-models: fit the per-context rr models.
    - --compute-rr: compute the normalized rr per window and context.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.regional_mutation_rate import (
    choose_regional_mutation_rate_penalty,
    compute_regional_mutation_rate_ht,
    fit_regional_mutation_rate_models,
    normalize_regional_mutation_rate_expr,
    regional_mutation_rate_models_to_ht,
    regional_mutation_rate_normalization_factors_agg_expr,
    tune_regional_mutation_rate_penalty,
)
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    CURRENT_VERSION,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_CONTEXT_FEATURES,
    RR_CONTEXTS,
    RR_RIDGE_PENALTIES,
    RR_RIDGE_PENALTY,
    RR_TRAINING_N_PARTITIONS,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import filter_for_test
from gnomad_constraint.utils.constraint import get_transcript_filter_expr
from gnomad_constraint.utils.regional_rate import (
    create_rr_training_ht,
    remove_qc_flagged_genes,
    rr_training_ht_to_pandas,
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
        if args.create_rr_training_set:
            logger.info("Creating the rr training set...")
            res = resources.create_rr_training_set
            res.check_resource_existence()

            context_ht = res.annotated_context_ht.ht()
            if test:
                context_ht = filter_for_test(context_ht)
            ht = create_rr_training_ht(
                context_ht,
                res.dnm_ht.ht(),
                res.genomic_features_ht.ht(),
                RR_BACKGROUND_FRACTION,
                RR_BACKGROUND_SEED,
            )
            ht = remove_qc_flagged_genes(ht, res.gene_flags_ht.ht())
            ht = ht.checkpoint(new_temp_file("rr_training", "ht"))
            ht = ht.naive_coalesce(RR_TRAINING_N_PARTITIONS)
            ht.write(res.rr_training_ht.path, overwrite=overwrite)
            logger.info("Done creating the rr training set.")

        if args.tune_rr_penalty:
            logger.info("Tuning the rr ridge penalty over %s...", args.rr_penalties)
            res = resources.tune_rr_penalty
            res.check_resource_existence()

            df = rr_training_ht_to_pandas(res.rr_training_ht.ht())
            tuning_df = tune_regional_mutation_rate_penalty(
                df, RR_CONTEXT_FEATURES, args.rr_penalties
            )
            penalty = choose_regional_mutation_rate_penalty(tuning_df, RR_RIDGE_PENALTY)
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
            models = fit_regional_mutation_rate_models(df, RR_CONTEXT_FEATURES, penalty)
            ht = regional_mutation_rate_models_to_ht(models)
            ht = ht.annotate_globals(
                ridge_penalty=penalty,
                non_converged_contexts=hl.literal(
                    [m["context"] for m in models if not m["converged"]],
                    hl.tarray(hl.tstr),
                ),
            )
            ht.write(res.rr_models_ht.path, overwrite=overwrite)
            logger.info("Done training the rr models.")

        if args.compute_rr:
            logger.info("Computing rr per window and context...")
            res = resources.compute_rr
            res.check_resource_existence()

            models_ht = res.rr_models_ht.ht()
            rr_ht = compute_regional_mutation_rate_ht(
                res.genomic_features_ht.ht(), models_ht.collect(), RR_CONTEXTS
            )
            rr_ht = rr_ht.checkpoint(new_temp_file("rr_unnormalized", "ht"))

            # Normalize rr on the synonymous Ensembl MANE Select sites, weighted by
            # their adj expected count. The RefSeq MANE Select transcripts are the
            # same transcripts, so keeping them would count every site twice.
            per_variant_ht = res.per_variant_expected_ht.ht()
            adj_idx = hl.eval(per_variant_ht.exomes_freq_meta).index(ADJ_FREQ_META)
            syn_ht = per_variant_ht.filter(
                (per_variant_ht.annotation == "synonymous_variant")
                & get_transcript_filter_expr(per_variant_ht, mane_select_only=True)
            )
            syn_ht = syn_ht.select("context", exp=syn_ht.expected_variants[adj_idx])
            syn_ht = syn_ht.checkpoint(new_temp_file("rr_syn_mane_select", "ht"))
            syn_ht = syn_ht.naive_coalesce(200)
            factors = syn_ht.aggregate(
                regional_mutation_rate_normalization_factors_agg_expr(
                    syn_ht.context,
                    rr_ht[syn_ht.locus].adj_r.get(syn_ht.context),
                    syn_ht.exp,
                )
            )
            logger.info("rr normalization factors: %s", factors)

            rr_ht = rr_ht.annotate(
                adj_r=normalize_regional_mutation_rate_expr(rr_ht.adj_r, factors)
            )
            rr_ht = rr_ht.annotate_globals(
                rr_globals=hl.struct(
                    ridge_penalty=models_ht.index_globals().ridge_penalty,
                    normalization_factors=hl.literal(
                        factors, hl.tdict(hl.tstr, hl.tfloat64)
                    ),
                    # --apply-regional-rr checks that the per-variant Table it
                    # aggregates still has these models.
                    per_variant_apply_models_globals=(
                        per_variant_ht.index_globals().apply_models_globals
                    ),
                )
            )
            rr_ht.write(res.regional_rr_ht.path, overwrite=overwrite)
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
            "Whether to filter the context Table to only chromosome 20, chromosome X,"
            " and chromosome Y for testing. Reads the import_data.py --test outputs"
            " and the constraint_pipeline.py --test per-variant expected Table, and"
            " writes the outputs to the testing bucket."
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
