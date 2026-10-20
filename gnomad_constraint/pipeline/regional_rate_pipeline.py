"""Script to compute the regional mutation rate (rr) adjustment of expected counts.

Port of the Gnocchi regional adjustment (https://github.com/atgu/gnomad_nc_constraint)
retrained on synonymous de novo mutations (DNMs). For each strand-collapsed
trinucleotide context, a logistic model separates synonymous DNM sites from
non-mutated synonymous sites using regional genomic features, and the rr of a 1kb
window is the model's predicted probability for the window relative to the average
window. The constraint pipeline multiplies per-variant expected counts by rr after
model calibration (constraint_pipeline.py --aggregate-per-variant-expected).

The constraint pipeline's annotated context Table and per-variant expected Table are
read-only inputs; the constraint models are not re-fit or re-applied.

The rr pipeline consists of the following steps:
    - --import-genomic-features: import the 1kb regional genomic feature matrices.
    - --create-rr-training-set: create the DNM / non-mutated background training set.
    - --train-rr-models: run the univariate feature screen and fit the per-context
      models for each feature selection cutoff (optionally chromosome-split).
    - --compute-rr: compute rr per window and context (optionally normalized, with
      chrX/chrY feature masking, and out-of-fold rr).
    - --evaluate-rr-cutoffs: compare synonymous calibration across feature selection
      cutoffs using out-of-fold rr.
    - --compare-rr-variants: compare synonymous calibration of rr variants (e.g. the
      legacy syn_adj_r and the new rr with or without normalization, masking, or
      out-of-fold models) to isolate each change from the legacy rr.
"""

import argparse
import logging
from typing import Any, Dict, List, Optional

import hail as hl
import pandas as pd
import statsmodels
from gnomad.resources.resource_utils import TableResource
from hail.utils.misc import new_temp_file

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    RR_BACKGROUND_FRACTION,
    RR_BACKGROUND_SEED,
    RR_BONFERRONI_CUTOFF,
    RR_COMPARISON_VARIANTS,
    RR_STATSMODELS_VERSION,
    RR_TRAINING_N_PARTITIONS,
    RR_XY_MASKED_FEATURES,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import (
    filter_for_test,
    get_adj_r_ht,
    get_syn_adj_r_ht,
)
from gnomad_constraint.utils.regional_rate import (
    MODEL_FOLD_ALL,
    build_rr_variant_ht,
    compute_rr_ht,
    compute_rr_normalization_factors,
    create_rr_training_ht,
    evaluate_rr_cutoffs,
    evaluate_rr_variants,
    fit_rr_models,
    get_syn_mane_select_ht,
    get_syn_sites_ht,
    import_genomic_features,
    normalize_rr_ht,
    rr_models_to_ht,
    rr_training_ht_to_pandas,
    run_mask_ablation,
    run_rr_feature_screen,
    summarize_genomic_features,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_rate_pipeline")
logger.setLevel(logging.INFO)


def _write_tsv(df: pd.DataFrame, path: str) -> None:
    """
    Write a pandas DataFrame as a TSV to a local or GCS path.

    :param df: DataFrame to write.
    :param path: Output path.
    :return: None.
    """
    with hl.hadoop_open(path, "w") as f:
        df.to_csv(f, sep="\t", index=False)


def _read_tsv(path: str) -> pd.DataFrame:
    """
    Read a TSV from a local or GCS path into a pandas DataFrame.

    :param path: Input path.
    :return: DataFrame.
    """
    with hl.hadoop_open(path, "r") as f:
        return pd.read_csv(f, sep="\t")


def _get_rr_globals(
    cutoff: str,
    models: List[hl.Struct],
    args: argparse.Namespace,
    masked_features: Optional[List[str]],
    out_of_fold: bool,
    factors: Optional[Dict[str, float]],
    models_ht_path: str,
    features_ht_path: str,
) -> hl.expr.StructExpression:
    """
    Get the ``rr_globals`` global of an rr Table.

    :param cutoff: Feature selection cutoff token.
    :param models: rr model rows used for the cutoff.
    :param args: Parsed arguments.
    :param masked_features: Masked feature names, or None if not masked.
    :param out_of_fold: Whether rr was computed from the chromosome-split models.
    :param factors: Normalization factors, or None if rr was not normalized.
    :param models_ht_path: Path of the rr model Table.
    :param features_ht_path: Path of the genomic feature Table.
    :return: Struct expression of rr parameters.
    """
    calibrated = factors is not None
    return hl.struct(
        cutoff=cutoff,
        n_contexts_modeled=len({m.context for m in models if m.fold == MODEL_FOLD_ALL}),
        corrected_window_coords=args.corrected_window_coords,
        xy_masked=masked_features is not None,
        masked_features=hl.literal(masked_features or [], hl.tarray(hl.tstr)),
        out_of_fold=out_of_fold,
        calibrated=calibrated,
        normalization_factors=(
            hl.literal(factors, hl.tdict(hl.tstr, hl.tfloat64))
            if calibrated
            else hl.missing(hl.tdict(hl.tstr, hl.tfloat64))
        ),
        models_ht_path=models_ht_path,
        features_ht_path=features_ht_path,
    )


def _compute_and_write_rr(
    features_ht: hl.Table,
    models: List[hl.Struct],
    cutoff: str,
    out_path: str,
    args: argparse.Namespace,
    syn_ht: Optional[hl.Table],
    masked_features: Optional[List[str]],
    out_of_fold: bool,
    models_ht_path: str,
    features_ht_path: str,
) -> None:
    """
    Compute the rr Table of one cutoff, normalize it if requested, and write it.

    :param features_ht: Genomic feature Table.
    :param models: rr model rows of the cutoff.
    :param cutoff: Feature selection cutoff token.
    :param out_path: Output path of the rr Table.
    :param args: Parsed arguments.
    :param syn_ht: Synonymous MANE Select per-variant rows used for normalization,
        or None if safeguards are not applied.
    :param masked_features: Masked feature names, or None if not masked.
    :param out_of_fold: Whether to use the chromosome-split models.
    :param models_ht_path: Path of the rr model Table.
    :param features_ht_path: Path of the genomic feature Table.
    :return: None.
    """
    rr_ht = compute_rr_ht(features_ht, models, masked_features, out_of_fold)
    factors = None
    if args.apply_rr_safeguards:
        rr_ht = rr_ht.checkpoint(new_temp_file("rr_unnormalized", "ht"))
        factors = compute_rr_normalization_factors(syn_ht, rr_ht)
        logger.info(
            "rr normalization factors for cutoff %s%s: %s",
            cutoff,
            " (out-of-fold)" if out_of_fold else "",
            factors,
        )
        rr_ht = normalize_rr_ht(rr_ht, factors)

    rr_ht = rr_ht.annotate_globals(
        rr_globals=_get_rr_globals(
            cutoff,
            models,
            args,
            masked_features,
            out_of_fold,
            factors,
            models_ht_path,
            features_ht_path,
        )
    )
    rr_ht = rr_ht.checkpoint(out_path, overwrite=args.overwrite)

    # Log rr quantiles per context from the written Table.
    stats = rr_ht.aggregate(
        hl.agg.explode(
            lambda kv: hl.agg.group_by(
                kv[0], hl.agg.approx_quantiles(kv[1], [0.0, 0.05, 0.5, 0.95, 1.0])
            ),
            rr_ht.adj_r.items(),
        )
    )
    logger.info(
        "rr quantiles (min, 5%%, median, 95%%, max) per context for cutoff %s:\n%s",
        cutoff,
        "\n".join(f"    {c}: {q}" for c, q in sorted(stats.items())),
    )


def main(args):
    """Compute the regional mutation rate (rr) adjustment."""
    hl.init(
        log="/regional_rate_pipeline.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test_gene_list = args.test_gene_list
    test = args.test or test_gene_list
    overwrite = args.overwrite
    cutoffs = args.feature_pval_cutoffs
    masked_features = args.masked_features if args.mask_features else None

    if version not in VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    resources = constraint_res.get_regional_rate_resources(
        version,
        overwrite,
        test,
        directory_post_fix=args.directory_post_fix,
        path_post_fix=args.path_post_fix,
        output_root=args.output_root,
        corrected_window_coords=args.corrected_window_coords,
        xy_masked=args.mask_features,
        calibrated=args.apply_rr_safeguards,
        out_of_fold=args.out_of_fold,
        cutoffs=cutoffs,
        include_xy_features=args.include_xy_features,
        extra_features_paths=args.extra_features_paths,
        loeuf_ht=TableResource(args.loeuf_ht_path) if args.loeuf_ht_path else None,
    )

    syn_hts: Dict[str, hl.Table] = {}

    def _get_syn_ht(ht: hl.Table) -> hl.Table:
        """
        Get the checkpointed synonymous MANE Select rows of the per-variant Table.

        The per-variant expected Table is filtered (and filtered for testing if
        requested) and checkpointed once, then reused by --compute-rr and
        --evaluate-rr-cutoffs.

        :param ht: Per-variant expected Table.
        :return: Synonymous MANE Select per-variant rows.
        """
        if "syn" not in syn_hts:
            if test:
                ht = filter_for_test(ht, use_gene_list=test_gene_list)
            syn_hts["syn"] = get_syn_mane_select_ht(ht).checkpoint(
                new_temp_file("rr_syn_mane_select", "ht")
            )

        return syn_hts["syn"]

    try:
        if args.import_genomic_features:
            logger.info("Importing regional genomic feature matrices...")
            res = resources.import_genomic_features
            res.check_resource_existence()

            paths = res.input_resources["genomic feature matrices"]
            ht = import_genomic_features(paths, args.corrected_window_coords)
            ht = ht.checkpoint(res.genomic_features_ht.path, overwrite=overwrite)

            summary = summarize_genomic_features(ht)
            if summary["n_duplicated_windows"]:
                raise ValueError(
                    f"{summary['n_duplicated_windows']} windows are duplicated across"
                    " the imported feature matrices!"
                )
            for contig, s in sorted(summary["contigs"].items()):
                logger.info(
                    "    %s: %d windows; features missing in all windows: %s",
                    contig,
                    s["n_windows"],
                    s["missing_features"],
                )
            logger.info("Done importing regional genomic feature matrices.")

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
                args.background_fraction,
                args.background_seed,
            ).checkpoint(new_temp_file("rr_training", "ht"))
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
            logger.info("Training rr models for cutoffs: %s...", cutoffs)
            res = resources.train_rr_models
            res.check_resource_existence()

            # The legacy syn_adj_r models were fit with statsmodels 0.14.4 (which
            # requires scipy < 1.16). fit_regularized results differ slightly
            # (~1e-5 relative rr) with statsmodels >= 0.14.5.
            if statsmodels.__version__ != RR_STATSMODELS_VERSION:
                logger.warning(
                    "statsmodels %s is installed. Use statsmodels==%s (with scipy <"
                    " 1.16) to exactly reproduce the legacy syn_adj_r rr models.",
                    statsmodels.__version__,
                    RR_STATSMODELS_VERSION,
                )
            df = rr_training_ht_to_pandas(res.rr_training_ht.ht())
            screen_df = run_rr_feature_screen(df)
            _write_tsv(screen_df, res.rr_feature_screen_tsv)

            models, metrics_df = fit_rr_models(
                df, screen_df, cutoffs, cv_chrom_split=args.cv_chrom_split
            )
            rr_models_to_ht(models).write(res.rr_models_ht.path, overwrite=overwrite)
            _write_tsv(metrics_df, res.rr_cv_metrics_tsv)
            if not metrics_df.empty:
                logger.info(
                    "Held-out metrics per cutoff (mean over contexts):\n%s",
                    metrics_df.assign(
                        ll_gain=(metrics_df["null_log_loss"] - metrics_df["log_loss"])
                        * metrics_df["n_test"]
                    )
                    .groupby("cutoff")
                    .agg(
                        n_contexts=("context", "nunique"),
                        heldout_ll_gain=("ll_gain", "sum"),
                        mean_auc=("auc", "mean"),
                    )
                    .to_string(),
                )
            logger.info("Done training rr models.")

        if args.compute_rr:
            logger.info("Computing rr per window and context...")
            res = resources.compute_rr
            res.check_resource_existence()

            models_ht = res.rr_models_ht.ht()
            all_models = models_ht.collect()
            features_ht = res.genomic_features_ht.ht()
            syn_ht = (
                _get_syn_ht(res.per_variant_expected_ht.ht())
                if args.apply_rr_safeguards
                else None
            )

            for cutoff in cutoffs:
                models = [m for m in all_models if m.cutoff == cutoff]
                if not models:
                    raise ValueError(
                        f"No rr models were trained for cutoff {cutoff}! Run"
                        " --train-rr-models with this cutoff first."
                    )
                kwargs: Dict[str, Any] = {
                    "features_ht": features_ht,
                    "models": models,
                    "cutoff": cutoff,
                    "args": args,
                    "syn_ht": syn_ht,
                    "masked_features": masked_features,
                    "models_ht_path": res.rr_models_ht.path,
                    "features_ht_path": res.genomic_features_ht.path,
                }
                _compute_and_write_rr(
                    out_path=getattr(res, f"rr_ht_{cutoff}").path,
                    out_of_fold=False,
                    **kwargs,
                )
                if args.out_of_fold:
                    _compute_and_write_rr(
                        out_path=getattr(res, f"rr_oof_ht_{cutoff}").path,
                        out_of_fold=True,
                        **kwargs,
                    )
                if args.mask_ablation:
                    ablation_df = run_mask_ablation(
                        features_ht, models, args.masked_features
                    )
                    logger.info(
                        "Masking ablation on autosomal windows for cutoff %s:\n%s",
                        cutoff,
                        ablation_df.to_string(index=False),
                    )
            logger.info("Done computing rr.")

        if args.evaluate_rr_cutoffs:
            logger.info("Evaluating rr feature selection cutoffs: %s...", cutoffs)
            res = resources.evaluate_rr_cutoffs
            res.check_resource_existence()

            rr_hts = {c: getattr(res, f"rr_oof_ht_{c}").ht() for c in cutoffs}
            n_contexts_modeled = {
                c: hl.eval(ht.rr_globals.n_contexts_modeled) for c, ht in rr_hts.items()
            }
            df = evaluate_rr_cutoffs(
                _get_syn_ht(res.per_variant_expected_ht.ht()),
                rr_hts,
                cv_metrics_df=_read_tsv(res.rr_cv_metrics_tsv),
                n_contexts_modeled=n_contexts_modeled,
                loeuf_ht=res.loeuf_ht.ht() if args.loeuf_ht_path else None,
                syn_z_outlier_threshold=args.syn_z_outlier_threshold,
            )
            _write_tsv(df, res.rr_cutoff_evaluation_tsv)
            logger.info("rr cutoff evaluation:\n%s", df.to_string(index=False))
            logger.info("Done evaluating rr feature selection cutoffs.")

        if args.compare_rr_variants:
            logger.info("Comparing rr variants: %s...", args.rr_variants)
            res = resources.compare_rr_variants
            res.check_resource_existence()

            syn_ht = _get_syn_ht(res.per_variant_expected_ht.ht())
            models = res.rr_models_ht.ht().collect()
            features_ht = res.genomic_features_ht.ht()
            legacy_hts = {
                "legacy": get_syn_adj_r_ht(),
                "legacy_noncoding": get_adj_r_ht(),
            }
            rr_hts = {}
            for spec in args.rr_variants:
                rr_hts[spec], factors = build_rr_variant_ht(
                    spec,
                    syn_ht,
                    features_ht=features_ht,
                    models=models,
                    legacy_hts=legacy_hts,
                    masked_features=args.masked_features,
                )
                if factors is not None:
                    logger.info("rr normalization factors for %s: %s", spec, factors)

            df, context_df = evaluate_rr_variants(
                syn_ht,
                rr_hts,
                loeuf_ht=res.loeuf_ht.ht() if args.loeuf_ht_path else None,
                syn_z_outlier_threshold=args.syn_z_outlier_threshold,
            )
            _write_tsv(df, res.rr_variant_comparison_tsv)
            _write_tsv(context_df, res.rr_variant_context_rr_tsv)
            logger.info("rr variant comparison:\n%s", df.to_string(index=False))
            logger.info(
                "Per-context expected-weighted mean rr:\n%s",
                context_df.to_string(index=False),
            )
            logger.info("Done comparing rr variants.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            constraint_res.get_logging_path("regional_rate_pipeline", version=version)
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
        "--output-root",
        help=(
            "Optional bucket path (e.g. gs://gnomad-tmp/rp) to use as the root of the"
            " outputs instead of the production root. The constraint pipeline inputs"
            " (annotated context and per-variant expected Tables) are always read"
            " from their production paths."
        ),
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
            " Tables to only a list of genes on chromosome 20, chromosome X, and"
            " chromosome Y for testing."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--corrected-window-coords",
        help=(
            "Use the corrected 1kb window coordinates ([start + 1, end] for BED"
            " windows) instead of the legacy coordinates ([start, end), 1bp left of"
            " the BED window). The legacy coordinates are required to reproduce the"
            " existing syn_adj_r Tables. Applies to all steps."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--feature-pval-cutoffs",
        help=(
            "Feature selection cutoff tokens: 'bonferroni' (legacy rule) and/or raw"
            " p-value cutoffs (e.g. 1e-3 0.01). Used by --train-rr-models,"
            f" --compute-rr, and --evaluate-rr-cutoffs. Default is"
            f" {RR_BONFERRONI_CUTOFF}."
        ),
        nargs="+",
        default=[RR_BONFERRONI_CUTOFF],
    )

    import_args = parser.add_argument_group(
        "Import genomic features args",
        "Arguments used for importing the regional genomic feature matrices.",
    )
    import_args.add_argument(
        "--import-genomic-features",
        help="Import the 1kb regional genomic feature matrices.",
        action="store_true",
    )
    import_args.add_argument(
        "--include-xy-features",
        help=(
            "Import the chrX/chrY feature matrix built by build_xy_features.py with"
            " the published (chr1-22) feature matrix."
        ),
        action="store_true",
    )
    import_args.add_argument(
        "--extra-features-paths",
        help="Additional feature matrix paths to import.",
        nargs="+",
        default=None,
    )

    training_set_args = parser.add_argument_group(
        "rr training set args", "Arguments used for creating the rr training set."
    )
    training_set_args.add_argument(
        "--create-rr-training-set",
        help=(
            "Create the rr training set of synonymous DNM sites and a sampled"
            " background of non-mutated synonymous sites."
        ),
        action="store_true",
    )
    training_set_args.add_argument(
        "--background-fraction",
        help=(
            "Fraction of non-mutated synonymous sites sampled as background. Default"
            f" is {RR_BACKGROUND_FRACTION}."
        ),
        type=float,
        default=RR_BACKGROUND_FRACTION,
    )
    training_set_args.add_argument(
        "--background-seed",
        help=f"Seed used to sample the background. Default is {RR_BACKGROUND_SEED}.",
        type=int,
        default=RR_BACKGROUND_SEED,
    )

    train_args = parser.add_argument_group(
        "Train rr models args", "Arguments used for training the rr models."
    )
    train_args.add_argument(
        "--train-rr-models",
        help=(
            "Run the univariate feature screen and fit the per-context rr models for"
            " each feature selection cutoff."
        ),
        action="store_true",
    )
    train_args.add_argument(
        "--cv-chrom-split",
        help=(
            "Also fit models with odd or even chromosomes held out and record the"
            " held-out log loss and AUC. Required for --compute-rr --out-of-fold."
        ),
        action="store_true",
    )

    compute_args = parser.add_argument_group(
        "Compute rr args", "Arguments used for computing rr per window and context."
    )
    compute_args.add_argument(
        "--compute-rr",
        help="Compute rr per 1kb window and context for each cutoff.",
        action="store_true",
    )
    compute_args.add_argument(
        "--apply-rr-safeguards",
        help=(
            "Divide rr by its per-context expected-weighted mean over exome"
            " synonymous MANE Select sites. Without this flag rr is written as"
            " trained (reproducing the legacy syn_adj_r)."
        ),
        action="store_true",
    )
    compute_args.add_argument(
        "--mask-features",
        help=(
            "On chrX/chrY windows, pin --masked-features to the genome-average"
            " window (z = 0) instead of requiring them. Autosomes are not masked."
        ),
        action="store_true",
    )
    compute_args.add_argument(
        "--masked-features",
        help=f"Features masked with --mask-features. Default is {RR_XY_MASKED_FEATURES}.",
        nargs="+",
        default=RR_XY_MASKED_FEATURES,
    )
    compute_args.add_argument(
        "--mask-ablation",
        help=(
            "Log, for autosomal windows, how rr changes when --masked-features are"
            " masked."
        ),
        action="store_true",
    )
    compute_args.add_argument(
        "--out-of-fold",
        help=(
            "Also compute out-of-fold rr from the chromosome-split models (requires"
            " --train-rr-models --cv-chrom-split). Used by --evaluate-rr-cutoffs."
        ),
        action="store_true",
    )

    evaluate_args = parser.add_argument_group(
        "Evaluate rr cutoffs args",
        "Arguments used for evaluating the feature selection cutoffs.",
    )
    evaluate_args.add_argument(
        "--evaluate-rr-cutoffs",
        help=(
            "Evaluate synonymous calibration (dispersion, outliers, o/e by LOEUF"
            " decile) of the out-of-fold rr of each cutoff."
        ),
        action="store_true",
    )
    evaluate_args.add_argument(
        "--loeuf-ht-path",
        help=(
            "Constraint Table (metrics or release format) used to bin transcripts by"
            " LOEUF decile. Default is the v4.1.2 constraint metrics Table"
            f" ({constraint_res.get_loeuf_decile_ht().path}), which fixes the v4.1.1"
            " LOEUF bins."
        ),
        type=str,
        default=constraint_res.get_loeuf_decile_ht().path,
    )
    evaluate_args.add_argument(
        "--syn-z-outlier-threshold",
        help="Absolute raw syn z-score defining an outlier. Default is 8.0.",
        type=float,
        default=8.0,
    )

    compare_args = parser.add_argument_group(
        "Compare rr variants args",
        "Arguments used for comparing rr variants. Uses --loeuf-ht-path,"
        " --syn-z-outlier-threshold, and --masked-features.",
    )
    compare_args.add_argument(
        "--compare-rr-variants",
        help=(
            "Compare synonymous calibration (dispersion, outliers, o/e by LOEUF"
            " decile, per-context mean rr) of the rr variants in --rr-variants and"
            " the unadjusted expected."
        ),
        action="store_true",
    )
    compare_args.add_argument(
        "--rr-variants",
        help=(
            "rr variant specs to compare: a source ('legacy', 'legacy_noncoding',"
            " 'path=<rr Table>', or 'cutoff=<token>') followed by optional"
            " comma-separated modifiers ('norm', and, for 'cutoff=', 'mask' and"
            " 'oof'), e.g. 'cutoff=bonferroni,norm,mask'. Default is"
            f" {RR_COMPARISON_VARIANTS}."
        ),
        nargs="+",
        default=RR_COMPARISON_VARIANTS,
    )

    args = parser.parse_args()
    main(args)
