"""Export the per-context rr logistic fits of the background runs, for the report figure.

The pipeline version of the local ``rr_logistic_fits.py`` background tables, from the
runs of ``run_backgrounds.sh``. Every site is scored in Hail with the chromosome-split
model that didn't see it (ridge, default penalty):

- bins are deciles of the possible synonymous sites (the all-sites run's background, one
  row per locus) by their score under the current background's model (absent from
  gnomAD, one row per locus);
- per context, background and bin: DNMs, possible sites, their mean score under the
  current background's model, the background's training rows, the observed and fitted
  share of DNMs among DNMs plus that background, the background's mean rr over the
  possible sites (normalized to 1 per context), and the DNM rate per possible site
  relative to the context's.

Writes ``<prefix>.bins.tsv`` and ``<prefix>.contexts.tsv`` in the local version's format.
"""

import argparse
from typing import Dict, List, Tuple

import hail as hl
import numpy as np
import pandas as pd

from gnomad_constraint.resources.constants import CPG_CONTEXTS, GENOMIC_FEATURE_COLUMNS
from gnomad_constraint.utils.regional_rate import RR_FOLDS, _chrom_fold_expr

ROOT = "gs://gnomad-tmp/jg/rr_bg"
RD = "gnomad_v4.1.1/constraint/regional_rate"
BACKGROUNDS = {
    "absent": "absent_locus",
    "absent_context_matched": "ctx_proportional",
    "absent_equal_per_context": "ctx_even",
    "all_sites": "all_locus",
}
"""Background name in the figure -> run in ``ROOT``."""
REFERENCE = "absent"
UNIVERSE = "all_sites"
N_BINS = 10


def class_of(context: str) -> str:
    """Get the context class of a strand-collapsed trinucleotide context."""
    if context[1] == "A":
        return "A/T"
    return "CpG" if context in CPG_CONTEXTS else "C non-CpG"


def read_models(run: str, cutoff: str) -> Dict[Tuple[str, str], hl.Struct]:
    """Get a run's chromosome-split models, keyed by (context, fold)."""
    ht = hl.read_table(f"{ROOT}/{run}/{RD}/gnomad.v4.1.1.rr_models.ht")
    return {
        (m.context, m.fold): m
        for m in ht.collect()
        if m.cutoff == cutoff and m.fold in RR_FOLDS
    }


def oof_exprs(
    x: hl.expr.ArrayExpression,
    key: hl.expr.StringExpression,
    models: Dict[Tuple[str, str], hl.Struct],
) -> Tuple[hl.expr.Float64Expression, hl.expr.Float64Expression]:
    """
    Get each row's out-of-fold score (β·z) and intercept.

    The models are looked up from a literal keyed by "context:fold", with each
    model's features as indices into ``x``, so the generated code stays small. The
    score is missing where the row's context has no model or a model feature is
    missing or NaN.

    :param x: Each row's ``GENOMIC_FEATURE_COLUMNS`` values.
    :param key: Each row's "context:fold".
    :param models: Output of :func:`read_models`.
    :return: (score, intercept) expressions.
    """
    lit = hl.literal(
        {
            f"{c}:{f}": hl.Struct(
                idx=[GENOMIC_FEATURE_COLUMNS.index(feat) for feat in m.features],
                coef=list(m.coefs),
                mean=list(m.mean),
                sd=list(m.sd),
                intercept=m.intercept,
            )
            for (c, f), m in models.items()
        },
        dtype=hl.tdict(
            hl.tstr,
            hl.tstruct(
                idx=hl.tarray(hl.tint32),
                coef=hl.tarray(hl.tfloat64),
                mean=hl.tarray(hl.tfloat64),
                sd=hl.tarray(hl.tfloat64),
                intercept=hl.tfloat64,
            ),
        ),
    )
    m = lit.get(key)
    s = hl.sum(
        hl.range(hl.len(m.idx)).map(
            lambda j: m.coef[j] * (x[m.idx[j]] - m.mean[j]) / m.sd[j]
        ),
        filter_missing=False,
    )
    return hl.or_missing(hl.is_finite(s), s), m.intercept


def scored(run: str, models: Dict[str, Dict], names: List[str]) -> hl.Table:
    """
    Read a run's training Table and score it with the given models.

    :param run: Run in ``ROOT``.
    :param models: Background name -> output of :func:`read_models`.
    :param names: Backgrounds whose models to score with.
    :return: Table with ``context``, ``y`` and ``s_<name>`` / ``b0_<name>`` per name.
    """
    ht = hl.read_table(f"{ROOT}/{run}/{RD}/gnomad.v4.1.1.rr_training.ht")
    ht = ht.annotate(
        feat_values=hl.array([hl.float64(ht[f]) for f in GENOMIC_FEATURE_COLUMNS]),
        context_fold=ht.context + ":" + _chrom_fold_expr(ht.locus.contig),
    )
    fields = {}
    for name in names:
        s, b0 = oof_exprs(ht.feat_values, ht.context_fold, models[name])
        fields[f"s_{name}"] = s
        fields[f"b0_{name}"] = b0
    ht = ht.select("context", "y", **fields).key_by()
    return ht.checkpoint(hl.utils.new_temp_file(f"bg_fits.{run}", "ht"))


def sigmoid(x: hl.expr.Float64Expression) -> hl.expr.Float64Expression:
    """Get the logistic function of x."""
    return 1.0 / (1.0 + hl.exp(-x))


def main(args: argparse.Namespace) -> None:
    """Score the background runs and write the tables."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    models = {b: read_models(run, args.cutoff) for b, run in BACKGROUNDS.items()}

    # Possible sites and DNMs, scored by every background's model.
    univ = scored(BACKGROUNDS[UNIVERSE], models, list(BACKGROUNDS))
    univ = univ.filter(hl.is_defined(univ[f"s_{REFERENCE}"]))
    edges = univ.aggregate(
        hl.agg.filter(
            univ.y == 0,
            hl.agg.group_by(
                univ.context,
                hl.agg.approx_quantiles(
                    univ[f"s_{REFERENCE}"],
                    [float(q) for q in np.linspace(0, 1, N_BINS + 1)[1:-1]],
                ),
            ),
        )
    )
    edges_lit = hl.literal(edges)

    def bin_of(ht: hl.Table) -> hl.expr.Int32Expression:
        s = ht[f"s_{REFERENCE}"]
        return hl.len(edges_lit[ht.context].filter(lambda e: s >= e))

    univ = univ.annotate(bin=bin_of(univ))
    aggs = {
        "dnms": hl.agg.count_where(univ.y == 1),
        "possible": hl.agg.count_where(univ.y == 0),
        "score_sum": hl.agg.filter(univ.y == 0, hl.agg.sum(univ[f"s_{REFERENCE}"])),
    }
    for b in BACKGROUNDS:
        eta = univ[f"b0_{b}"] + univ[f"s_{b}"]
        aggs[f"rr_sum_{b}"] = hl.agg.filter(
            univ.y == 0, hl.agg.sum(sigmoid(eta) / sigmoid(univ[f"b0_{b}"]))
        )
        aggs[f"rr_n_{b}"] = hl.agg.filter(
            (univ.y == 0) & hl.is_defined(eta), hl.agg.count()
        )
        aggs[f"fit_dnm_{b}"] = hl.agg.filter(univ.y == 1, hl.agg.sum(sigmoid(eta)))
    base = univ.group_by("context", "bin").aggregate(**aggs).to_pandas()

    # Each background's own training rows, binned by the current background's model.
    rows = []
    for b, run in BACKGROUNDS.items():
        ht = univ if b == UNIVERSE else scored(run, models, [REFERENCE, b])
        ht = ht.filter((ht.y == 0) & hl.is_defined(ht[f"s_{REFERENCE}"]))
        ht = ht.annotate(bin=bin_of(ht))
        bg = (
            ht.group_by("context", "bin")
            .aggregate(
                background_weight=hl.agg.count(),
                fit_bg=hl.agg.sum(sigmoid(ht[f"b0_{b}"] + ht[f"s_{b}"])),
            )
            .to_pandas()
        )
        d = base.merge(bg, on=["context", "bin"], how="left").fillna(
            {"background_weight": 0, "fit_bg": 0}
        )
        tot = d.groupby("context")[
            [f"rr_sum_{b}", f"rr_n_{b}", "dnms", "possible"]
        ].sum()
        rr_norm = d.context.map(tot[f"rr_sum_{b}"] / tot[f"rr_n_{b}"])
        rows.append(
            pd.DataFrame(
                {
                    "context": d.context,
                    "background": b,
                    "bin": d.bin,
                    "score": d.score_sum / d.possible,
                    "dnms": d.dnms,
                    "possible": d.possible,
                    "background_weight": d.background_weight,
                    "observed_p": d.dnms / (d.dnms + d.background_weight),
                    "fitted_p": (d[f"fit_dnm_{b}"] + d.fit_bg)
                    / (d.dnms + d.background_weight),
                    "rr": d[f"rr_sum_{b}"] / d[f"rr_n_{b}"] / rr_norm,
                    "dnm_rate_relative": (d.dnms / d.possible)
                    / d.context.map(tot.dnms / tot.possible),
                }
            )
        )
    bins = pd.concat(rows).sort_values(["context", "background", "bin"])
    contexts = pd.DataFrame({"context": sorted(bins.context.unique())}).assign(
        **{"class": lambda c: c.context.map(class_of)}
    )
    for name, table in (("bins", bins), ("contexts", contexts)):
        with hl.hadoop_open(f"{args.out_prefix}.{name}.tsv", "w") as f:
            table.to_csv(f, sep="\t", index=False)
    print(
        bins[bins.context.isin(["CAG", "CCT", "CCG"])].round(4).to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff", default="ridge")
    parser.add_argument(
        "--out-prefix",
        default=f"{ROOT}/logistic/rr_logistic_fits.pipeline_backgrounds",
        help="GCS path prefix.",
    )
    main(parser.parse_args())
