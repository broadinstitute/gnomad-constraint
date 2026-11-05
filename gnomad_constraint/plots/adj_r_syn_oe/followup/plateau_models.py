"""Refit the v4.1.1 plateau models a few ways, to see where the context offsets come from.

The plateau models turn a mutation rate into the expected share of possible variants
seen, with one straight line for CpG transitions and one for every other mutation
type (per genomic region), fit by unweighted least squares over the mutation-rate
groups (context, ref, alt, methylation level) of the high-coverage synonymous
training sites. Inputs: local copies of the training Table
(``gs://gnomad/v4.1.1/constraint/training_data/gnomad.v4.1.1.constraint_training.ht``)
and the production plateau models (``.../models/gnomad.v4.1.1.plateau.he``).

On the autosomes and PARs, the versions compared are:

- ``stored``: the production models as stored;
- ``refit``: the same two lines refit from the training Table as the pipeline does
  (unweighted). The CpG line matches the stored one exactly; the other comes out
  slightly different (intercept 0.031 against 0.034), so the training Table has
  probably changed a little since the models were built;
- ``weighted``: the same two lines weighted by possible variants
  (``build_models(weighted=True)``, the pipeline's ``--use-weights``), which makes
  total observed equal total expected within each line;
- ``by_class``: one weighted line per context class (A/T, other C, CpG transitions;
  CpG transversions go with other C, as in the non-CpG line now);
- ``by_substitution``: one weighted line per substitution type (A>C, A>G, A>T, C>A,
  C>G, C>T outside CpG, and CpG transitions);
- ``per_group``: each mutation-rate group calibrated to its own observed share (no
  line; expected = observed on the training sites by construction).

Writes ``<prefix>.groups.tsv`` (each group's counts and expected count under each
version), ``<prefix>.fits.tsv`` (every line's intercept and slope) and
``<prefix>.contexts.tsv`` (each strand-collapsed context's obs/exp). With ``--cells``
and the other evaluation-export inputs, it also scales each exported cell's expected
count by its context's ratio of new to stored expected (training sites) and writes
the per-gene noise with no rr and with the ridge rr to ``<prefix>.noise.tsv``.
"""

import argparse
import logging

import hail as hl
import numpy as np
import pandas as pd
from full_run_scores import cell_rr, tx_expected
from rr_shrinkage_experiments import Evaluator

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

REGION = "autosome_or_par"
FITS = ("stored", "refit", "weighted", "by_class", "by_substitution", "per_group")
LINES = {
    "refit": ("cpg", False),
    "weighted": ("cpg", True),
    "by_class": ("cls", True),
    "by_substitution": ("sub", True),
}
"""Line versions: grouping column and whether the fit is weighted."""
COMP = str.maketrans("ACGT", "TGCA")


def export_groups(training_ht: str) -> pd.DataFrame:
    """Sum the high-coverage training sites per mutation-rate group."""
    hl.init(quiet=True)
    ht = hl.read_table(training_ht)
    ht = ht.filter(
        (ht.build_model.high_or_low_coverage == "high") & (ht.genomic_region == REGION)
    )
    ht = ht.group_by("context", "ref", "alt", "methylation_level", "cpg").aggregate(
        observed=hl.agg.sum(ht.observed_variants[0]),
        possible=hl.agg.sum(ht.possible_variants),
        mu=hl.agg.mean(ht.mu_snp),
    )
    return ht.to_pandas()


def annotate_groups(g: pd.DataFrame) -> pd.DataFrame:
    """Add the observed share, line groups and strand-collapsed context."""
    g["prop"] = g.observed / g.possible
    flip = ~g.ref.isin(["A", "C"])
    ref = g.ref.where(~flip, g.ref.str.translate(COMP))
    alt = g.alt.where(~flip, g.alt.str.translate(COMP))
    g["ctx"] = g.context.where(~flip, g.context.str.translate(COMP).str[::-1])
    g["cls"] = np.where(g.cpg, "CpG transition", np.where(ref == "A", "A/T", "other C"))
    g["sub"] = np.where(g.cpg, "CpG transition", ref + ">" + alt)
    return g


def fit_line(x: np.ndarray, y: np.ndarray, w: np.ndarray = None) -> tuple:
    """Fit y = intercept + slope * x by (weighted) least squares."""
    w = np.ones_like(x) if w is None else w
    X = np.column_stack([np.ones_like(x), x])
    beta = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (y * w))
    return float(beta[0]), float(beta[1])


def fit_versions(g: pd.DataFrame, stored: dict) -> pd.DataFrame:
    """Add each version's expected count to the groups and return the lines."""
    fits = []
    g["exp_stored"] = np.nan
    for k, line in stored.items():
        if k.genomic_region == REGION:
            b0, b1 = line[0]
            label = "CpG transition" if k.cpg else "everything else"
            fits.append({"fit": "stored", "line": label, "intercept": b0, "slope": b1})
            m = g.cpg == k.cpg
            g.loc[m, "exp_stored"] = g.possible[m] * (b0 + b1 * g.mu[m])
    for fit, (by, weighted) in LINES.items():
        g[f"exp_{fit}"] = np.nan
        for key, d in g.groupby(by):
            b0, b1 = fit_line(
                d.mu.to_numpy(),
                d.prop.to_numpy(),
                d.possible.to_numpy().astype(float) if weighted else None,
            )
            if by == "cpg":
                key = "CpG transition" if key else "everything else"
            fits.append({"fit": fit, "line": key, "intercept": b0, "slope": b1})
            g.loc[d.index, f"exp_{fit}"] = d.possible * (b0 + b1 * d.mu)
    g["exp_per_group"] = g.observed.astype(float)
    return pd.DataFrame(fits)


def downstream_noise(args: argparse.Namespace, ctx: pd.DataFrame) -> pd.DataFrame:
    """Score per-gene noise on the exported cells under each plateau version."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    exp_cell = ev.exp_cell.copy()
    outcome_free = ev.no_qc_flag & (ev.exp0 >= 20)
    rows = []
    for fit in FITS:
        factor = (ctx[f"exp_{fit}"] / ctx.exp_stored).reindex(ev.cells.context)
        ev.exp_cell = exp_cell * factor.fillna(1.0).to_numpy()
        for rr_name, col in (("none", None), ("ridge_qcexcl", "rr_ridge_qcexcl")):
            exp = tx_expected(ev, cell_rr(ev, col))
            rows.append(
                {
                    "plateau": fit,
                    "rr": rr_name,
                    "core_noise": ev.noise(exp, ev.core),
                    "outcome_free_noise": ev.noise(exp, outcome_free),
                    "obs_over_exp": ev.obs[ev.core].sum() / exp[ev.core].sum(),
                }
            )
    ev.exp_cell = exp_cell
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    """Refit the plateau models and write the tables."""
    g = annotate_groups(export_groups(args.training_ht))
    fits = fit_versions(g, hl.eval(hl.experimental.read_expression(args.plateau_he)))
    g.to_csv(f"{args.out_prefix}.groups.tsv", sep="\t", index=False)
    fits.to_csv(f"{args.out_prefix}.fits.tsv", sep="\t", index=False)
    ctx = g.groupby("ctx")[["observed", "possible", *[f"exp_{f}" for f in FITS]]].sum()
    for f in FITS:
        ctx[f"oe_{f}"] = ctx.observed / ctx[f"exp_{f}"]
    ctx.reset_index().to_csv(f"{args.out_prefix}.contexts.tsv", sep="\t", index=False)

    pd.set_option("display.width", 200)
    logger.info("Lines:\n%s", fits.to_string(index=False))
    for f in FITS:
        tot = g.observed.sum() / g[f"exp_{f}"].sum()
        oe = ctx[f"oe_{f}"] / tot
        gr = g.observed / g[f"exp_{f}"] / tot
        logger.info(
            "%s: total obs/exp %.4f; contexts %.3f-%.3f, SD %.4f; groups SD %.4f",
            f,
            tot,
            ctx[f"oe_{f}"].min(),
            ctx[f"oe_{f}"].max(),
            np.sqrt(np.average((oe - 1) ** 2, weights=ctx.exp_stored)),
            np.sqrt(np.average((gr - 1) ** 2, weights=g.exp_stored)),
        )
    if args.cells:
        noise = downstream_noise(args, ctx)
        noise.to_csv(f"{args.out_prefix}.noise.tsv", sep="\t", index=False)
        logger.info("Per-gene noise:\n%s", noise.round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--training-ht", default="data/plateau/gnomad.v4.1.1.constraint_training.ht"
    )
    parser.add_argument("--plateau-he", default="data/plateau/gnomad.v4.1.1.plateau.he")
    parser.add_argument("--cells", help="Cells TSV with rr_ridge_qcexcl (optional).")
    parser.add_argument("--windows")
    parser.add_argument("--transcripts")
    parser.add_argument("--release-tsv")
    parser.add_argument("--out-prefix", default="plateau_models")
    main(parser.parse_args())
