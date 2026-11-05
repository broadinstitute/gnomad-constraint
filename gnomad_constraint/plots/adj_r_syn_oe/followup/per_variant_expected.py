"""Numbers for the companion note on per-variant expected counts and plateau curves.

Inputs: ``plateau_models.groups.tsv`` (the v4.1.1 high-coverage synonymous training
sites per mutation-rate group, from ``plateau_models.py``), and optionally Roulette's
rates for high-quality synonymous variants with the hg38 reference.

1. **Which fits make total expected equal total observed.** On the training groups,
   for each of the two plateau groups (CpG transitions and everything else), fits
   share seen = f(rate) five ways: a straight line by unweighted least squares (the
   stored models), the same weighted by possible variants, a saturating curve
   1 - exp(-(a + b * rate)) by unweighted least squares on the shares and by binomial
   likelihood, and a logistic curve (logit share = a + b * rate) by binomial likelihood.
   Writes each fit's parameters and total obs/exp to ``<prefix>.totals.tsv``.
2. **What applying a curve to an averaged rate would cost.** If the variants of a group
   had different rates (as Roulette's per-variant rates do within a trinucleotide
   change), a curve applied to the group's mean rate gives a different expected count
   from the curve applied to each variant and summed (Jensen's inequality). With
   ``--roulette``, takes each Roulette variant's rate relative to its (strand-collapsed
   trinucleotide context, alt) mean as the within-group spread, puts it around each
   training group's rate, and writes the relative error of the averaged-rate shortcut
   under the saturating curve to ``<prefix>.jensen.tsv``. A straight line has no
   such error. CpG transitions are only approximate here: Roulette's spread there
   includes methylation, which our groups already split into three levels.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from roulette_combo import BASE_CODE, LETTERS, read_chromosome
from scipy.optimize import minimize

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

SCALE = 1e8
"""Rates are fit in units of 1e-8 per generation."""


def saturating(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Share seen under the saturating curve."""
    return 1 - np.exp(-(a + b * x))


def logistic(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Share seen under the logistic curve."""
    return 1 / (1 + np.exp(-(a + b * x)))


def binomial_fit(f, x: np.ndarray, k: np.ndarray, n: np.ndarray, start: list) -> tuple:
    """Fit a two-parameter share curve by binomial likelihood."""

    def nll(t: np.ndarray) -> float:
        p = np.clip(f(x, *t), 1e-12, 1 - 1e-12)
        return float(-(k * np.log(p) + (n - k) * np.log(1 - p)).sum())

    r = minimize(
        nll,
        start,
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-10, "maxiter": 20000},
    )
    return tuple(float(v) for v in r.x)


def ls_fit(f, x: np.ndarray, y: np.ndarray, start: list) -> tuple:
    """Fit a two-parameter share curve by unweighted least squares on the shares."""
    r = minimize(
        lambda t: float(((y - f(x, *t)) ** 2).sum()),
        start,
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000},
    )
    return tuple(float(v) for v in r.x)


def fit_all(g: pd.DataFrame) -> pd.DataFrame:
    """Fit every form to each plateau group and get its total obs/exp."""
    rows = []
    for cpg, d in g.groupby("cpg"):
        x = d.mu.to_numpy() * SCALE
        k, n = d.observed.to_numpy(), d.possible.to_numpy().astype(float)
        y = k / n
        line = "CpG transition" if cpg else "everything else"
        lin = lambda x_, a, b: a + b * x_  # noqa: E731
        stored = d.exp_stored.sum()
        X = np.column_stack([np.ones_like(x), x])
        wls = np.linalg.solve(X.T @ (X * n[:, None]), X.T @ (y * n))
        sat_start = [0.03, 0.4] if not cpg else [0.4, 0.3]
        fits = {
            "line, unweighted (stored)": (lin, None, stored),
            "line, weighted": (lin, tuple(wls), None),
            "saturating, least squares": (
                saturating,
                ls_fit(saturating, x, y, sat_start),
                None,
            ),
            "saturating, likelihood": (
                saturating,
                binomial_fit(saturating, x, k, n, sat_start),
                None,
            ),
            "logistic, likelihood": (
                logistic,
                binomial_fit(logistic, x, k, n, [-2.0, 1.0] if not cpg else [0.0, 0.3]),
                None,
            ),
        }
        for name, (f, params, exp_total) in fits.items():
            if exp_total is None:
                exp_total = float((n * f(x, *params)).sum())
                a, b = params
            else:
                a, b = np.nan, np.nan
            rows.append(
                {
                    "line": line,
                    "fit": name,
                    "a": a,
                    "b": b,
                    "observed": float(k.sum()),
                    "expected": exp_total,
                    "obs_over_exp": float(k.sum()) / exp_total,
                }
            )
    out = pd.DataFrame(rows)
    tot = out.groupby("fit", sort=False)[["observed", "expected"]].sum().reset_index()
    tot["line"] = "both"
    tot["obs_over_exp"] = tot.observed / tot.expected
    return pd.concat([out, tot], ignore_index=True)


def roulette_spread(path: str, twobit: str, n_per_group: int, seed: int) -> dict:
    """Get a sample of each Roulette variant's rate relative to its group's mean."""
    d = pd.read_csv(
        path,
        sep="\t",
        usecols=["CHROM", "POS", "REF", "ALT", "mu"],
        dtype={"CHROM": str},
    )
    d = d[d.REF.str.len().eq(1) & d.ALT.str.len().eq(1)]
    parts = []
    for chrom, c in d.groupby("CHROM", sort=False):
        seq = read_chromosome(twobit, f"chr{chrom}")
        p = c.POS.to_numpy()
        codes = BASE_CODE[seq[(p - 1)[:, None] + np.arange(-1, 2)]]
        alt = BASE_CODE[np.frombuffer("".join(c.ALT).encode(), dtype=np.uint8)]
        ok = (codes < 4).all(axis=1)
        codes, alt, mu = codes[ok], alt[ok], c.mu.to_numpy()[ok]
        flip = codes[:, 1] >= 2
        codes[flip] = (3 - codes[flip])[:, ::-1]
        alt[flip] = 3 - alt[flip]
        ctx = LETTERS[codes].astype(np.uint8).copy().view("S3").ravel()
        parts.append(
            pd.DataFrame(
                {"ctx": ctx, "alt": LETTERS[alt].astype(np.uint8).view("S1"), "mu": mu}
            )
        )
    r = pd.concat(parts)
    r["ctx"] = r.ctx.str.decode("ascii")
    r["alt"] = r.alt.str.decode("ascii")
    rng = np.random.default_rng(seed)
    out = {}
    for (ctx, alt), x in r.groupby(["ctx", "alt"]).mu:
        x = x.to_numpy()
        x = x[x > 0]
        if len(x) > n_per_group:
            x = rng.choice(x, n_per_group, replace=False)
        out[(ctx, alt)] = x / x.mean()
    return out


def jensen(g: pd.DataFrame, totals: pd.DataFrame, spread: dict) -> pd.DataFrame:
    """Get the error of applying the saturating curve to each group's mean rate."""
    comp = str.maketrans("ACGT", "TGCA")
    rows = []
    for _, r in g.iterrows():
        flip = r.ref not in "AC"
        alt = r.alt.translate(comp) if flip else r.alt
        rel = spread.get((r.ctx, alt))
        if rel is None:
            continue
        line = "CpG transition" if r.cpg else "everything else"
        a, b = (
            totals[(totals.line == line) & (totals.fit == "saturating, likelihood")]
            .iloc[0][["a", "b"]]
            .astype(float)
        )
        x = r.mu * SCALE
        per_variant = saturating(x * rel, a, b).mean()
        averaged = saturating(np.array([x]), a, b)[0]
        rows.append(
            {
                "ctx": r.ctx,
                "alt": alt,
                "methylation_level": r.methylation_level,
                "cls": r.cls,
                "possible": r.possible,
                "rate_cv": float(rel.std()),
                "averaged_over_per_variant": averaged / per_variant,
            }
        )
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    """Write the tables."""
    g = pd.read_csv(args.groups, sep="\t")
    totals = fit_all(g)
    totals.to_csv(f"{args.out_prefix}.totals.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    logger.info("Totals:\n%s", totals.round(5).to_string(index=False))
    if args.roulette:
        spread = roulette_spread(
            args.roulette, args.twobit, args.n_per_group, args.seed
        )
        j = jensen(g, totals, spread)
        j.to_csv(f"{args.out_prefix}.jensen.tsv", sep="\t", index=False)
        summ = pd.DataFrame(
            {
                cls: {
                    "groups": len(d),
                    "median_cv": d.rate_cv.median(),
                    "weighted_error": np.average(
                        d.averaged_over_per_variant - 1, weights=d.possible
                    ),
                    "max_error": (d.averaged_over_per_variant - 1).abs().max(),
                }
                for cls, d in j.groupby("cls")
            }
        ).T
        logger.info(
            "Averaged-rate error under the saturating curve:\n%s", summ.round(4)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", default="plateau_models.groups.tsv")
    parser.add_argument(
        "--roulette", help="Roulette all_hq_synonymous_variants.tsv.gz (optional)."
    )
    parser.add_argument("--twobit", default="data/ext_hg38/hg38.2bit")
    parser.add_argument("--n-per-group", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-prefix", default="per_variant_expected")
    main(parser.parse_args())
