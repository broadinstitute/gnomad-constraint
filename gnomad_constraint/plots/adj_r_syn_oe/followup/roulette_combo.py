"""Try combining Roulette's per-site mutation rates with the DNM-trained rr, locally.

Roulette (Seplyarskiy et al. 2023) estimates a mutation rate for every possible SNV
from gnomAD rare variants, using the pentanucleotide context, wider sequence context,
methylation, transcription and an empirical regional mutation rate. Its rates for all
high-quality synonymous variants are in ``all_hq_synonymous_variants.tsv.gz``.

For each (1kb window, trinucleotide context) cell of Ruchit's evaluation export, the
Roulette multiplier is the cell's summed Roulette rate divided by the sum expected from
Roulette's own average rate for each (trinucleotide context, alt): how much more or less
mutable Roulette thinks the cell's synonymous sites are than an average site with the
same trinucleotide change. It splits into two parts through Roulette's average rate for
each (heptamer context, alt):

- a sequence-context part (heptamer average over trinucleotide average): what the three
  bases on either side add, which our trinucleotide mutation rates leave out;
- a regional part (Roulette rate over heptamer average): everything else Roulette sees
  at the site, mostly its empirical regional rate, plus methylation and transcription.

Versions scored on the synonymous outcome, like the rest of the report:

- the ridge rr alone (DNMs against the current background, penalty 10 per DNM);
- Roulette alone, all contexts and outside CpG contexts only (our expected counts
  already model CpG methylation), and each of its two parts outside CpG;
- ridge × Roulette (or × one of its parts) outside CpG;
- ridge with Roulette as a 53rd feature: log of the multiplier, so the DNMs decide how
  much of Roulette to use on top of the regional features;
- ridge on a Roulette offset: Roulette fixed at full strength outside CpG, and the ridge
  fitted to the DNMs on top of it, so it only learns what Roulette misses.

Each version is also checked against the DNMs alone, which never touch gnomAD: a
Poisson model of DNM counts per cell with all possible sites as the exposure and one
slope on log rr per context class. A slope of 1 means the rr is right about where DNMs
land; the fitted rr versions use the models from the other chromosome fold.

The DNMs and backgrounds are the local approximation of ``rr_background_experiments.py``.
"""

import argparse
import logging
import struct

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rr_background_experiments import (
    cell_background,
    class_of,
    class_slopes,
    fit_weighted,
    load_dnms,
)
from rr_shrinkage_experiments import CONTEXTS, FOLDS, Evaluator, allowed_columns
from scipy.optimize import minimize, minimize_scalar

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

BACKGROUND = "absent"
PENALTY = 10.0
CODE = np.frombuffer(b"TCAG", dtype=np.uint8)
LETTERS = np.frombuffer(b"ACGT", dtype=np.uint8)
BASE_CODE = np.full(256, 4, dtype=np.int64)
BASE_CODE[LETTERS] = np.arange(4)
FLANK = 3
CLASSES = ["A/T", "C non-CpG", "CpG"]


def read_chromosome(path: str, name: str) -> np.ndarray:
    """
    Read one sequence of a .2bit file as an array of ASCII bases (N blocks as N).

    :param path: .2bit path.
    :param name: Sequence name, e.g. "chr1".
    :return: uint8 array of bases.
    """
    with open(path, "rb") as f:
        _, _, n_seq, _ = struct.unpack("<IIII", f.read(16))
        offsets = {}
        for _ in range(n_seq):
            nm = f.read(ord(f.read(1))).decode()
            offsets[nm] = struct.unpack("<I", f.read(4))[0]
        f.seek(offsets[name])
        size, n_n = struct.unpack("<II", f.read(8))
        n_starts = np.frombuffer(f.read(4 * n_n), dtype="<u4")
        n_sizes = np.frombuffer(f.read(4 * n_n), dtype="<u4")
        (n_mask,) = struct.unpack("<I", f.read(4))
        f.seek(8 * n_mask + 4, 1)
        packed = np.frombuffer(f.read((size + 3) // 4), dtype=np.uint8)
    codes = np.empty(len(packed) * 4, dtype=np.uint8)
    for i, shift in enumerate((6, 4, 2, 0)):
        codes[i::4] = (packed >> shift) & 3
    seq = CODE[codes[:size]].copy()
    for s, n in zip(n_starts, n_sizes):
        seq[s : s + n] = ord("N")
    return seq


def roulette_variants(path: str, twobit: str, keys: pd.DataFrame) -> pd.DataFrame:
    """
    Read Roulette's synonymous variants with strand-collapsed contexts.

    :param path: Roulette ``all_hq_synonymous_variants.tsv.gz``.
    :param twobit: hg38.2bit path.
    :param keys: (``win``, ``context``) cells to keep.
    :return: Rows with ``win``, ``context``, ``alt``, ``ext`` (heptamer and alt as an
        integer) and ``mu``.
    """
    d = pd.read_csv(
        path,
        sep="\t",
        usecols=["CHROM", "POS", "REF", "ALT", "mu"],
        dtype={"CHROM": str},
    )
    d = d[d.REF.str.len().eq(1) & d.ALT.str.len().eq(1)]
    out = []
    for chrom, g in d.groupby("CHROM", sort=False):
        name = f"chr{chrom}"
        seq = read_chromosome(twobit, name)
        pos = g.POS.to_numpy()
        codes = BASE_CODE[seq[(pos - 1)[:, None] + np.arange(-FLANK, FLANK + 1)]]
        ref = BASE_CODE[np.frombuffer("".join(g.REF).encode(), dtype=np.uint8)]
        alt = BASE_CODE[np.frombuffer("".join(g.ALT).encode(), dtype=np.uint8)]
        ok = (codes < 4).all(axis=1) & (codes[:, FLANK] == ref)
        codes, alt, g = codes[ok], alt[ok], g[ok]
        # G and T references go to the other strand, so the middle base is A or C.
        flip = codes[:, FLANK] >= 2
        codes[flip] = (3 - codes[flip])[:, ::-1]
        alt[flip] = 3 - alt[flip]
        ext = (codes * 4 ** np.arange(2 * FLANK, -1, -1)).sum(axis=1) * 4 + alt
        tri = LETTERS[codes[:, FLANK - 1 : FLANK + 2]].astype(np.uint8).copy()
        start = g.POS.to_numpy() // 1000 * 1000
        out.append(
            pd.DataFrame(
                {
                    "win": (
                        name
                        + "-"
                        + pd.Series(start).astype(str)
                        + "-"
                        + pd.Series(start + 1000).astype(str)
                    ).to_numpy(),
                    "context": tri.view("S3").ravel().astype(str),
                    "alt": alt,
                    "ext": ext,
                    "mu": g.mu.to_numpy(),
                }
            )
        )
        logger.info("%s: %d Roulette synonymous variants", name, ok.sum())
    return pd.concat(out).merge(keys, on=["win", "context"])


def roulette_cells(r: pd.DataFrame) -> pd.DataFrame:
    """
    Get the Roulette multiplier of each (window, context) cell and its two parts.

    :param r: Variants from ``roulette_variants``.
    :return: Rows with ``win``, ``context``, ``multiplier`` (Roulette over its
        trinucleotide average), ``sequence_part`` (heptamer over trinucleotide average)
        and ``regional_part`` (Roulette over heptamer average).
    """
    r = r.assign(
        tri_mean=r.groupby(["context", "alt"]).mu.transform("mean"),
        ext_mean=r.groupby("ext").mu.transform("mean"),
    )
    c = (
        r.groupby(["win", "context"])
        .agg(
            roulette_sum=("mu", "sum"),
            tri_sum=("tri_mean", "sum"),
            ext_sum=("ext_mean", "sum"),
            n_variants=("mu", "size"),
        )
        .reset_index()
    )
    c["multiplier"] = c.roulette_sum / c.tri_sum
    c["sequence_part"] = c.ext_sum / c.tri_sum
    c["regional_part"] = c.roulette_sum / c.ext_sum
    return c


def matrix_from_cells(ev: Evaluator, cells: pd.DataFrame, col: str) -> np.ndarray:
    """
    Put a per-(window, context) value into a (window x context) matrix (NaN = none).

    :param ev: Evaluator with the windows.
    :param cells: Rows with ``win``, ``context`` and ``col``.
    :param col: Column to place.
    :return: Matrix.
    """
    m = np.full((len(ev.windows), len(CONTEXTS)), np.nan)
    wi = (
        pd.Series(np.arange(len(ev.windows)), index=ev.windows.win)
        .reindex(cells.win)
        .to_numpy()
    )
    ci = cells.context.map({c: i for i, c in enumerate(CONTEXTS)}).to_numpy()
    ok = ~np.isnan(wi)
    m[wi[ok].astype(int), ci[ok]] = cells[col].to_numpy()[ok]
    return m


def fit_offset(
    x: np.ndarray, y: np.ndarray, w: np.ndarray, offset: np.ndarray, penalty: float
) -> dict:
    """
    Fit a ridge logistic model with a fixed offset, penalized like ``fit_weighted``.

    Minimizes the weighted log loss plus (n_DNM * penalty / 2) * ||coef||^2 on z-scored
    features, which is ``fit_weighted``'s objective with an offset added to the logit.

    :param x: Features.
    :param y: Outcome (1 = DNM).
    :param w: Row weights.
    :param offset: Fixed logit offset per row.
    :param penalty: Per-DNM penalty.
    :return: Model with ``mean``, ``sd``, ``coef`` and ``intercept``.
    """
    mean = np.average(x, axis=0, weights=w)
    sd = np.sqrt(np.average((x - mean) ** 2, axis=0, weights=w))
    sd[sd == 0] = 1.0
    z = (x - mean) / sd
    lam = y.sum() * penalty

    def objective(theta: np.ndarray) -> tuple:
        eta = theta[0] + offset + z @ theta[1:]
        loss = np.sum(w * (np.logaddexp(0, eta) - y * eta))
        resid = w * (1 / (1 + np.exp(-eta)) - y)
        grad = np.concatenate([[resid.sum()], z.T @ resid + lam * theta[1:]])
        return loss + lam / 2 * theta[1:] @ theta[1:], grad

    p0 = np.average(y, weights=w)
    theta0 = np.concatenate([[np.log(p0 / (1 - p0))], np.zeros(z.shape[1])])
    res = minimize(objective, theta0, jac=True, method="L-BFGS-B")
    return {"mean": mean, "sd": sd, "coef": res.x[1:], "intercept": res.x[0]}


def fit_context(
    context: str,
    dnm: pd.DataFrame,
    bg: pd.DataFrame,
    cols: list,
    offset_col: str = None,
) -> tuple:
    """
    Fit one context's ridge models (per DNM penalty) for each fold.

    :param context: Trinucleotide context.
    :param dnm: DNM rows of the context.
    :param bg: Background rows of the context with the weight column.
    :param cols: Feature columns.
    :param offset_col: Column holding a fixed logit offset, or None for no offset.
    :return: (context, {fold: model}).
    """
    d = dnm.dropna(subset=cols)
    b = bg[bg[BACKGROUND] > 0].dropna(subset=cols)
    x = np.vstack([d[cols].to_numpy(), b[cols].to_numpy()])
    y = np.concatenate([np.ones(len(d)), np.zeros(len(b))])
    w = np.concatenate([np.ones(len(d)), b[BACKGROUND].to_numpy()])
    fold = np.concatenate([d.fold.to_numpy(), b.fold.to_numpy()])
    if offset_col:
        off = np.concatenate([d[offset_col].to_numpy(), b[offset_col].to_numpy()])
    out = {}
    for f in (*FOLDS, "all"):
        k = fold != f
        if offset_col:
            m = fit_offset(x[k], y[k], w[k], off[k], PENALTY)
        else:
            m = fit_weighted(x[k], y[k], w[k], PENALTY)
        m["cols"] = cols
        out[f] = m
    return context, out


def predict(ev: Evaluator, fits: dict, extra: np.ndarray, oof: bool) -> np.ndarray:
    """
    Build an rr matrix from per-context fits, with an optional per-cell extra feature.

    :param ev: Evaluator with the windows.
    :param fits: Context -> {fold: model}.
    :param extra: (window x context) values of the ``roulette`` feature, or None.
    :param oof: Whether to use the other-fold models.
    :return: rr matrix.
    """
    rr = np.full((len(ev.windows), len(CONTEXTS)), np.nan)
    w = ev.windows
    for context, f in fits.items():
        j = CONTEXTS.index(context)
        for fold in FOLDS if oof else ("all",):
            rows = np.where(w.fold == fold)[0] if oof else np.arange(len(w))
            m = f[fold]
            cols = [c for c in m["cols"] if c != "roulette"]
            x = w.iloc[rows][cols].to_numpy()
            if "roulette" in m["cols"]:
                x = np.column_stack([x, np.nan_to_num(extra[rows, j], nan=0.0)])
            z = (x - m["mean"]) / m["sd"]
            eta0 = m["intercept"]
            rr[rows, j] = (1 / (1 + np.exp(-(eta0 + z @ m["coef"])))) / (
                1 / (1 + np.exp(-eta0))
            )
    return rr


def fit_all(
    dnm: pd.DataFrame,
    agg: pd.DataFrame,
    extra_feature: bool,
    offset_col: str,
    n_jobs: int,
) -> dict:
    """
    Fit every context's ridge models.

    :param dnm: DNM rows with the features.
    :param agg: Background cells with the features.
    :param extra_feature: Whether to add the ``roulette`` feature.
    :param offset_col: Offset column, or None.
    :param n_jobs: Parallel jobs.
    :return: Context -> {fold: model}.
    """
    jobs = []
    for c in CONTEXTS:
        cols = allowed_columns(c) + (["roulette"] if extra_feature else [])
        keep = ["fold", *cols] + ([offset_col] if offset_col else [])
        jobs.append(
            delayed(fit_context)(
                c,
                dnm[dnm.context == c][keep],
                agg[agg.context == c][[BACKGROUND, *keep]],
                cols,
                offset_col,
            )
        )
    return dict(Parallel(n_jobs=n_jobs)(jobs))


def dnm_calibration(agg: pd.DataFrame, log_rr: np.ndarray) -> dict:
    """
    Fit DNM counts per cell to log rr: one slope per context class.

    Poisson model: DNMs in a cell ~ Poisson(n * exp(a_context + b * log rr)), with n the
    possible synonymous variants of the cell. The context intercepts are profiled out,
    so the log likelihood in b is b * sum(d * l) - sum_c D_c * log(sum_c n * exp(b * l)).

    :param agg: Cells with ``context``, ``n`` and ``d`` (DNM count).
    :param log_rr: log rr per cell (0 = none).
    :return: {class: slope, SE, and log-likelihood gain per 1,000 DNMs at slope 1 and
        at the fitted slope}.
    """
    out = {}
    cls = class_of(agg.context)
    for c in CLASSES:
        m = cls == c
        d, n, l = agg.d.to_numpy()[m], agg.n.to_numpy()[m], log_rr[m]
        ctx = pd.factorize(agg.context.to_numpy()[m])[0]
        dc = np.bincount(ctx, weights=d)
        s1 = np.sum(d * l)

        def ll(b: float) -> float:
            return b * s1 - np.sum(
                dc * np.log(np.bincount(ctx, weights=n * np.exp(b * l)))
            )

        b = minimize_scalar(lambda v: -ll(v), bounds=(-2, 5), method="bounded").x
        h = 1e-3
        curv = (ll(b + h) - 2 * ll(b) + ll(b - h)) / h**2
        out[c] = {
            "dnm_slope": b,
            "dnm_slope_se": float(np.sqrt(-1 / curv)) if curv < 0 else np.nan,
            "dnm_gain_at_1": 1000 * (ll(1) - ll(0)) / dc.sum(),
            "dnm_gain_at_fit": 1000 * (ll(b) - ll(0)) / dc.sum(),
        }
    return out


def main(args: argparse.Namespace) -> None:
    """Build the Roulette multiplier, fit the ridge versions, and score everything."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    feat_cols = sorted({c for ctx in CONTEXTS for c in allowed_columns(ctx)})
    agg = cell_background(ev).merge(ev.windows[["win", "fold", *feat_cols]], on="win")
    rv = roulette_variants(args.roulette, args.twobit, agg[["win", "context"]])
    rcells = roulette_cells(rv)
    del rv
    rcells.to_csv(f"{args.out_prefix}.cells.tsv.gz", sep="\t", index=False)
    for col in ("multiplier", "sequence_part", "regional_part"):
        logger.info(
            "Roulette %s over %d cells: %s",
            col,
            len(rcells),
            rcells[col].describe().round(3).to_dict(),
        )
    mats = {
        col: matrix_from_cells(ev, rcells, col)
        for col in ("multiplier", "sequence_part", "regional_part")
    }
    cpg = [CONTEXTS.index(c) for c in CONTEXTS if c[1] == "C" and c[2] == "G"]

    def non_cpg(m: np.ndarray) -> np.ndarray:
        m = m.copy()
        m[:, cpg] = np.nan
        return m

    dnm = load_dnms(args.dnm_files, args.twobit)
    dnm = dnm.merge(agg[["win", "context"]], on=["win", "context"]).merge(
        ev.windows[["win", "fold", *feat_cols]], on="win"
    )
    agg = agg.merge(
        dnm.groupby(["win", "context"]).size().rename("d").reset_index(),
        on=["win", "context"],
        how="left",
    )
    agg["d"] = agg.d.fillna(0)
    lr = rcells.assign(
        roulette=np.log(rcells.multiplier),
        roulette_offset=np.where(
            class_of(rcells.context) == "CpG", 0.0, np.log(rcells.multiplier)
        ),
    )[["win", "context", "roulette", "roulette_offset"]]
    agg = agg.merge(lr, on=["win", "context"], how="left")
    dnm = dnm.merge(lr, on=["win", "context"], how="left")
    for col in ("roulette", "roulette_offset"):
        agg[col] = agg[col].fillna(0.0)
        dnm[col] = dnm[col].fillna(0.0)

    log_roul = np.log(mats["multiplier"])
    versions = {}
    fits = fit_all(dnm, agg, False, None, args.n_jobs)
    versions["ridge"] = (
        predict(ev, fits, None, False),
        predict(ev, fits, None, True),
    )
    fits = fit_all(dnm, agg, True, None, args.n_jobs)
    versions["ridge_with_roulette_feature"] = (
        predict(ev, fits, log_roul, False),
        predict(ev, fits, log_roul, True),
    )
    coef = pd.DataFrame(
        [
            {
                "context": c,
                "class": class_of(pd.Series([c]))[0],
                "roulette_coef_per_sd": f["all"]["coef"][-1],
                "roulette_sd": f["all"]["sd"][-1],
            }
            for c, f in fits.items()
        ]
    )
    coef["roulette_coef"] = coef.roulette_coef_per_sd / coef.roulette_sd
    coef.to_csv(f"{args.out_prefix}.feature_coef.tsv", sep="\t", index=False)
    logger.info(
        "Roulette feature: coefficient on log multiplier by class: %s",
        coef.groupby("class").roulette_coef.mean().round(2).to_dict(),
    )
    fits = fit_all(dnm, agg, False, "roulette_offset", args.n_jobs)
    roul_off = np.where(np.isnan(non_cpg(mats["multiplier"])), 1.0, mats["multiplier"])
    versions["ridge_on_roulette_offset"] = (
        roul_off * predict(ev, fits, None, False),
        roul_off * predict(ev, fits, None, True),
    )

    versions["roulette"] = (mats["multiplier"],) * 2
    for col, name in (
        ("multiplier", "roulette_non_cpg"),
        ("sequence_part", "roulette_sequence_non_cpg"),
        ("regional_part", "roulette_regional_non_cpg"),
    ):
        m = non_cpg(mats[col])
        versions[name] = (m, m)
        r_all, r_oof = versions["ridge"]
        f = np.where(np.isnan(m), 1.0, m)
        versions[f"ridge_x_{name}"] = (r_all * f, r_oof * f)

    wi = (
        pd.Series(np.arange(len(ev.windows)), index=ev.windows.win)
        .reindex(agg.win)
        .to_numpy()
        .astype(int)
    )
    ci = agg.context.map({c: i for i, c in enumerate(CONTEXTS)}).to_numpy()

    def cell_log_rr(rr: np.ndarray) -> np.ndarray:
        return np.nan_to_num(np.log(rr[wi, ci]), nan=0.0)

    base = ev.score(None)
    rows = [
        {
            "version": "none",
            "core_noise": base["core_noise"],
            "oof_noise": base["core_noise"],
            "outcome_free_noise": ev.noise(ev.exp0, ev.no_qc_flag & (ev.exp0 >= 20)),
            "outliers": base["outliers"],
            "outliers_no_qc_flag": base["outliers_no_qc_flag"],
        }
    ]
    for name, (ra, ro) in versions.items():
        sc = ev.score(ra)
        row = {
            "version": name,
            "core_noise": sc["core_noise"],
            "oof_noise": ev.score(ro)["core_noise"],
            "outcome_free_noise": ev.noise(
                ev.expected(ra), ev.no_qc_flag & (ev.exp0 >= 20)
            ),
            "outliers": sc["outliers"],
            "outliers_no_qc_flag": sc["outliers_no_qc_flag"],
            "slope": sc["slope"],
            **{f"slope_{c}": s for c, s in class_slopes(ev, ra, False).items()},
        }
        for c, v in dnm_calibration(agg, cell_log_rr(ro)).items():
            row.update({f"{k}_{c}": x for k, x in v.items()})
        rows.append(row)
        logger.info(
            "%s",
            {k: round(v, 3) if isinstance(v, float) else v for k, v in row.items()},
        )
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.metrics.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, help="Decompressed syn_cells TSV.")
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--dnm-files", nargs="+", required=True)
    parser.add_argument("--twobit", required=True)
    parser.add_argument("--roulette", required=True)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
