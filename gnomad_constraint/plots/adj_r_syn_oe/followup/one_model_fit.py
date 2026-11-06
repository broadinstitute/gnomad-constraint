"""Try a single variant-level model of expected counts in place of the stacked models.

For the notes ``notes_one_model_for_expected.md``, linked from the rr report.

Inputs: the synonymous variant cells from ``cluster/export_model_cells.py`` (MANE
Select synonymous SNVs on the autosomes, grouped by transcript, 1-kb window, context,
ref, alt, methylation level and exome AN %), its per-window mapping-quality summary,
Ruchit's window features, the full run's rr cells (for the current stack with the
ridge rr), and the export's transcript table and the release TSV (for the gene sets).

Every fitted model is a binomial model of each possible variant being seen with a
complementary log-log link,

    P(seen) = 1 - exp(-exp(eta)),

fit on cells (all variants in a cell share eta), and scored on chromosomes it wasn't
fit on (odd/even halves). Terms in ``eta``, by model:

- group terms: an intercept and a slope on log mu per substitution type (A>C, A>G,
  A>T, C>A, C>G, C>T outside CpG, CpG transitions); replaces the plateau lines;
- coverage: a linear spline in exome AN % below 90 (flat above); replaces the
  coverage model;
- MQ (optional): the window's mean AS_MQ of called SNVs (spline below 60) and the
  share of them below 50;
- regional (optional): the 52 window features per context class, either fit to DNMs
  first (a Poisson model of DNM counts with per-context intercepts, ridge-penalized)
  and entered with one free strength, or fit directly on gnomAD (ridge-penalized);
- Roulette (optional): the cell's mean log Roulette multiplier as an offset.

Compared with the current stack (stored plateau lines and coverage model) with and
without the full run's ridge rr. Writes ``<prefix>.metrics.tsv`` (noise on the core
and outcome-free sets, outliers, total obs/exp, obs/exp by context class and by
coverage), ``<prefix>.coverage.tsv`` and ``<prefix>.per_transcript.tsv.gz``.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from rr_shrinkage_experiments import chrom_fold
from scipy.optimize import minimize

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

COMP = str.maketrans("ACGT", "TGCA")
SUBS = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "CpG transition"]
CLASSES = ["A/T", "C non-CpG", "CpG"]
COV_KNOTS = [90, 80, 70, 60, 50, 40, 30, 20, 10]
MQ_KNOTS = [60, 55, 50, 40]
COV_BINS = [0, 20, 50, 80, 90, 101]


def load(args: argparse.Namespace) -> tuple:
    """Load the cells, window features and gene sets."""
    c = pd.read_csv(args.cells, sep="\t", compression="gzip")
    c = c[(c.coverage > 0) & (c.possible > 0)].reset_index(drop=True)
    flip = ~c.ref.isin(["A", "C"])
    ref = c.ref.where(~flip, c.ref.str.translate(COMP))
    alt = c.alt.where(~flip, c.alt.str.translate(COMP))
    c["ctx"] = c.context.where(~flip, c.context.str.translate(COMP).str[::-1])
    c["sub"] = np.where(c.cpg, "CpG transition", ref + ">" + alt)
    c["cls"] = np.where(
        c.ctx.str[1] == "A", "A/T", np.where(c.ctx.str[2] == "G", "CpG", "C non-CpG")
    )
    c["contig"] = c.win.str.split("-").str[0]
    c["fold"] = chrom_fold(c.contig)
    if args.test:
        # One chromosome: alternate 1-Mb blocks instead of chromosome halves.
        mb = c.win.str.split("-").str[1].astype(int) // 1_000_000
        c["fold"] = np.where(mb % 2 == 1, "odd", "even")
    c["log_mu"] = np.log(c.mu)

    w = pd.read_csv(args.windows, sep="\t", compression="gzip")
    feats = [x for x in w.columns if x not in ("win", "contig")]
    mq = pd.read_csv(args.window_mq, sep="\t", compression="gzip")
    w = w.merge(mq, on="win", how="outer")
    w = w[w.win.isin(set(c.win))].reset_index(drop=True)
    F = w[feats].to_numpy(float)
    F = np.where(np.isnan(F), np.nanmean(F, axis=0), F)
    F = (F - F.mean(axis=0)) / np.where(F.std(axis=0) > 0, F.std(axis=0), 1)
    win_idx = pd.Series(np.arange(len(w)), index=w.win)
    c["wi"] = win_idx.reindex(c.win).to_numpy()
    c = c[c.wi.notna()].reset_index(drop=True)
    c["wi"] = c.wi.astype(int)

    # Window mapping quality, mean-imputed where the window has no called SNV.
    mean_mq = w.mean_as_mq.fillna(w.mean_as_mq.mean()).to_numpy()
    low_mq = w.frac_mq_below_50.fillna(w.frac_mq_below_50.mean()).to_numpy()
    c["mq"] = mean_mq[c.wi]
    c["low_mq"] = low_mq[c.wi]

    rr = pd.read_csv(
        args.rr_cells, sep="\t", usecols=["win", "context", "rr_ridge_qcexcl"]
    ).drop_duplicates(["win", "context"])
    c = c.merge(rr.rename(columns={"context": "ctx"}), on=["win", "ctx"], how="left")
    c["rr_ridge"] = c.rr_ridge_qcexcl.fillna(1.0)

    tx = pd.read_csv(args.transcripts, sep="\t", compression="gzip").set_index(
        "transcript"
    )
    rel = pd.read_csv(
        args.release_tsv,
        sep="\t",
        usecols=["transcript", "gene_flags"],
        compression="gzip",
        na_values=["NA"],
    ).drop_duplicates("transcript")
    flags = rel.set_index("transcript").gene_flags.astype(str)
    return c, F, tx, flags


def design(c: pd.DataFrame, mq: bool) -> tuple:
    """Get the dense design columns: group intercepts and slopes, coverage, MQ."""
    cols, names = [], []
    lm = c.log_mu.to_numpy() - c.log_mu.mean()
    for s in SUBS:
        m = (c["sub"] == s).to_numpy().astype(float)
        cols += [m, m * lm]
        names += [f"{s}: intercept", f"{s}: log mu"]
    cov = c.coverage.to_numpy(float)
    for k in COV_KNOTS:
        cols.append(np.maximum(0, k - cov) / 10)
        names.append(f"coverage below {k}")
    if mq:
        q = c.mq.to_numpy(float)
        for k in MQ_KNOTS:
            cols.append(np.maximum(0, k - q) / 5)
            names.append(f"window MQ below {k}")
        cols.append(c.low_mq.to_numpy(float))
        names.append("window share MQ < 50")
    return np.column_stack(cols), names


def fit_cloglog(
    n, y, D, offset, F=None, wi=None, ki=None, n_k=0, lam=0.0, start=None
) -> np.ndarray:
    """
    Fit a binomial model with a cloglog link by penalized likelihood (L-BFGS).

    eta = D @ beta + offset (+ F[wi] @ gamma[ki] if F is given); only gamma is penalized.
    """
    p_d = D.shape[1]
    p_g = 0 if F is None else F.shape[1] * n_k

    def unpack(t):
        beta = t[:p_d]
        gamma = None if F is None else t[p_d:].reshape(n_k, F.shape[1])
        return beta, gamma

    def f(t):
        beta, gamma = unpack(t)
        eta = D @ beta + offset
        if gamma is not None:
            eta = eta + (F @ gamma.T)[wi, ki]
        eta = np.clip(eta, -30, 5)
        lam_ = np.exp(eta)
        logp = np.log(-np.expm1(-lam_))
        nll = -(y * logp - (n - y) * lam_).sum()
        g = (n - y) * lam_ - y * lam_ / np.expm1(lam_)
        grad = [D.T @ g]
        if gamma is not None:
            gg = np.stack(
                [
                    F.T @ np.bincount(wi[ki == k], g[ki == k], F.shape[0])
                    for k in range(n_k)
                ]
            )
            nll = nll + 0.5 * lam * (gamma**2).sum()
            grad.append((gg + lam * gamma).ravel())
        return nll, np.concatenate(grad)

    if start is None:
        # Each group's intercept at its mean share seen; every other weight 0.
        t0 = np.zeros(p_d + p_g)
        for j in range(p_d):
            col = D[:, j]
            if set(np.unique(col)) <= {0.0, 1.0} and col.sum() > 0:
                m = col == 1
                share = np.clip(y[m].sum() / n[m].sum(), 1e-6, 1 - 1e-6)
                t0[j] = np.log(-np.log(1 - share))
    else:
        t0 = start
    r = minimize(f, t0, jac=True, method="L-BFGS-B", options={"maxiter": 5000})
    if not r.success:
        logger.warning("cloglog fit: %s", r.message)
    return r.x


def predict_cloglog(t, D, offset, F=None, wi=None, ki=None, n_k=0) -> np.ndarray:
    """Get each cell's probability of a variant being seen."""
    eta = D @ t[: D.shape[1]] + offset
    if F is not None:
        gamma = t[D.shape[1] :].reshape(n_k, F.shape[1])
        eta = eta + (F @ gamma.T)[wi, ki]
    return -np.expm1(-np.exp(np.clip(eta, -30, 5)))


def fit_dnm(c: pd.DataFrame, F: np.ndarray, m: np.ndarray, lam: float) -> np.ndarray:
    """
    Fit DNM counts: Poisson with E = possible * mu * exp(a_ctx + F[w] @ gamma_class).

    Returns the regional score F[w] @ gamma_class for every cell (intercepts left
    out), and the full linear predictor.
    """
    ctxs = sorted(c.ctx.unique())
    ci = pd.Categorical(c.ctx, categories=ctxs).codes
    ki = pd.Categorical(c.cls, categories=CLASSES).codes
    wi = c.wi.to_numpy()
    off = np.log(c.possible.to_numpy() * c.mu.to_numpy())
    d = c.dnm.to_numpy(float)
    n_c, n_k, p = len(ctxs), len(CLASSES), F.shape[1]

    def f(t):
        a, gamma = t[:n_c], t[n_c:].reshape(n_k, p)
        eta = off + a[ci] + (F @ gamma.T)[wi, ki]
        e = np.exp(eta)
        g = (e - d) * m
        nll = ((e - d * eta) * m).sum() + 0.5 * lam * (gamma**2).sum()
        ga = np.bincount(ci, g, n_c)
        gg = np.stack(
            [F.T @ np.bincount(wi[ki == k], g[ki == k], F.shape[0]) for k in range(n_k)]
        )
        return nll, np.concatenate([ga, (gg + lam * gamma).ravel()])

    t0 = np.zeros(n_c + n_k * p)
    t0[:n_c] = np.log(
        np.bincount(ci, d * m, n_c).clip(0.5) / np.bincount(ci, np.exp(off) * m, n_c)
    )
    r = minimize(f, t0, jac=True, method="L-BFGS-B", options={"maxiter": 3000})
    a, gamma = r.x[:n_c], r.x[n_c:].reshape(n_k, p)
    reg = (F @ gamma.T)[wi, ki]
    return reg, off + a[ci] + reg


def score(
    name: str, c: pd.DataFrame, e: np.ndarray, tx: pd.DataFrame, flags: pd.Series
) -> tuple:
    """Get the report's metrics for one model's expected counts."""
    t = (
        c.assign(e=e)
        .groupby("transcript")[["observed", "e", "e_old"]]
        .sum()
        .join(tx[["unflagged"]], how="left")
    )
    unflagged = t.unflagged.astype(str).str.lower().eq("true")
    noflag = (
        flags.reindex(t.index)
        .fillna("[]")
        .replace("nan", "[]")
        .str.strip("[]")
        .str.replace('"', "")
        .eq("")
    )

    def noise(m):
        o, ee = t.observed[m], t.e[m]
        ee = ee * o.sum() / ee.sum()
        return float(((o - ee) ** 2 / ee).mean())

    core = unflagged & (t.e_old >= 20)
    free = noflag & (t.e_old >= 20)
    z = (t.observed - t.e) / np.sqrt(t.e)
    row = {
        "model": name,
        "core_noise": noise(core),
        "outcome_free_noise": noise(free),
        "outliers": int((z.abs() > 8).sum()),
        "outliers_no_qc_flag": int(((z.abs() > 8) & noflag).sum()),
        "obs_over_exp": c.observed.sum() / e.sum(),
    }
    for cls in CLASSES:
        m = (c.cls == cls).to_numpy()
        row[f"oe_{cls}"] = c.observed[m].sum() / e[m].sum()
    cov = pd.cut(c.coverage, COV_BINS, right=False)
    byc = (
        c.assign(e=e, bin=cov.astype(str))
        .groupby("bin")[["observed", "e"]]
        .sum()
        .assign(oe=lambda d: d.observed / d.e)
        .reset_index()
        .assign(model=name)
    )
    return row, byc, t.e.rename(name)


def main(args: argparse.Namespace) -> None:
    """Fit and score every model."""
    c, F, tx, flags = load(args)
    logger.info(
        "%d cells, %d windows, %d transcripts.",
        len(c),
        F.shape[0],
        c.transcript.nunique(),
    )
    n = c.possible.to_numpy(float)
    y = c.observed.to_numpy(float)
    wi = c.wi.to_numpy()
    ki = pd.Categorical(c.cls, categories=CLASSES).codes
    folds = c.fold.to_numpy()
    rou = np.log(c.m_sum.to_numpy() / n)

    # The regional score from DNMs, out of fold, with the penalty that predicts
    # held-out DNMs best.
    d = c.dnm.to_numpy(float)
    best = None
    for lam in args.dnm_lambdas:
        reg, nll = np.zeros(len(c)), 0.0
        for f in ("odd", "even"):
            train = folds != f
            r_, eta = fit_dnm(c, F, train.astype(float), lam)
            reg[~train] = r_[~train]
            nll += float((np.exp(eta[~train]) - d[~train] * eta[~train]).sum())
        logger.info("DNM penalty %g: held-out Poisson NLL %.2f", lam, nll)
        if best is None or nll < best[0]:
            best = (nll, lam, reg)
    _, lam_dnm, reg_dnm = best
    logger.info("DNM penalty %g chosen; regional score sd %.3f", lam_dnm, reg_dnm.std())

    models = {
        "current stack, no rr": c.e_old.to_numpy(),
        "current stack + ridge rr": (c.e_old * c.rr_ridge).to_numpy(),
        "current stack x Roulette": c.e_old.to_numpy() * np.exp(rou),
        "current stack x Roulette + ridge rr": (c.e_old * c.rr_ridge).to_numpy()
        * np.exp(rou),
    }
    specs = [
        ("one model: groups + coverage", False, None, False),
        ("one model: + window MQ", True, None, False),
        ("one model: + DNM regional", False, "dnm", False),
        ("one model: + DNM regional + window MQ", True, "dnm", False),
        ("one model: + gnomAD regional + window MQ", True, "gnomad", False),
        ("one model: + window MQ + Roulette", True, None, True),
        ("one model: + DNM regional + window MQ + Roulette", True, "dnm", True),
    ]
    for name, use_mq, regional, use_rou in specs:
        D, names = design(c, use_mq)
        if regional == "dnm":
            D = np.column_stack([D, reg_dnm])
        offset = np.zeros(len(c))
        if use_rou:
            offset = offset + rou
        e = np.zeros(len(c))
        for f in ("odd", "even"):
            tr, te = folds != f, folds == f
            if regional == "gnomad":
                t = fit_cloglog(
                    n[tr],
                    y[tr],
                    D[tr],
                    offset[tr],
                    F,
                    wi[tr],
                    ki[tr],
                    len(CLASSES),
                    lam=args.gnomad_lambda,
                )
                p = predict_cloglog(
                    t, D[te], offset[te], F, wi[te], ki[te], len(CLASSES)
                )
            else:
                t = fit_cloglog(n[tr], y[tr], D[tr], offset[tr])
                p = predict_cloglog(t, D[te], offset[te])
            e[te] = n[te] * p
            if regional == "dnm":
                logger.info(
                    "%s (%s fold): DNM regional strength %.3f",
                    name,
                    f,
                    t[D.shape[1] - 1],
                )
        models[name] = e
        logger.info("Fitted %s", name)

    rows, covs, per_tx = [], [], []
    for name, e in models.items():
        r, b, t = score(name, c, e, tx, flags)
        rows.append(r)
        covs.append(b)
        per_tx.append(t)
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}.metrics.tsv", sep="\t", index=False)
    pd.concat(covs).to_csv(f"{args.out_prefix}.coverage.tsv", sep="\t", index=False)
    pd.concat(per_tx, axis=1).reset_index().to_csv(
        f"{args.out_prefix}.per_transcript.tsv.gz", sep="\t", index=False
    )
    pd.set_option("display.width", 250)
    logger.info("\n%s", out.round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", default="data/one_model/syn_model_cells.tsv.bgz")
    parser.add_argument(
        "--window-mq", default="data/one_model/syn_model_cells.window_mq.tsv.bgz"
    )
    parser.add_argument("--windows", default="data/rr_eval_export/windows.tsv.bgz")
    parser.add_argument("--rr-cells", default="data/full_run/syn_cells_rr_full_run.tsv")
    parser.add_argument(
        "--transcripts", default="data/rr_eval_export/transcripts.tsv.bgz"
    )
    parser.add_argument(
        "--release-tsv", default="../gnomad.v4.1.2.constraint_metrics.tsv.bgz"
    )
    parser.add_argument(
        "--dnm-lambdas",
        type=float,
        nargs="+",
        default=[300, 1000, 3000, 10000, 30000, 100000],
    )
    parser.add_argument("--gnomad-lambda", type=float, default=1000.0)
    parser.add_argument("--out-prefix", default="one_model")
    parser.add_argument(
        "--test",
        action="store_true",
        help="Fold by alternate 1-Mb blocks (one chromosome).",
    )
    main(parser.parse_args())
