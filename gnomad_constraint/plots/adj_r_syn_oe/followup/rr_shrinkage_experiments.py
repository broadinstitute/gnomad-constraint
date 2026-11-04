"""Try alternatives to a single rr strength (alpha) on Ruchit's rr data, locally.

Inputs are local copies of Ruchit Panchal's outputs under
gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/:

- ``regional_rate{,_v1_ddd_asd}/gnomad.v4.1.1.rr_training.ht`` exported to TSV: rr
  training sites (DNM ``y`` = 1, background ``y`` = 0) with the 52 regional feature
  columns, for the Satterstrom and old 17k DNM sets.
- ``regional_rate{,_v1_ddd_asd}/gnomad.v4.1.1.rr_feature_screen.tsv``: his univariate
  feature screen.
- ``regional_rate/rr_eval_export/``: per (transcript, 1kb window, context) synonymous
  expected and observed counts (``syn_cells``, here the copy with rr columns written by
  ``remap_rr_to_shet.py --annotate-rr``), the window features (``windows``) and the
  transcript table (``transcripts``).

Every rr is built out of fold: models applied to odd-chromosome windows are trained
on even-chromosome sites and vice versa (chrX counts as odd, chrY as even). rr is then
normalized per context over all exported synonymous cells and scored with the same
metrics as Ruchit's report.
"""

import argparse
import logging
import zlib

import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression, PoissonRegressor

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

FEATURES = [
    "dist2telo",
    "dist2cent",
    "GC_content",
    "LCR",
    "SINE",
    "LINE",
    "recomb_male",
    "recomb_female",
    "met_sperm",
    "Nucleosome",
    "cDNM_maternal_05M",
    "cDNM_paternal_05M",
    "CpG_island",
]
SCREEN_ORDER = [
    "dist2telo",
    "dist2cent",
    "LCR",
    "SINE",
    "LINE",
    "GC_content",
    "recomb_male",
    "recomb_female",
    "met_sperm",
    "Nucleosome",
    "CpG_island",
    "cDNM_maternal_05M",
    "cDNM_paternal_05M",
]
SCALES = ["1k", "10k", "100k", "1M"]
COLUMNS = [f"{f}_{s}" for s in SCALES for f in FEATURES]
CPG_CONTEXTS = ["ACG", "CCG", "GCG", "TCG"]
METHYLATION_FEATURES = ["GC_content", "SINE", "met_sperm", "Nucleosome", "CpG_island"]
CONTEXTS = sorted(f"{a}{c}{b}" for a in "ACGT" for c in "AC" for b in "ACGT")
FOLDS = ("odd", "even")
ALPHAS = [0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.25]


def chrom_fold(contig: pd.Series) -> pd.Series:
    """
    Get the chromosome-parity fold of each contig (chrX odd, chrY even).

    :param contig: Contig names, e.g. "chr12".
    :return: "odd" or "even" per contig.
    """
    name = contig.str.replace("chr", "", regex=False)
    num = pd.to_numeric(name.replace({"X": "1", "Y": "2"}), errors="coerce")
    return np.where(num % 2 == 1, "odd", "even")


def allowed_columns(context: str) -> list:
    """
    Get the feature columns a context may use (no methylation features for CpG).

    :param context: Trinucleotide context.
    :return: Feature column names.
    """
    feats = [
        f
        for f in FEATURES
        if not (context in CPG_CONTEXTS and f in METHYLATION_FEATURES)
    ]
    return [f"{f}_{s}" for s in SCALES for f in feats]


def bonferroni_select(screen: pd.DataFrame) -> dict:
    """
    Apply the legacy Bonferroni rule to a feature screen.

    :param screen: Screen with ``context``, ``window``, ``feature`` and ``pval``.
    :return: Dictionary of context to selected feature columns (first passing scale).
    """
    out = {}
    for context in CONTEXTS:
        n = 8 if context in CPG_CONTEXTS else 13
        sel = screen[(screen.context == context) & (screen.pval <= 0.05 / 4 / n)]
        sel = sel.drop_duplicates("feature")
        if context in CPG_CONTEXTS:
            sel = sel[~sel.feature.isin(METHYLATION_FEATURES)]
        cols = list(sel.feature + "_" + sel.window)
        if cols:
            out[context] = cols
    return out


def run_screen(df: pd.DataFrame) -> pd.DataFrame:
    """
    Run Ruchit's univariate logistic feature screen (one fit per context and column).

    :param df: Training rows with ``context``, ``y`` and the feature columns.
    :return: Screen with ``context``, ``window``, ``feature`` and ``pval``.
    """
    rows = []
    for context in CONTEXTS:
        ctx = df[df.context == context]
        for scale in SCALES:
            for feature in SCREEN_ORDER:
                col = f"{feature}_{scale}"
                d = ctx[["y", col]].dropna()
                x = scipy.stats.zscore(d[col].to_numpy())
                try:
                    fit = sm.Logit(d.y.to_numpy(), sm.add_constant(x)).fit(disp=0)
                    pval = fit.pvalues[1]
                except Exception:  # pylint: disable=broad-except
                    pval = np.nan
                rows.append(
                    {
                        "context": context,
                        "window": scale,
                        "feature": feature,
                        "pval": pval,
                    }
                )
    return pd.DataFrame(rows)


def fit_logistic(
    d: pd.DataFrame, cols: list, c_value: float = None, weight: str = None
) -> dict:
    """
    Fit a z-scored logistic model of ``y`` and return what rr needs.

    :param d: Training rows of one context.
    :param cols: Feature columns.
    :param c_value: Inverse ridge penalty (sklearn ``C``). None fits unpenalized.
    :param weight: Optional column of row weights.
    :return: Dictionary with ``cols``, ``mean``, ``sd`` and ``coef`` (per z unit).
    """
    d = d.dropna(subset=cols)
    x = d[cols].to_numpy()
    mean, sd = x.mean(axis=0), x.std(axis=0)
    sd[sd == 0] = 1.0
    z = (x - mean) / sd
    model = LogisticRegression(
        penalty=None if c_value is None else "l2",
        C=1.0 if c_value is None else c_value,
        max_iter=2000,
    )
    model.fit(z, d.y.to_numpy(), sample_weight=None if weight is None else d[weight])
    return {
        "cols": cols,
        "mean": mean,
        "sd": sd,
        "coef": model.coef_[0],
        "intercept": model.intercept_[0],
    }


def logistic_rr(model: dict, windows: pd.DataFrame) -> np.ndarray:
    """
    Get rr = p(window) / p(average window) on the logistic (sigma) scale.

    :param model: Output of :func:`fit_logistic`.
    :param windows: Window features.
    :return: rr per window (NaN where a feature is missing).
    """
    z = (windows[model["cols"]].to_numpy() - model["mean"]) / model["sd"]
    eta0 = model["intercept"]
    p = 1 / (1 + np.exp(-(eta0 + z @ model["coef"])))
    return p / (1 / (1 + np.exp(-eta0)))


class Evaluator:
    """Score rr matrices on Ruchit's evaluation export, as in his report."""

    def __init__(
        self,
        cells_path: str,
        windows_path: str,
        transcripts_path: str,
        release_tsv: str,
    ):
        """
        Load the evaluation export.

        :param cells_path: syn cells TSV (with or without rr columns).
        :param windows_path: windows TSV.
        :param transcripts_path: transcripts TSV.
        :param release_tsv: v4.1.2 release TSV, for release QC gene flags.
        """
        self.windows = pd.read_csv(windows_path, sep="\t", compression="gzip")
        self.windows["fold"] = chrom_fold(self.windows.contig)
        win_idx = pd.Series(np.arange(len(self.windows)), index=self.windows.win)
        cells = pd.read_csv(cells_path, sep="\t", na_values=["NA"])
        self.cells = cells
        self.cell_win = win_idx.reindex(cells.win).to_numpy()
        self.cell_ctx = cells.context.map(
            {c: i for i, c in enumerate(CONTEXTS)}
        ).to_numpy()
        tx = pd.read_csv(transcripts_path, sep="\t", compression="gzip").set_index(
            "transcript"
        )
        self.tx_codes, self.tx_names = pd.factorize(cells.transcript)
        tx = tx.reindex(self.tx_names)
        self.exp_cell = cells.exp.to_numpy()
        self.obs = np.bincount(self.tx_codes, weights=cells.obs.to_numpy())
        self.exp0 = np.bincount(self.tx_codes, weights=self.exp_cell)
        unflagged = tx.unflagged.astype(str).str.lower().eq("true").to_numpy()
        self.unflagged = unflagged
        self.core = unflagged & (self.exp0 >= 20)
        self.slope_set = unflagged & (self.exp0 >= 50)
        self.tx_fold = chrom_fold(tx.chrom.astype(str))
        rel = pd.read_csv(
            release_tsv,
            sep="\t",
            usecols=["transcript", "gene_flags"],
            compression="gzip",
            na_values=["NA"],
        ).drop_duplicates("transcript")
        flags = rel.set_index("transcript").gene_flags.reindex(self.tx_names)
        self.no_qc_flag = (
            flags.fillna("[]").str.strip("[]").str.replace('"', "").eq("").to_numpy()
        )

    def cell_rr(self, rr: np.ndarray) -> np.ndarray:
        """
        Map a (window x context) rr matrix to cells and normalize per context.

        :param rr: Array of shape (n_windows, 32); NaN means no rr (1).
        :return: Normalized rr per cell.
        """
        ok = ~np.isnan(self.cell_win)
        r = np.ones(len(self.exp_cell))
        r[ok] = rr[self.cell_win[ok].astype(int), self.cell_ctx[ok]]
        r = np.where(np.isnan(r), 1.0, r)
        num = np.bincount(self.cell_ctx, weights=self.exp_cell * r, minlength=32)
        den = np.bincount(self.cell_ctx, weights=self.exp_cell, minlength=32)
        return r / (num / den)[self.cell_ctx]

    def expected(self, rr: np.ndarray, alpha: float = 1.0) -> np.ndarray:
        """
        Get per-transcript expected counts under rr^alpha (normalized after alpha).

        :param rr: (n_windows, 32) rr matrix, or None for no correction.
        :param alpha: Strength.
        :return: Expected count per transcript.
        """
        if rr is None:
            return self.exp0
        r = self.cell_rr(np.power(rr, alpha))
        return np.bincount(self.tx_codes, weights=self.exp_cell * r)

    def noise(self, exp: np.ndarray, mask: np.ndarray) -> float:
        """
        Get level-free mean z^2 over transcripts in ``mask`` (expected rescaled).

        :param exp: Expected per transcript.
        :param mask: Transcripts to use.
        :return: Mean z^2.
        """
        o, e = self.obs[mask], exp[mask]
        e = e * o.sum() / e.sum()
        return float(np.mean((o - e) ** 2 / e))

    def score(self, rr: np.ndarray, alpha: float = 1.0) -> dict:
        """
        Score an rr matrix with the report's metrics.

        :param rr: (n_windows, 32) rr matrix, or None for no correction.
        :param alpha: Strength.
        :return: Dictionary of metrics.
        """
        exp = self.expected(rr, alpha)
        z = (self.obs - exp) / np.sqrt(exp)
        out = {
            "core_noise": self.noise(exp, self.core),
            "core_noise_odd": self.noise(exp, self.core & (self.tx_fold == "odd")),
            "core_noise_even": self.noise(exp, self.core & (self.tx_fold == "even")),
            "outliers": int((np.abs(z) > 8).sum()),
            "outliers_no_qc_flag": int(((np.abs(z) > 8) & self.no_qc_flag).sum()),
        }
        if rr is not None:
            m = self.slope_set
            x = np.log(exp[m] / self.exp0[m])
            y = np.log(self.obs[m] / self.exp0[m])
            w = self.exp0[m]
            xm, ym = np.average(x, weights=w), np.average(y, weights=w)
            out["slope"] = float(
                np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2)
            )
        return out

    def best_alpha(self, rr: np.ndarray) -> dict:
        """
        Choose alpha on one chromosome half and score it on the other.

        :param rr: (n_windows, 32) rr matrix.
        :return: Chosen alpha per half and the cross-fit core noise.
        """
        curves = {a: self.score(rr, a) for a in ALPHAS}
        best_odd = min(ALPHAS, key=lambda a: curves[a]["core_noise_odd"])
        best_even = min(ALPHAS, key=lambda a: curves[a]["core_noise_even"])
        n_odd = (self.core & (self.tx_fold == "odd")).sum()
        n_even = (self.core & (self.tx_fold == "even")).sum()
        crossfit = (
            curves[best_odd]["core_noise_even"] * n_even
            + curves[best_even]["core_noise_odd"] * n_odd
        ) / (n_odd + n_even)
        return {
            "alpha_odd": best_odd,
            "alpha_even": best_even,
            "crossfit_noise": crossfit,
        }


def oof_rr(windows: pd.DataFrame, fit_fold) -> np.ndarray:
    """
    Build an out-of-fold rr matrix.

    :param windows: Window features with a ``fold`` column.
    :param fit_fold: Function (fold) -> {context: rr function of window rows}, with
        models trained without ``fold``.
    :return: (n_windows, 32) rr matrix.
    """
    rr = np.full((len(windows), len(CONTEXTS)), np.nan)
    for fold in FOLDS:
        rows = np.where(windows.fold == fold)[0]
        sub = windows.iloc[rows]
        for context, predict in fit_fold(fold).items():
            rr[rows, CONTEXTS.index(context)] = predict(sub)
    return rr


def write_recommended(ev: Evaluator, train: pd.DataFrame, penalty: float, args) -> None:
    """
    Fit the recommended rr (ridge, all features, all chromosomes) and score it.

    Writes the cells TSV with an added ``rr_recommended`` column (normalized per
    context; 1 where a window lacks a feature, which covers every chrX/Y window) and
    a metrics TSV: core noise, slope and outliers at each strength alpha, and chrX/Y
    noise.

    :param ev: Evaluator.
    :param train: Training rows of the chosen DNM set.
    :param penalty: Per-site ridge penalty (mean log loss + penalty / 2 * sum(coef^2)),
        as ``RR_RIDGE_PENALTY`` in ``regional_rate.py``.
    :param args: Parsed arguments with ``recommended_cells_out`` and
        ``recommended_metrics_out``.
    :return: None.
    """
    rr = np.full((len(ev.windows), len(CONTEXTS)), np.nan)
    for context in CONTEXTS:
        cols = allowed_columns(context)
        d = train[train.context == context]
        n = len(d.dropna(subset=cols))
        m = fit_logistic(d, cols, 1.0 / (n * penalty))
        rr[:, CONTEXTS.index(context)] = logistic_rr(m, ev.windows)
    ev.cells.assign(rr_recommended=ev.cell_rr(rr)).to_csv(
        args.recommended_cells_out, sep="\t", index=False
    )
    tx = pd.read_csv(args.transcripts, sep="\t", compression="gzip").set_index(
        "transcript"
    )
    xy = tx.reindex(ev.tx_names).is_xy.astype(str).str.lower().eq("true").to_numpy()
    rows = []
    for alpha in [0.0, *ALPHAS]:
        sc = ev.score(rr, alpha) if alpha > 0 else ev.score(None)
        exp = ev.expected(rr, alpha) if alpha > 0 else ev.exp0
        rows.append({"alpha": alpha, **sc, "xy_noise": ev.noise(exp, xy)})
    pd.DataFrame(rows).to_csv(args.recommended_metrics_out, sep="\t", index=False)
    logger.info("Recommended rr: %s", rows)


def main(args: argparse.Namespace) -> None:
    """Run the experiments and write one row per variant."""
    ev = Evaluator(args.cells, args.windows, args.transcripts, args.release_tsv)
    if "recommended" in args.experiments:
        path = (
            args.train_satt if args.recommended_set == "Satterstrom" else args.train_17k
        )
        train = pd.read_csv(path, sep="\t", compression="gzip", na_values=["NA"])
        write_recommended(ev, train, args.recommended_penalty, args)
        if args.experiments == ["recommended"]:
            return
    results = []

    def record(name: str, train_set: str, rr: np.ndarray, note: str = "") -> None:
        row = {"variant": name, "dnm_set": train_set, "note": note, **ev.score(rr)}
        if rr is not None:
            row.update(ev.best_alpha(rr))
        logger.info("%s", row)
        results.append(row)

    record("no correction", "", None)
    # Check the evaluator against Ruchit's own calibrated tables (all-chromosome).
    for name, col in (("ruchit table (17k)", "rr_ruchit_17k"),):
        if col in ev.cells:
            r = ev.cells[col].fillna(1.0).to_numpy()
            rr_cells = np.full((len(ev.windows), 32), np.nan)
            ok = ~np.isnan(ev.cell_win)
            rr_cells[ev.cell_win[ok].astype(int), ev.cell_ctx[ok]] = r[ok]
            record(name, "17k", rr_cells, "his calibrated table, all-chromosome models")

    wfeat = ev.windows
    for train_set, path, screen_path in (
        ("17k", args.train_17k, args.screen_17k),
        ("Satterstrom", args.train_satt, args.screen_satt),
    ):
        train = pd.read_csv(path, sep="\t", compression="gzip", na_values=["NA"])
        train["fold"] = chrom_fold(train.locus.str.split(":").str[0])
        screen = pd.read_csv(screen_path, sep="\t")
        selected = bonferroni_select(screen)

        # 1. Ruchit's method: Bonferroni features (screened on all data), unpenalized.
        def fit_baseline(fold):
            out = {}
            for context, cols in selected.items():
                d = train[(train.context == context) & (train.fold != fold)]
                m = fit_logistic(d, cols)
                out[context] = lambda w, m=m: logistic_rr(m, w)
            return out

        if "baseline" in args.experiments:
            record(
                "Bonferroni, unpenalized (his method)",
                train_set,
                oof_rr(wfeat, fit_baseline),
            )

        # 2. Ridge on every allowed feature, one penalty for all contexts. Also record
        #    the held-out DNM log loss, so the penalty can be chosen without gnomAD.
        for c_value in args.ridge_c if "ridge" in args.experiments else []:
            heldout = []

            def fit_ridge(fold, c_value=c_value, heldout=heldout):
                out = {}
                for context in CONTEXTS:
                    cols = allowed_columns(context)
                    d = train[(train.context == context) & (train.fold != fold)]
                    m = fit_logistic(d, cols, c_value)
                    out[context] = lambda w, m=m: logistic_rr(m, w)
                    t = train[(train.context == context) & (train.fold == fold)]
                    t = t.dropna(subset=cols)
                    z = (t[cols].to_numpy() - m["mean"]) / m["sd"]
                    p = 1 / (1 + np.exp(-(m["intercept"] + z @ m["coef"])))
                    p = np.clip(p, 1e-12, 1 - 1e-12)
                    y = t.y.to_numpy()
                    heldout.append(-np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))
                return out

            rr = oof_rr(wfeat, fit_ridge)
            record(
                f"ridge, all features, C={c_value:g}",
                train_set,
                rr,
                f"held-out DNM log loss {sum(heldout):.1f}",
            )

        # 3. Sample splitting: screen on one half of the training sites, fit on
        # the other.
        half = train.locus.map(lambda s: zlib.crc32(s.encode()) % 2)

        def fit_split(fold):
            tr = train[train.fold != fold]
            sel = bonferroni_select(run_screen(tr[half[tr.index] == 0]))
            out = {}
            for context, cols in sel.items():
                d = tr[(tr.context == context) & (half[tr.index] == 1)]
                m = fit_logistic(d, cols)
                out[context] = lambda w, m=m: logistic_rr(m, w)
            return out

        if "split" in args.experiments:
            record(
                "sample split: screen half A, fit half B",
                train_set,
                oof_rr(wfeat, fit_split),
            )

        # 4. Background design: weight background windows from the exported cells, either
        #    sites absent from gnomAD (current design) or all possible sites.
        agg = (
            ev.cells.groupby(["win", "context"], as_index=False)[
                ["n_sites", "obs"]
            ].sum()
            if "n_sites" in ev.cells and "background" in args.experiments
            else None
        )
        if agg is not None:
            bg = agg.merge(wfeat[["win", "fold", *COLUMNS]], on="win")
            dnm = train[train.y == 1][["context", "fold", *COLUMNS]].assign(y=1, w=1.0)
            for label, weight in (
                ("absent from gnomAD (current design)", 0.05 * (bg.n_sites - bg.obs)),
                ("all possible sites", 0.05 * bg.n_sites),
            ):
                bgw = bg[["context", "fold", *COLUMNS]].assign(y=0, w=weight.to_numpy())
                both = pd.concat([dnm, bgw[bgw.w > 0]], ignore_index=True)

                def fit_bg(fold, both=both):
                    out = {}
                    for context, cols in selected.items():
                        d = both[(both.context == context) & (both.fold != fold)]
                        m = fit_logistic(d, cols, weight="w")
                        out[context] = lambda w, m=m: logistic_rr(m, w)
                    return out

                record(
                    f"background: {label}",
                    train_set,
                    oof_rr(wfeat, fit_bg),
                    "background from exported MANE Select cells, Bonferroni features",
                )

    # 5. No DNMs: fit regional effects on gnomAD synonymous counts themselves
    #    (Poisson, offset log expected, ridge), training on the other chromosome half.
    agg = ev.cells.groupby(["win", "context"], as_index=False)[["exp", "obs"]].sum()
    agg = agg.merge(wfeat[["win", "fold", *COLUMNS]], on="win")
    for alpha_pen in args.poisson_alpha if "gnomad" in args.experiments else []:

        def fit_gnomad(fold, alpha_pen=alpha_pen):
            out = {}
            for context in CONTEXTS:
                cols = allowed_columns(context)
                d = agg[(agg.context == context) & (agg.fold != fold)].dropna(
                    subset=cols
                )
                x = d[cols].to_numpy()
                mean, sd = x.mean(axis=0), x.std(axis=0)
                sd[sd == 0] = 1.0
                model = PoissonRegressor(alpha=alpha_pen, max_iter=1000)
                model.fit((x - mean) / sd, d.obs / d.exp, sample_weight=d.exp)
                out[context] = lambda w, mo=model, mean=mean, sd=sd, cols=cols: np.exp(
                    ((w[cols].to_numpy() - mean) / sd) @ mo.coef_
                )
            return out

        record(
            f"no DNMs: Poisson on gnomAD syn, alpha={alpha_pen:g}",
            "",
            oof_rr(wfeat, fit_gnomad),
            "trained on the other chromosome half's exome synonymous counts",
        )

    pd.DataFrame(results).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--windows", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument("--train-17k", required=True)
    parser.add_argument("--screen-17k", required=True)
    parser.add_argument("--train-satt", required=True)
    parser.add_argument("--screen-satt", required=True)
    parser.add_argument(
        "--ridge-c", type=float, nargs="+", default=[1e-4, 1e-3, 1e-2, 1e-1]
    )
    parser.add_argument(
        "--poisson-alpha", type=float, nargs="+", default=[1e-4, 1e-3, 1e-2]
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=["baseline", "ridge", "split", "background", "gnomad"],
    )
    parser.add_argument(
        "--recommended-set", default="Satterstrom", choices=["17k", "Satterstrom"]
    )
    parser.add_argument("--recommended-penalty", type=float, default=0.3)
    parser.add_argument("--recommended-cells-out")
    parser.add_argument("--recommended-metrics-out")
    parser.add_argument("--out")
    main(parser.parse_args())
