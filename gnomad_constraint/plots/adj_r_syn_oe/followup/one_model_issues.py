"""Score every tried fix against each problem of the current expected counts.

For the notes ``notes_one_model_for_expected.md``, linked from the rr report.

"Current" is the v4.1.1 expected count without any rr: possible * stored plateau line *
coverage correction. The fixes are the full run's ridge rr, Roulette's per-variant
rates (with and without the ridge), four plateau versions (weighted lines, the
saturating and logistic curves, one line per substitution type; each as
``f_new(mu) / f_line(mu)`` on the current count, as in the pipeline test), and the
single cloglog models. Everything is scored on the exported MANE Select synonymous
cells (the one models out of fold, from ``one_model_fit.py``), except the pLoF and NDD
side effects, which come from ``one_model_release_eval.py`` (one models, rr, Roulette)
and the plateau pipeline runs (``score_release.py``, ``disease_gene_pr.py``).

Writes ``data/one_model/issues.{summary,contexts,coverage,deciles,hotspots}.tsv``, and
``data/one_model/issues.table.md``: each problem against each fix, judged by how much of
the current gap to the target the fix closes (yes: at least 80%; partly: 20-80%; no:
-10% to 20%; worse: the gap grows by more than 10%).
"""

import json

import numpy as np
import pandas as pd
from one_model_fit import prep_cells
from outlier_drivers import HISTONE, MATERNAL_HOTSPOTS
from remap_rr_to_shet import read_release_flags

CURRENT = "current stack, no rr"
METHODS = [
    (CURRENT, "Current (v4.1.1, no rr)"),
    ("current stack + ridge rr", "Ridge rr"),
    ("current stack x Roulette", "× Roulette"),
    ("current stack x Roulette + ridge rr", "× Roulette + ridge rr"),
    ("plateau: weighted", "Plateau: weighted lines"),
    ("plateau: saturating", "Plateau: saturating curve"),
    ("plateau: logistic", "Plateau: logistic curve"),
    ("plateau: by_substitution", "Plateau: 1 line per substitution"),
    ("one model: groups + coverage", "One model: groups + coverage"),
    ("one model: + window MQ", "One model: + MQ"),
    ("one model: + DNM regional + window MQ", "One model: + DNM regional + MQ"),
    ("one model: + gnomAD regional + window MQ", "One model: + gnomAD regional + MQ"),
    ("one model: + DNM regional + window MQ + Roulette", "One model: full"),
]
"""Methods in the table: (name in the inputs, label)."""
# The v4 refit only added the per-group outputs; its per-transcript output is
# identical to v3's.
PREFIX = "data/one_model/one_model_full_v4"
PER_TRANSCRIPT = "data/one_model/one_model_full_v3.by_csq.per_transcript.tsv.gz"
OUT = "data/one_model/issues"
PLATEAU_RUNS = {
    "plateau: saturating": "plateau_saturating",
    "plateau: logistic": "plateau_logistic",
    "plateau: by_substitution": "plateau_by_substitution",
}
ONE_MODEL_KEYS = {
    "current stack + ridge rr": "stack_ridge",
    "current stack x Roulette": "stack_roulette",
    "current stack x Roulette + ridge rr": "stack_roulette_ridge",
    "one model: groups + coverage": "om_groups_coverage",
    "one model: + window MQ": "om_mq",
    "one model: + DNM regional + window MQ": "om_dnm_mq",
    "one model: + gnomAD regional + window MQ": "om_gnomad_mq",
    "one model: + DNM regional + window MQ + Roulette": "om_full",
}


def plateau_ratios(c: pd.DataFrame) -> dict:
    """Get f_new(mu) / f_line(mu) per cell for each plateau version."""
    params = json.load(open("plateau_params.json"))
    fits = pd.read_csv("plateau_models.fits.tsv", sep="\t")
    x = c.mu.to_numpy() * 1e8
    group = np.where(c.cpg, "CpG transition", "everything else")

    def lin(d, key):
        a = np.array([d[k][0] for k in key])
        b = np.array([d[k][1] for k in key])
        return a, b

    a, b = lin(params["stored"], group)
    f_line = a + b * x
    w = fits[fits.fit == "weighted"].set_index("line")
    weighted = {k: [w.loc[k, "intercept"], w.loc[k, "slope"] * 1e-8] for k in w.index}
    out = {}
    a, b = lin(weighted, group)
    out["plateau: weighted"] = (a + b * x) / f_line
    a, b = lin(params["saturating"], group)
    out["plateau: saturating"] = (1 - np.exp(-(a + b * x))) / f_line
    a, b = lin(params["logistic"], group)
    out["plateau: logistic"] = 1 / (1 + np.exp(-(a + b * x))) / f_line
    a, b = lin(params["by_substitution"], c["sub"])
    out["plateau: by_substitution"] = (a + b * x) / f_line
    return out


def main() -> None:
    """Score every method."""
    names = [m for m, _ in METHODS]
    # Synonymous cells with the plateau versions' expected counts.
    c = prep_cells(
        pd.read_csv(
            "data/one_model/syn_model_cells.tsv.bgz",
            sep="\t",
            compression="gzip",
            usecols=[
                "transcript",
                "context",
                "ref",
                "alt",
                "methylation_level",
                "cpg",
                "coverage",
                "possible",
                "observed",
                "mu",
                "e_old",
                "win",
            ],
        ),
        test=False,
    )
    # The cells the one-model fit keeps: those in a window with features or MQ.
    wins = set(
        pd.read_csv(
            "data/rr_eval_export/windows.tsv.bgz",
            sep="\t",
            compression="gzip",
            usecols=["win"],
        ).win
    ) | set(
        pd.read_csv(
            "data/one_model/syn_model_cells.window_mq.tsv.bgz",
            sep="\t",
            compression="gzip",
            usecols=["win"],
        ).win
    )
    c = c[c.win.isin(wins)].reset_index(drop=True)
    for k, r in plateau_ratios(c).items():
        c[k] = c.e_old.to_numpy() * r
    plat = list(PLATEAU_RUNS) + ["plateau: weighted"]
    plat_tx = c.groupby("transcript")[plat].sum()
    plat_groups = c.groupby(["ctx", "coverage"])[["observed", "e_old", *plat]].sum()

    # The one models, rr and Roulette, out of fold.
    pt = pd.read_csv(PER_TRANSCRIPT, sep="\t")
    pt = pt[pt.csq == "syn"].set_index("transcript")
    t = pt.join(plat_tx, how="inner")
    bg = pd.read_csv(f"{PREFIX}.by_group.tsv.gz", sep="\t")
    bg = prep_cells(bg.assign(win="chr1-0-1000", mu=1.0), test=False)
    groups = bg.groupby(["ctx", "coverage"])[
        ["observed"] + [n for n in names if n not in plat]
    ].sum()
    groups = groups.join(plat_groups[plat])

    # Gene sets.
    tx = pd.read_csv(
        "data/rr_eval_export/transcripts.tsv.bgz", sep="\t", compression="gzip"
    ).set_index("transcript")
    rel = read_release_flags("../gnomad.v4.1.2.constraint_metrics.tsv.bgz")
    rel = rel.drop_duplicates("transcript").set_index("transcript")
    t = t.join(tx[["unflagged", "decile"]]).join(rel[["gene", "gene_flags"]])
    t["unflagged"] = t.unflagged.astype(str).str.lower().eq("true")
    flags = t.gene_flags.map(lambda f: f if isinstance(f, list) else [])
    no_qc = flags.map(len).eq(0)
    low_mq = flags.map(lambda f: "low_exome_mapping_quality" in f)
    histone = t.gene.fillna("").str.match(HISTONE)
    genes = pd.read_csv("outlier_drivers.all_genes.tsv", sep="\t").set_index(
        "transcript"
    )
    g = genes.reindex(t.index)
    outlier = g.outlier.astype(str).str.lower().eq("true")
    e0 = t[CURRENT]
    core = t.unflagged & (e0 >= 20)
    free = no_qc & (e0 >= 20)
    control = no_qc & ~histone & ~outlier & (e0 >= 20)

    def noise(o, e):
        e = e * o.sum() / e.sum()
        return float(np.mean((o - e) ** 2 / e))

    # Side effects from the pLoF runs.
    om_sc = pd.read_csv(
        "data/one_model/one_model_release.scores.tsv", sep="\t"
    ).set_index("model")
    om_pr = pd.read_csv("data/one_model/one_model_pr.auprc.tsv", sep="\t")
    pl_sc = pd.read_csv("release_plateau.scores.tsv", sep="\t").set_index("version")
    pl_pr = pd.read_csv("disease_gene_pr_plateau.auprc.tsv", sep="\t")

    def ndd(pr, key):
        r = pr[(pr.test == "NDD") & (pr.gene_set == "all genes") & (pr.version == key)]
        return r.iloc[0] if len(r) else None

    rows, ctx_rows, cov_rows, dec_rows, hot_rows = [], [], [], [], []
    tot = groups.sum()
    for m, label in METHODS:
        o, e = t.observed, t[m]
        lvl = o[core].sum() / e[core].sum()
        z = (o - e * lvl) / np.sqrt(e * lvl)
        big = z.abs() > 8
        cols = list(dict.fromkeys(["observed", m, CURRENT]))
        by_ctx = groups.groupby(level="ctx")[cols].sum()
        oe_ctx = by_ctx.observed / by_ctx[m] / (tot.observed / tot[m])
        wts = by_ctx[CURRENT] / by_ctx[CURRENT].sum()
        by_cov = groups.groupby(level="coverage")[["observed", m]].sum()
        oe_cov = by_cov.observed / by_cov[m]

        def cov_bin(lo, hi):
            s = by_cov.loc[lo : hi - 1]
            return float(s.observed.sum() / s[m].sum())

        d = t.dropna(subset=["decile"]).groupby("decile")[["observed", m]].sum()
        oe_dec = d.observed / d[m]
        row = {
            "method": m,
            "label": label,
            "syn_oe": float(tot.observed / tot[m]),
            "context_sd": float(np.sqrt(np.sum(wts * (oe_ctx - 1) ** 2))),
            "cov_20_49": cov_bin(20, 50),
            "cov_80_89": cov_bin(80, 90),
            "core_noise": noise(o[core], e[core]),
            "outcome_free_noise": noise(o[free], e[free]),
            "loeuf_top_rise": float(oe_dec.iloc[-1] - oe_dec.iloc[:6].mean()),
            "outliers_no_flag_excess": int((big & (z > 0) & no_qc & ~histone).sum()),
            "outliers_no_flag_deficit": int((big & (z < 0) & no_qc & ~histone).sum()),
            "low_mq_deficit": int((big & (z < 0) & low_mq).sum()),
            "low_mq_excess": int((big & (z > 0) & low_mq).sum()),
            "histone_outliers": int((big & histone).sum()),
            "histone_oe": float(o[histone].sum() / e[histone].sum() / lvl),
        }
        for name, (chrom, start, end) in MATERNAL_HOTSPOTS.items():
            r = g.chromosome.eq(chrom) & g.mid.between(start, end) & control
            v = float(o[r].sum() / e[r].sum() / lvl)
            hot_rows.append({"method": m, "region": name, "oe": v})
        row["hotspot_max"] = max(h["oe"] for h in hot_rows if h["method"] == m)
        # Side effects: pLoF and missense levels relative to synonymous, NDD AUPRC.
        if m == CURRENT:
            row.update(lof_vs_syn=0.0, mis_vs_syn=0.0, ndd_change=0.0)
        elif m in om_sc.index:
            ref = om_sc.loc[CURRENT]
            row["lof_vs_syn"] = 100 * (
                om_sc.loc[m, "lof_over_syn"] / ref.lof_over_syn - 1
            )
            row["mis_vs_syn"] = 100 * (
                om_sc.loc[m, "mis_over_syn"] / ref.mis_over_syn - 1
            )
            r = ndd(om_pr, ONE_MODEL_KEYS[m])
            row.update(ndd_change=r.diff_vs_none, ndd_lo=r.diff_lo, ndd_hi=r.diff_hi)
        elif m in PLATEAU_RUNS:
            k, ref = PLATEAU_RUNS[m], pl_sc.loc["none"]
            row["lof_vs_syn"] = 100 * (
                pl_sc.loc[k, "lof_oe"]
                / pl_sc.loc[k, "syn_oe"]
                / (ref.lof_oe / ref.syn_oe)
                - 1
            )
            row["mis_vs_syn"] = 100 * (
                pl_sc.loc[k, "mis_oe"]
                / pl_sc.loc[k, "syn_oe"]
                / (ref.mis_oe / ref.syn_oe)
                - 1
            )
            r = ndd(pl_pr, k)
            row.update(ndd_change=r.diff_vs_none, ndd_lo=r.diff_lo, ndd_hi=r.diff_hi)
        rows.append(row)
        ctx_rows += [{"method": m, "ctx": k, "oe": v} for k, v in oe_ctx.items()]
        cov_rows += [{"method": m, "coverage": k, "oe": v} for k, v in oe_cov.items()]
        dec_rows += [
            {"method": m, "decile": int(k), "oe": v} for k, v in oe_dec.items()
        ]
    for name, data in (
        ("summary", rows),
        ("contexts", ctx_rows),
        ("coverage", cov_rows),
        ("deciles", dec_rows),
        ("hotspots", hot_rows),
    ):
        pd.DataFrame(data).to_csv(f"{OUT}.{name}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).drop(columns="method").round(3).to_string(index=False))


ISSUES = [
    (
        "Synonymous level is off 1",
        "total obs/exp",
        lambda r: abs(r.syn_oe - 1),
        "{:.3f}",
        "syn_oe",
    ),
    (
        "Contexts are offset from each other",
        "weighted SD of the 32 contexts' obs/exp",
        lambda r: r.context_sd,
        "{:.1%}",
        "context_sd",
    ),
    (
        "Low coverage is miscalibrated",
        "obs/exp at AN 20–49% / 80–89%",
        lambda r: max(abs(r.cov_20_49 - 1), abs(r.cov_80_89 - 1)),
        None,
        None,
    ),
    (
        "Regional variation: per-gene noise",
        "noise, core / no QC flag (1 = Poisson)",
        lambda r: (r.core_noise - 1 + r.outcome_free_noise - 1) / 2,
        None,
        None,
    ),
    (
        "Maternal hotspot regions run high",
        "control genes' obs/exp, highest region (8p)",
        lambda r: abs(r.hotspot_max - 1),
        "{:.2f}",
        "hotspot_max",
    ),
    (
        "Least constrained LOEUF decile runs high",
        "obs/exp in decile 10 minus deciles 1–6",
        lambda r: r.loeuf_top_rise,
        "{:+.3f}",
        "loeuf_top_rise",
    ),
    (
        "Outliers with no QC flag",
        "genes with \\|z\\| > 8",
        lambda r: r.outliers_no_flag_excess + r.outliers_no_flag_deficit,
        "{:.0f}",
        None,
    ),
    (
        "Low-mapping-quality genes look depleted",
        "flagged genes with z < −8",
        lambda r: r.low_mq_deficit,
        "{:.0f}",
        "low_mq_deficit",
    ),
    (
        "Histone genes run high",
        "outliers",
        lambda r: r.histone_outliers,
        "{:.0f}",
        "histone_outliers",
    ),
]
"""Problems: (name, metric, gap to the target, format of the current value, column)."""


def verdict(f: float) -> str:
    """Get the verdict for the share of the current gap a fix closes."""
    if f >= 0.8:
        return "yes"
    if f >= 0.2:
        return "partly"
    return "no" if f > -0.1 else "worse"


def table() -> None:
    """Write the problems-by-fixes table as markdown."""
    sm = pd.read_csv(f"{OUT}.summary.tsv", sep="\t").set_index("method")
    cur = sm.loc[CURRENT]
    fixes = [m for m, _ in METHODS if m != CURRENT]
    short = {
        "current stack + ridge rr": "Ridge rr",
        "current stack x Roulette": "× Roulette",
        "current stack x Roulette + ridge rr": "× Roulette + ridge",
        "plateau: weighted": "Plateau: weighted",
        "plateau: saturating": "Plateau: saturating",
        "plateau: logistic": "Plateau: logistic",
        "plateau: by_substitution": "Plateau: per substitution",
        "one model: groups + coverage": "1M: groups + coverage",
        "one model: + window MQ": "1M: + MQ",
        "one model: + DNM regional + window MQ": "1M: + DNM regional + MQ",
        "one model: + gnomAD regional + window MQ": "1M: + gnomAD regional + MQ",
        "one model: + DNM regional + window MQ + Roulette": "1M: full",
    }
    lines = [
        "| Problem (metric) | Current | " + " | ".join(short[m] for m in fixes) + " |",
        "|---|---|" + "---|" * len(fixes),
    ]
    rows = []
    for name, metric, gap, fmt, col in ISSUES:
        if name.startswith("Low coverage"):
            now = f"{cur.cov_20_49:.2f} / {cur.cov_80_89:.2f}"
        elif name.startswith("Regional"):
            now = f"{cur.core_noise:.2f} / {cur.outcome_free_noise:.2f}"
        elif name.startswith("Outliers"):
            now = f"{cur.outliers_no_flag_excess + cur.outliers_no_flag_deficit:.0f}"
        else:
            now = fmt.format(cur[col])
        g0 = gap(cur)
        cells = []
        for m in fixes:
            v = verdict(1 - gap(sm.loc[m]) / g0)
            if (
                name.startswith("Low-mapping")
                and v == "yes"
                and sm.loc[m, "low_mq_excess"] >= 5
            ):
                v = "yes, overshoots"
            cells.append(v)
            rows.append(
                {
                    "problem": name,
                    "method": m,
                    "share_closed": 1 - gap(sm.loc[m]) / g0,
                    "verdict": v,
                }
            )
        lines.append(f"| {name} ({metric}) | {now} | " + " | ".join(cells) + " |")
    lof = [
        (
            "n/t"
            if pd.isna(sm.loc[m, "lof_vs_syn"])
            else f"{sm.loc[m, 'lof_vs_syn']:+.1f}%"
        )
        for m in fixes
    ]
    ndd = [
        "n/t" if pd.isna(sm.loc[m, "ndd_change"]) else f"{sm.loc[m, 'ndd_change']:+.3f}"
        for m in fixes
    ]
    lines.append(
        "| *Side effect:* pLoF expected vs synonymous | 0 | " + " | ".join(lof) + " |"
    )
    lines.append("| *Outcome:* NDD AUPRC change | 0 | " + " | ".join(ndd) + " |")
    open(f"{OUT}.table.md", "w").write("\n".join(lines) + "\n")
    pd.DataFrame(rows).to_csv(f"{OUT}.verdicts.tsv", sep="\t", index=False)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
    table()
