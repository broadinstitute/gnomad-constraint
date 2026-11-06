"""Score the single cloglog models with the report's other synonymous checks.

For the notes ``notes_one_model_for_expected.md``, linked from the rr report.

From ``one_model_fit.py``'s per-transcript expected counts (every model scored on the
chromosome half it wasn't fit on) and the exported cells' observed counts, on the same
transcripts as the fit:

- outliers (|z| > 8) by direction and category (low exome mapping quality, low exome
  coverage, histone genes, unflagged), as ``full_run_scores.outlier_categories``;
- the calibration slope against the current stack without rr: log obs/e0 on log e/e0,
  weighted by e0, on unflagged genes with e0 >= 50 (1 = the model's departures from the
  current stack are the right size), as ``full_run_scores.slope``;
- synonymous obs/exp by LOEUF decile (v4.1.2) and by s_het decile (Kar et al. 2026),
  decile 0 = most constrained for both;
- the maternal hotspot regions: pooled obs/exp, relative to the core set's, of the
  control genes and the excess outliers in each region, and how many of those are
  still outliers, as ``score_downstream.hotspots``.

Writes ``data/one_model/one_model_eval.{outliers,scores,deciles,hotspots}.tsv``.
"""

import numpy as np
import pandas as pd
from outlier_drivers import HISTONE, MATERNAL_HOTSPOTS
from remap_rr_to_shet import read_release_flags, read_shet

PREFIX = "data/one_model/one_model_full_v2"
OUT = "data/one_model/one_model_eval"


def main() -> None:
    """Score every model."""
    c = pd.read_csv(
        "data/one_model/syn_model_cells.tsv.bgz",
        sep="\t",
        compression="gzip",
        usecols=["transcript", "coverage", "possible", "observed"],
    )
    c = c[(c.coverage > 0) & (c.possible > 0)]
    obs = c.groupby("transcript")["observed"].sum()
    e = pd.read_csv(f"{PREFIX}.per_transcript.tsv.gz", sep="\t").set_index("transcript")
    models = list(e.columns)
    t = e.join(obs.rename("obs"), how="inner")

    tx = pd.read_csv(
        "data/rr_eval_export/transcripts.tsv.bgz", sep="\t", compression="gzip"
    ).set_index("transcript")
    rel = read_release_flags("../gnomad.v4.1.2.constraint_metrics.tsv.bgz").set_index(
        "transcript"
    )
    shet = read_shet("data/kar2026/supp_table3.xlsx").drop_duplicates("gene_id")
    t = t.join(tx[["unflagged", "decile"]]).join(rel[["gene", "gene_id", "gene_flags"]])
    t["unflagged"] = t.unflagged.astype(str).str.lower().eq("true")
    flags = t.gene_flags.map(lambda f: f if isinstance(f, list) else [])
    t["no_qc_flag"] = flags.map(len).eq(0)
    t = t.reset_index().merge(shet, on="gene_id", how="left").set_index("transcript")
    # Decile 0 = highest s_het, the most constrained, as for LOEUF.
    t["shet_decile"] = pd.qcut((-t.shet).rank(method="first"), 10, labels=False)
    histone = t.gene.fillna("").str.match(HISTONE)
    t["category"] = np.where(
        flags.map(lambda f: "low_exome_mapping_quality" in f),
        "low_exome_mapping_quality",
        np.where(
            flags.map(lambda f: "low_exome_coverage" in f),
            "low_exome_coverage",
            np.where(histone, "histone genes", "unflagged"),
        ),
    )

    e0 = t["current stack, no rr"]
    core = t.unflagged & (e0 >= 20)
    slope_set = t.unflagged & (e0 >= 50) & (t.obs > 0)
    genes = pd.read_csv("outlier_drivers.all_genes.tsv", sep="\t").set_index(
        "transcript"
    )
    g = t.join(genes[["chromosome", "mid", "outlier", "group"]], how="left")
    g["outlier"] = g.outlier.astype(str).str.lower().eq("true")
    control = g.no_qc_flag & ~histone & ~g.outlier & (e0 >= 20)
    hot = g.outlier & g.group.eq("maternal hotspot region")

    outliers, scores, deciles, hotspots = [], [], [], []
    for m in models:
        ex = t[m]
        lvl = t.obs[core].sum() / ex[core].sum()
        z = (t.obs - ex * lvl) / np.sqrt(ex * lvl)
        big = z.abs() > 8
        for (direction, cat), n in (
            t[big]
            .assign(direction=np.where(z[big] > 0, "excess", "deficit"))
            .groupby(["direction", "category"])
            .size()
            .items()
        ):
            outliers.append(
                {"model": m, "direction": direction, "category": cat, "n": n}
            )
        s = slope_set
        x, y, w = np.log(ex[s] / e0[s]), np.log(t.obs[s] / e0[s]), e0[s]
        xm, ym = np.average(x, weights=w), np.average(y, weights=w)
        scores.append(
            {
                "model": m,
                "slope_vs_current": (
                    float(np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2))
                    if x.std() > 0
                    else np.nan
                ),
                "sd_log_change": float(np.sqrt(np.average((x - xm) ** 2, weights=w))),
            }
        )
        for binning in ("decile", "shet_decile"):
            d = t.dropna(subset=[binning]).groupby(binning)
            oe = d.obs.sum() / d[m].sum()
            for dec, val in oe.items():
                deciles.append(
                    {"model": m, "binning": binning, "decile": int(dec), "oe": val}
                )
        for name, (chrom, start, end) in MATERNAL_HOTSPOTS.items():
            r = g.chromosome.eq(chrom) & g.mid.between(start, end)
            for label, k in (
                ("control genes", r & control),
                ("excess outliers", r & hot),
            ):
                hotspots.append(
                    {
                        "model": m,
                        "region": name,
                        "set": label,
                        "n": int(k.sum()),
                        "oe": t.obs[k].sum() / ex[k].sum() / lvl,
                        "still_outliers": int(z[k].gt(8).sum()),
                    }
                )
    for name, rows in (
        ("outliers", outliers),
        ("scores", scores),
        ("deciles", deciles),
        ("hotspots", hotspots),
    ):
        pd.DataFrame(rows).to_csv(f"{OUT}.{name}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(scores).round(3).to_string(index=False))
    o = pd.DataFrame(outliers).pivot_table(
        index="model", columns=["direction", "category"], values="n", fill_value=0
    )
    print(o.to_string())
    h = pd.DataFrame(hotspots)
    print(
        h[h.set == "excess outliers"]
        .pivot(index="model", columns="region", values="oe")
        .round(2)
        .to_string()
    )


if __name__ == "__main__":
    main()
