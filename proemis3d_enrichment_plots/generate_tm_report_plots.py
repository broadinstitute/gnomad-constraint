"""Generate all plots for the transmembrane constraint analysis report.

Reads:
    - full_forward_tm_genes.tsv (per-residue constraint from Proemis3D, full dataset)
    - full_genes_uniprot_features.json (UniProt TM/intramembrane annotations)
    - autism_dd_variants.tsv (DD de novo case variants)

Writes (to the same directory as this script):
    - fig1_tm_vs_nontm_distributions.png  (Figure 1)
    - fig2_ion_channel_segment_types.png   (Figure 2)
    - fig3_segment_type_histograms.png     (Figure 3)
    - fig4_category_summary_strip.png      (Figure 4)
    - fig5_plddt_bias_check.png            (Figure 5)
"""

import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, wilcoxon

matplotlib.use("Agg")

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent
FORWARD_TSV = REPO_ROOT / "full_forward_tm_genes.tsv"
FEATURES_JSON = REPO_ROOT / "full_genes_uniprot_features.json"
DD_VARIANTS_TSV = REPO_ROOT / "autism_dd_variants.tsv"

METHOD_GROUPS = [
    ("No filter", "AIC (no pLDDT/PAE)", "LRT (no pLDDT/PAE)"),
    ("pLDDT exclude", "AIC, pLDDT exclude", "LRT, pLDDT exclude"),
    ("PAE filter center", "AIC, PAE filter_center", "LRT, PAE filter_center"),
]

CAT_COLORS = {
    "Non-TM": "#999999",
    "Generic TM helix": "#2166ac",
    "Voltage sensor (S4)": "#d6604d",
    "VSD other (S1-S3)": "#f4a582",
    "Pore-lining (S5)": "#4393c3",
    "Pore-lining (S6)": "#053061",
    "Transporter helix": "#7fbc41",
    "Named TM (M1/M2)": "#b2abd2",
    "Intramembrane": "#e08214",
}
CAT_MARKERS = {
    "Non-TM": "D",
    "Generic TM helix": "o",
    "Voltage sensor (S4)": "^",
    "VSD other (S1-S3)": "v",
    "Pore-lining (S5)": "s",
    "Pore-lining (S6)": "P",
    "Transporter helix": "h",
    "Named TM (M1/M2)": "*",
    "Intramembrane": "X",
}

C_AIC = "#4393c3"
C_LRT = "#b2182b"
C_TM = "#2166ac"
C_NONTM = "#ef8a62"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def filter_has_nonnull(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only rows from (gene, run_label) pairs that have at least one
    non-null residue.

    This includes catch-all (null) residues for genes where the forward
    algorithm found *some* constraint structure, but drops gene-method
    combinations that are entirely null (no constraint regions at all).
    """
    nonnull_pairs = (
        df[df["is_null"] == False]
        .groupby(["uniprot_id", "run_label"])
        .size()
        .reset_index()[["uniprot_id", "run_label"]]
    )
    return df.merge(nonnull_pairs, on=["uniprot_id", "run_label"])


def geomean(arr):
    arr = np.asarray(arr, dtype=float)
    arr = arr[(arr > 0) & np.isfinite(arr)]
    if len(arr) == 0:
        return np.nan
    return float(np.exp(np.mean(np.log(arr))))


def pval_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _classify_tm_note(note: str) -> str:
    """Classify a single TM feature note into a category."""
    note = note or ""
    if "voltage-sensor" in note.lower() or (
        "S4" in note and ("repeat" in note.lower() or "segment" in note.lower())
    ):
        return "Voltage sensor (S4)"
    if "S5" in note:
        return "Pore-lining (S5)"
    if "S6" in note:
        return "Pore-lining (S6)"
    if any(f"S{i}" in note for i in [1, 2, 3]):
        return "VSD other (S1-S3)"
    if any(f"Name={i}" in note for i in range(1, 13)):
        return "Transporter helix"
    if note.startswith("Helical; Name=M"):
        return "Named TM (M1/M2)"
    return "Generic TM helix"


def build_residue_annotations(all_features: dict, genes_in_data: set) -> pd.DataFrame:
    """Build a DataFrame of (uniprot_id, residue_index, tm_category) from features.

    Returns rows only for residues that fall within TM/intramembrane annotations.
    Residues not in this DataFrame are Non-TM.
    """
    records = []
    for gene in genes_in_data:
        feats = all_features.get(gene, {})
        seen = set()
        for key, entries in feats.items():
            kl = key.lower()
            is_tm = "transmembrane" in kl
            is_intra = "intramembrane" in kl
            if not is_tm and not is_intra:
                continue
            for e in entries:
                s, end = e.get("start"), e.get("end")
                if s is None or end is None:
                    continue
                s, end = int(s), int(end)
                note = e.get("note", "") or ""
                cat = "Intramembrane" if is_intra else _classify_tm_note(note)
                for r in range(s, end + 1):
                    if (gene, r) not in seen:
                        seen.add((gene, r))
                        records.append((gene, r, cat))

    return pd.DataFrame(records, columns=["uniprot_id", "residue_index", "tm_category"])


# ---------------------------------------------------------------------------
# Data loading (vectorized)
# ---------------------------------------------------------------------------


def load_data():
    print("Loading constraint data ...")
    df = pd.read_csv(
        FORWARD_TSV, sep="\t", low_memory=False,
        usecols=["uniprot_id", "residue_index", "oe_upper", "is_null", "run_label", "plddt"],
    )
    print(f"  {len(df):,d} rows, {df['uniprot_id'].nunique()} genes")

    print("Loading UniProt features ...")
    with open(FEATURES_JSON) as fh:
        all_features = json.load(fh)

    genes_in_data = set(df["uniprot_id"].unique())
    print("Building residue-level TM annotations ...")
    tm_ann = build_residue_annotations(all_features, genes_in_data)
    print(f"  {len(tm_ann):,d} annotated residue entries across {tm_ann['uniprot_id'].nunique()} genes")

    df = df.merge(tm_ann, on=["uniprot_id", "residue_index"], how="left")
    df["is_tm"] = df["tm_category"].notna()
    df.loc[~df["is_tm"], "tm_category"] = "Non-TM"

    print("Loading DD de novo variants ...")
    vdf = pd.read_csv(DD_VARIANTS_TSV, sep="\t", low_memory=False)
    dd_col = "variant_level_annotations.dd_denovo_no_transcript_match.case_control"
    dd_case = (
        vdf[vdf[dd_col].str.contains("DD", na=False)]
        .dropna(subset=["uniprot_id", "residue_index"])
        .copy()
    )
    dd_case["residue_index"] = dd_case["residue_index"].astype(int)
    dd_counts = dd_case.groupby("uniprot_id").size().rename("n_dd")

    # Build per-gene metadata: gene symbol and gnomAD gene-level metrics
    gene_meta_cols = {
        "uniprot_id": "uniprot_id",
        "gene_symbol": "gene_symbol",
        "gene_level_annotations.mis.oe_ci.upper": "gnomad_oe_upper",
        "gene_level_annotations.mis.z_score": "gnomad_z",
    }
    available = [c for c in gene_meta_cols if c in vdf.columns]
    gene_meta = (
        vdf[available]
        .rename(columns=gene_meta_cols)
        .dropna(subset=["uniprot_id"])
        .drop_duplicates(subset=["uniprot_id"])
        .set_index("uniprot_id")
    )
    print(f"  Gene metadata: {len(gene_meta)} genes with symbol/gnomAD metrics")

    return df, all_features, dd_counts, gene_meta


# ---------------------------------------------------------------------------
# Figure 1: TM vs non-TM paired distributions
# ---------------------------------------------------------------------------


def _series_geomean(s):
    """Compute geomean of a Series, returning NaN if fewer than 3 positive values."""
    vals = s.dropna()
    vals = vals[vals > 0]
    if len(vals) < 3:
        return np.nan
    return float(np.exp(np.mean(np.log(vals))))


def compute_per_gene_tm_stats(df, dd_counts):
    """Compute per-gene geomean OE upper for TM vs non-TM across all methods."""
    print("Computing per-gene TM stats (vectorized) ...")

    tm_residue_counts = (
        df[df["is_tm"]]
        .drop_duplicates(subset=["uniprot_id", "residue_index"])
        .groupby("uniprot_id")
        .size()
    )
    genes_with_enough_tm = set(tm_residue_counts[tm_residue_counts >= 10].index)

    valid = filter_has_nonnull(df[df["uniprot_id"].isin(genes_with_enough_tm)].copy())
    print(f"  {len(valid):,d} rows (incl. catch-all; excl. fully-null gene-methods) for {valid['uniprot_id'].nunique()} genes")

    agg = (
        valid.groupby(["uniprot_id", "run_label", "is_tm"])["oe_upper"]
        .agg(
            geomean_oe=_series_geomean,
            n_valid="count",
        )
        .reset_index()
    )
    agg = agg[agg["n_valid"] >= 3]

    label_to_group = {}
    for group_name, aic_label, lrt_label in METHOD_GROUPS:
        label_to_group[aic_label] = (group_name, "AIC")
        label_to_group[lrt_label] = (group_name, "LRT")

    agg["group"] = agg["run_label"].map(lambda x: label_to_group.get(x, (None, None))[0])
    agg["comp"] = agg["run_label"].map(lambda x: label_to_group.get(x, (None, None))[1])
    agg = agg[agg["group"].notna()]

    tm_agg = agg[agg["is_tm"]].rename(columns={"geomean_oe": "geomean_oe_tm"})
    nontm_agg = agg[~agg["is_tm"]].rename(columns={"geomean_oe": "geomean_oe_nontm"})

    pgdf = tm_agg[["uniprot_id", "run_label", "group", "comp", "geomean_oe_tm"]].merge(
        nontm_agg[["uniprot_id", "run_label", "geomean_oe_nontm"]],
        on=["uniprot_id", "run_label"],
    )
    pgdf = pgdf.rename(columns={"uniprot_id": "gene"})
    pgdf = pgdf.merge(dd_counts.reset_index().rename(columns={"uniprot_id": "gene"}), on="gene", how="left")
    pgdf["n_dd"] = pgdf["n_dd"].fillna(0).astype(int)

    print("Computing per-gene Mann-Whitney U ...")
    mwu_rows = []
    for (gene, run_label), grp in valid.groupby(["uniprot_id", "run_label"]):
        tm_oe = grp.loc[grp["is_tm"], "oe_upper"].dropna()
        nontm_oe = grp.loc[~grp["is_tm"], "oe_upper"].dropna()
        if len(tm_oe) < 3 or len(nontm_oe) < 3:
            continue
        try:
            _, pval = mannwhitneyu(tm_oe, nontm_oe, alternative="two-sided")
        except ValueError:
            pval = np.nan
        mwu_rows.append({"gene": gene, "run_label": run_label, "mwu_pval": pval})
    mwu_df = pd.DataFrame(mwu_rows)
    pgdf = pgdf.merge(mwu_df, on=["gene", "run_label"], how="left")

    print(f"  {len(pgdf)} gene-method combinations")
    return pgdf


def figure_1_paired_distributions(pgdf):
    """Paired violin + strip showing TM vs non-TM geomean OE distributions."""
    print("\n--- Figure 1: paired TM vs non-TM distributions ---")
    group_order = [g[0] for g in METHOD_GROUPS]

    wilcox_results = {}
    print("\nPaired Wilcoxon signed-rank (TM vs non-TM geomean OE upper):")
    for group_name, _, _ in METHOD_GROUPS:
        for comp in ["AIC", "LRT"]:
            sub = pgdf[(pgdf["group"] == group_name) & (pgdf["comp"] == comp)].copy()
            sub = sub.dropna(subset=["geomean_oe_tm", "geomean_oe_nontm"])
            diffs = sub["geomean_oe_tm"] - sub["geomean_oe_nontm"]
            diffs = diffs[diffs != 0]
            key = (group_name, comp)
            if len(diffs) >= 10:
                stat, p = wilcoxon(diffs)
                direction = "TM < non-TM" if diffs.median() < 0 else "TM > non-TM"
                n_sig = (sub["mwu_pval"] < 0.05).sum()
                n_tm_lower = (sub["geomean_oe_tm"] < sub["geomean_oe_nontm"]).sum()
                wilcox_results[key] = {
                    "p": p, "median_diff": diffs.median(), "direction": direction,
                    "n": len(diffs), "n_sig": int(n_sig), "n_tm_lower": int(n_tm_lower),
                    "n_total": len(sub),
                }
                print(
                    f"  {group_name:25s} ({comp}): p={p:.2e}, "
                    f"median diff={diffs.median():.4f}, {direction}, "
                    f"n={len(diffs)}, {n_sig} sig (p<0.05), "
                    f"{n_tm_lower}/{len(sub)} TM<non-TM"
                )
            else:
                wilcox_results[key] = None

    n_groups = len(group_order)
    fig, axes = plt.subplots(2, n_groups, figsize=(6 * n_groups, 12), sharey=True)
    fig.suptitle(
        "Distribution of per-gene geometric mean OE upper:\n"
        "TM vs non-TM residues across methods",
        fontsize=15, fontweight="bold", y=1.02,
    )

    rng = np.random.default_rng(42)
    for gi, gname in enumerate(group_order):
        for ci, comp in enumerate(["AIC", "LRT"]):
            ax = axes[ci, gi]
            csub = pgdf[(pgdf["group"] == gname) & (pgdf["comp"] == comp)].dropna(
                subset=["geomean_oe_tm", "geomean_oe_nontm"]
            )
            if len(csub) == 0:
                ax.set_visible(False)
                continue

            tm_vals = csub["geomean_oe_tm"].values
            nontm_vals = csub["geomean_oe_nontm"].values

            parts = ax.violinplot(
                [tm_vals, nontm_vals], positions=[0, 1],
                showmedians=False, showextrema=False, widths=0.7,
            )
            for i, body in enumerate(parts["bodies"]):
                body.set_facecolor(C_TM if i == 0 else C_NONTM)
                body.set_alpha(0.3)

            for i, (vals, color) in enumerate([(tm_vals, C_TM), (nontm_vals, C_NONTM)]):
                jitter = rng.normal(0, 0.06, len(vals))
                ax.scatter(
                    np.full_like(vals, i) + jitter, vals, c=color, s=8,
                    alpha=0.3, edgecolors="none", zorder=3,
                )
                ax.plot(
                    [i - 0.2, i + 0.2], [np.median(vals)] * 2,
                    color="black", linewidth=2.5, zorder=5,
                )

            key = (gname, comp)
            wr = wilcox_results.get(key)
            if wr:
                stars = pval_stars(wr["p"])
                ann = (
                    f"Wilcoxon p={wr['p']:.1e} {stars}\n"
                    f"med diff={wr['median_diff']:.3f}\n"
                    f"n={wr['n']}, {wr['n_tm_lower']}/{wr['n_total']} TM<nonTM"
                )
                ax.text(
                    0.5, 0.98, ann, transform=ax.transAxes, fontsize=8, ha="center",
                    va="top", bbox=dict(facecolor="white", edgecolor="gray", alpha=0.8),
                )

            color = C_AIC if comp == "AIC" else C_LRT
            ax.set_title(f"{gname}\n({comp})", fontsize=11, fontweight="bold", color=color)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["TM", "Non-TM"], fontsize=10)
            if gi == 0:
                ax.set_ylabel("Geometric mean OE upper", fontsize=11)
            ax.axhline(0.6, color="red", linestyle="--", alpha=0.3, linewidth=0.8)
            ax.axhline(1.0, color="gray", linestyle=":", alpha=0.3, linewidth=0.8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_tm_vs_nontm_distributions.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig1_tm_vs_nontm_distributions.png")
    return wilcox_results


# ---------------------------------------------------------------------------
# Figures 2–4: TM segment type breakdown (vectorized)
# ---------------------------------------------------------------------------


def figures_2_3_4_segment_types(df):
    print("\n--- Figures 2-4: TM segment type breakdown ---")
    METHOD = "AIC (no pLDDT/PAE)"
    msub = filter_has_nonnull(df[df["run_label"] == METHOD].copy())
    print(f"  Working with {len(msub):,d} rows (incl. catch-all; excl. fully-null gene-methods)")

    msub["log_oe"] = np.log(msub["oe_upper"].clip(lower=1e-10))

    gene_cat_stats = (
        msub.groupby(["uniprot_id", "tm_category"])
        .agg(
            n_valid=("oe_upper", "count"),
            mean_log_oe=("log_oe", "mean"),
        )
        .reset_index()
    )
    gene_cat_stats = gene_cat_stats[gene_cat_stats["n_valid"] >= 2]
    gene_cat_stats["geomean_oe"] = np.exp(gene_cat_stats["mean_log_oe"])

    genes_with_nontm = set(
        gene_cat_stats[
            (gene_cat_stats["tm_category"] == "Non-TM") & (gene_cat_stats["n_valid"] >= 3)
        ]["uniprot_id"]
    )
    cdf = gene_cat_stats[gene_cat_stats["uniprot_id"].isin(genes_with_nontm)].copy()
    cdf = cdf.rename(columns={"uniprot_id": "gene", "tm_category": "category"})

    print(f"  {cdf['gene'].nunique()} genes with both TM and non-TM data")

    # Per-gene MWU tests (vectorized by pre-grouping)
    print("  Computing Mann-Whitney U per gene/category ...")
    nontm_data = msub[msub["tm_category"] == "Non-TM"].groupby("uniprot_id")["oe_upper"].apply(list).to_dict()
    cat_data = msub[msub["tm_category"] != "Non-TM"].groupby(["uniprot_id", "tm_category"])["oe_upper"].apply(list).to_dict()

    pval_records = []
    for (gene, cat), cat_vals in cat_data.items():
        nontm_vals = nontm_data.get(gene)
        if nontm_vals is None or len(nontm_vals) < 3 or len(cat_vals) < 3:
            continue
        try:
            _, pval = mannwhitneyu(cat_vals, nontm_vals, alternative="two-sided")
        except ValueError:
            pval = np.nan
        pval_records.append({"gene": gene, "category": cat, "mwu_pval": pval})

    if pval_records:
        pval_df = pd.DataFrame(pval_records)
        cdf = cdf.merge(pval_df, on=["gene", "category"], how="left")
    else:
        cdf["mwu_pval"] = np.nan

    # ----- Figure 2: Ion channels -----
    ion_channel_cats = [
        "Voltage sensor (S4)", "VSD other (S1-S3)",
        "Pore-lining (S5)", "Pore-lining (S6)", "Non-TM",
    ]
    ion_genes = cdf[cdf["category"] == "Voltage sensor (S4)"]["gene"].unique()
    if len(ion_genes) > 0:
        ion_sort = cdf[
            (cdf["gene"].isin(ion_genes)) & (cdf["category"] == "Voltage sensor (S4)")
        ]
        ion_gene_order = ion_sort.sort_values("geomean_oe", ascending=False)["gene"].tolist()

        fig, ax = plt.subplots(figsize=(14, max(5, len(ion_gene_order) * 0.45)))
        fig.suptitle(
            f"Ion channel TM segment constraint by type ({len(ion_gene_order)} genes)\n"
            "(AIC, no filter; geometric mean OE upper)",
            fontsize=14, fontweight="bold",
        )

        for yi, gene in enumerate(ion_gene_order):
            gsub = cdf[cdf["gene"] == gene]
            vals = {row["category"]: row["geomean_oe"] for _, row in gsub.iterrows()}
            all_vals = [v for v in vals.values() if np.isfinite(v)]
            if all_vals:
                ax.plot(
                    [min(all_vals), max(all_vals)], [yi, yi],
                    color="#cccccc", linewidth=0.8, zorder=0,
                )
            for cat in ion_channel_cats:
                if cat not in vals or not np.isfinite(vals[cat]):
                    continue
                prow = gsub[gsub["category"] == cat]
                if len(prow) == 0:
                    continue
                prow = prow.iloc[0]
                pval = prow.get("mwu_pval", np.nan) if cat != "Non-TM" else np.nan
                ec = "black" if (not np.isnan(pval) and pval < 0.05) else CAT_COLORS[cat]
                lw = 1.5 if (not np.isnan(pval) and pval < 0.05) else 0.5
                sz = 80 if cat != "Non-TM" else 60
                ax.scatter(
                    vals[cat], yi, color=CAT_COLORS[cat], s=sz, zorder=2,
                    edgecolors=ec, linewidth=lw, marker=CAT_MARKERS[cat],
                )

        ax.set_yticks(range(len(ion_gene_order)))
        ax.set_yticklabels(ion_gene_order, fontsize=max(5, min(10, 500 // len(ion_gene_order))))
        ax.set_xlabel("Geometric mean OE upper (lower = more constrained)", fontsize=12)
        ax.axvline(0.6, color="red", linestyle="--", alpha=0.3, linewidth=0.8)
        ax.axvline(1.0, color="gray", linestyle=":", alpha=0.3, linewidth=0.8)
        legend_els = [
            Line2D([0], [0], marker=CAT_MARKERS[c], color="w", markerfacecolor=CAT_COLORS[c],
                   markeredgecolor="black", markersize=9, label=c)
            for c in ion_channel_cats
        ]
        legend_els.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                   markeredgecolor="black", markersize=9, linewidth=1.5,
                   label="p<0.05 vs non-TM (black edge)")
        )
        ax.legend(handles=legend_els, fontsize=9, loc="upper right", framealpha=0.9, edgecolor="gray")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "fig2_ion_channel_segment_types.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved fig2_ion_channel_segment_types.png ({len(ion_gene_order)} ion channel genes)")
    else:
        print("  No ion channel genes found, skipping Figure 2")

    # ----- Figure 3: Histograms by category -----
    all_cats = [
        "Non-TM", "Generic TM helix", "Voltage sensor (S4)", "VSD other (S1-S3)",
        "Pore-lining (S5)", "Pore-lining (S6)", "Transporter helix",
        "Named TM (M1/M2)", "Intramembrane",
    ]
    # Sort categories by median geomean OE descending (highest = least constrained first)
    cat_medians = {
        c: cdf[cdf["category"] == c]["geomean_oe"].median()
        for c in all_cats if c in cdf["category"].values
    }
    present_cats = sorted(cat_medians, key=cat_medians.get, reverse=True)
    n_cats = len(present_cats)

    fig, axes = plt.subplots(n_cats, 1, figsize=(12, 2.2 * n_cats), sharex=True)
    if n_cats == 1:
        axes = [axes]
    fig.suptitle(
        "Distribution of per-gene geomean OE upper by TM segment type\n"
        f"(AIC, no filter; {cdf['gene'].nunique()} genes)",
        fontsize=14, fontweight="bold",
    )

    for i, cat in enumerate(present_cats):
        ax = axes[i]
        vals = cdf[cdf["category"] == cat]["geomean_oe"].dropna().values
        ax.hist(vals, bins=50, color=CAT_COLORS[cat], alpha=0.7, edgecolor="white", linewidth=0.3)
        med = np.median(vals)
        ax.axvline(med, color="black", linewidth=1.5, linestyle="-", zorder=5)
        ax.axvline(0.6, color="red", linestyle="--", alpha=0.4, linewidth=0.8)
        ax.axvline(1.0, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
        ax.set_ylabel(f"{cat}\n(n={len(vals)})", fontsize=9, rotation=0, ha="right", va="center")
        ax.text(
            0.98, 0.85, f"median={med:.3f}",
            transform=ax.transAxes, fontsize=9, ha="right", va="top",
            bbox=dict(facecolor="white", edgecolor="gray", alpha=0.7),
        )
        ax.set_yticks([])

    axes[-1].set_xlabel("Geometric mean OE upper (lower = more constrained)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_segment_type_histograms.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig3_segment_type_histograms.png")

    # ----- Figure 4: Aggregate strip -----
    plot_cats = [c for c in present_cats if len(cdf[cdf["category"] == c]) >= 2]
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.suptitle(
        "Constraint by TM segment type — aggregate\n"
        "(geometric mean OE upper per gene; each dot = one gene)",
        fontsize=13, fontweight="bold",
    )

    rng = np.random.default_rng(42)
    for xi, cat in enumerate(plot_cats):
        vals = cdf[cdf["category"] == cat]["geomean_oe"].dropna().values
        jitter = rng.normal(0, 0.12, len(vals))
        ax.scatter(
            np.full_like(vals, xi) + jitter, vals, color=CAT_COLORS[cat],
            s=15, alpha=0.4, edgecolors="none", marker=CAT_MARKERS[cat],
        )
        ax.plot([xi - 0.3, xi + 0.3], [np.median(vals)] * 2, color="black", linewidth=2, zorder=5)
        q25, q75 = np.percentile(vals, [25, 75])
        ax.plot([xi, xi], [q25, q75], color="black", linewidth=1.5, zorder=4)

    ax.set_xticks(range(len(plot_cats)))
    ax.set_xticklabels(plot_cats, fontsize=9, rotation=30, ha="right")
    ax.set_ylabel("Geometric mean OE upper\n(lower = more constrained)", fontsize=11)
    ax.axhline(0.6, color="red", linestyle="--", alpha=0.3, linewidth=0.8)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.3, linewidth=0.8)
    ax.text(
        0.01, 0.98, "Horizontal bar = median\nVertical bar = IQR",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(facecolor="white", edgecolor="gray", alpha=0.8),
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_category_summary_strip.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig4_category_summary_strip.png")

    # Print summary table
    print("\n  Category summary:")
    for cat in present_cats:
        sub = cdf[cdf["category"] == cat]
        n = len(sub)
        med = sub["geomean_oe"].median()
        if cat != "Non-TM" and "mwu_pval" in sub.columns:
            n_sig = (sub["mwu_pval"] < 0.05).sum()
            print(f"    {cat:30s}: n={n:5d}, median geomean={med:.3f}, {n_sig}/{n} sig vs non-TM")
        else:
            print(f"    {cat:30s}: n={n:5d}, median geomean={med:.3f}")

    return cdf


# ---------------------------------------------------------------------------
# Figure 5: pLDDT bias check (vectorized)
# ---------------------------------------------------------------------------


def figure_5_plddt_bias(df):
    print("\n--- Figure 5: pLDDT bias check ---")
    METHOD = "AIC (no pLDDT/PAE)"
    msub = df[df["run_label"] == METHOD].copy()

    tm_residue_counts = (
        msub[msub["is_tm"]]
        .drop_duplicates(subset=["uniprot_id", "residue_index"])
        .groupby("uniprot_id")
        .size()
    )
    genes_with_enough_tm = set(tm_residue_counts[tm_residue_counts >= 10].index)
    msub = msub[msub["uniprot_id"].isin(genes_with_enough_tm)]

    plddt_valid = msub.dropna(subset=["plddt"])

    tm_plddt = plddt_valid[plddt_valid["is_tm"]]["plddt"].values
    nontm_plddt = plddt_valid[~plddt_valid["is_tm"]]["plddt"].values

    per_gene = (
        plddt_valid.groupby(["uniprot_id", "is_tm"])
        .agg(
            median_plddt=("plddt", "median"),
            pct_below_70=("plddt", lambda x: (x < 70).mean() * 100),
            n=("plddt", "count"),
        )
        .reset_index()
    )
    per_gene_tm = per_gene[per_gene["is_tm"]].set_index("uniprot_id")
    per_gene_nontm = per_gene[~per_gene["is_tm"]].set_index("uniprot_id")
    common_genes = per_gene_tm.index.intersection(per_gene_nontm.index)
    common_genes = common_genes[
        (per_gene_tm.loc[common_genes, "n"] >= 3) & (per_gene_nontm.loc[common_genes, "n"] >= 3)
    ]

    diffs = per_gene_tm.loc[common_genes, "median_plddt"] - per_gene_nontm.loc[common_genes, "median_plddt"]
    pct_diff = per_gene_nontm.loc[common_genes, "pct_below_70"] - per_gene_tm.loc[common_genes, "pct_below_70"]

    n_genes = len(common_genes)
    n_tm_higher = (diffs > 0).sum()
    print(f"  TM:     median pLDDT={np.median(tm_plddt):.1f}, %<70: {(tm_plddt < 70).mean() * 100:.1f}%  ({len(tm_plddt):,d} residues)")
    print(f"  Non-TM: median pLDDT={np.median(nontm_plddt):.1f}, %<70: {(nontm_plddt < 70).mean() * 100:.1f}%  ({len(nontm_plddt):,d} residues)")
    print(f"  {n_tm_higher}/{n_genes} genes have higher median pLDDT in TM")

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle(f"pLDDT bias check: TM vs non-TM residues ({n_genes} genes)", fontsize=14, fontweight="bold")

    ax = axes[0]
    ax.hist(nontm_plddt, bins=60, alpha=0.5, color=C_NONTM, density=True, label=f"Non-TM (n={len(nontm_plddt):,d})")
    ax.hist(tm_plddt, bins=60, alpha=0.5, color=C_TM, density=True, label=f"TM (n={len(tm_plddt):,d})")
    ax.axvline(70, color="red", linestyle="--", linewidth=1.5, label="pLDDT=70 cutoff")
    ax.set_xlabel("pLDDT")
    ax.set_ylabel("Density")
    ax.set_title("pLDDT distributions (all residues pooled)")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.hist(diffs.values, bins=50, color="#7570b3", alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linestyle="-", linewidth=1)
    ax.axvline(np.median(diffs), color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Median pLDDT difference (TM − non-TM)")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Per-gene pLDDT difference\n(median={np.median(diffs):.1f}, {n_tm_higher}/{n_genes} TM higher)")

    ax = axes[2]
    ax.hist(pct_diff.values, bins=50, color="#e7298a", alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linestyle="-", linewidth=1)
    ax.axvline(np.median(pct_diff), color="red", linestyle="--", linewidth=1.5)
    n_nontm_more_excluded = (pct_diff > 0).sum()
    ax.set_xlabel("% excluded difference (non-TM − TM)")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Differential pLDDT exclusion\n(median={np.median(pct_diff):.1f}pp, {n_nontm_more_excluded}/{n_genes} lose more non-TM)")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_plddt_bias_check.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig5_plddt_bias_check.png")


# ---------------------------------------------------------------------------
# Per-category Wilcoxon table (broken down by TM segment type × method)
# ---------------------------------------------------------------------------


CAT_ORDER = [
    "Pore-lining (S5)", "Pore-lining (S6)", "Named TM (M1/M2)",
    "Voltage sensor (S4)", "VSD other (S1-S3)", "Intramembrane",
    "Transporter helix", "Generic TM helix",
]


def compute_category_wilcoxon_table(df):
    """Paired Wilcoxon signed-rank test: each TM category vs Non-TM, per method."""
    print("\n--- Category × Method Wilcoxon table ---")

    valid = filter_has_nonnull(df.copy())

    label_to_method = {}
    for group_name, aic_label, lrt_label in METHOD_GROUPS:
        label_to_method[aic_label] = f"{group_name} (AIC)"
        label_to_method[lrt_label] = f"{group_name} (LRT)"

    valid["method"] = valid["run_label"].map(label_to_method)
    valid = valid[valid["method"].notna()]

    print("  Computing per-gene geomean by category ...")
    cat_agg = (
        valid.groupby(["uniprot_id", "method", "tm_category"])["oe_upper"]
        .agg(geomean_oe=_series_geomean, n_valid="count")
        .reset_index()
    )
    cat_agg = cat_agg[cat_agg["n_valid"] >= 2]

    nontm = cat_agg[cat_agg["tm_category"] == "Non-TM"][
        ["uniprot_id", "method", "geomean_oe"]
    ].rename(columns={"geomean_oe": "geomean_nontm"})

    rows = []
    method_order = []
    for group_name, aic_label, lrt_label in METHOD_GROUPS:
        method_order.append(f"{group_name} (AIC)")
        method_order.append(f"{group_name} (LRT)")

    for cat in CAT_ORDER:
        cat_sub = cat_agg[cat_agg["tm_category"] == cat][
            ["uniprot_id", "method", "geomean_oe"]
        ].rename(columns={"geomean_oe": "geomean_cat"})

        paired = cat_sub.merge(nontm, on=["uniprot_id", "method"])
        paired = paired.dropna(subset=["geomean_cat", "geomean_nontm"])

        for method in method_order:
            msub = paired[paired["method"] == method]
            diffs = msub["geomean_cat"] - msub["geomean_nontm"]
            diffs_nz = diffs[diffs != 0]
            n = len(diffs_nz)
            n_total = len(msub)
            if n >= 10:
                stat, p = wilcoxon(diffs_nz)
                med_diff = diffs.median()
                direction = "cat < non-TM" if med_diff < 0 else "cat ≥ non-TM"
                n_cat_lower = (diffs < 0).sum()
            else:
                p = np.nan
                med_diff = diffs.median() if n_total > 0 else np.nan
                direction = "—"
                n_cat_lower = (diffs < 0).sum() if n_total > 0 else 0

            rows.append({
                "TM category": cat,
                "Method": method,
                "n genes": n_total,
                "Median diff": med_diff,
                "Direction": direction,
                "n cat<nonTM": int(n_cat_lower),
                "% cat<nonTM": f"{n_cat_lower / n_total * 100:.0f}%" if n_total > 0 else "—",
                "Wilcoxon p": p,
                "Sig": pval_stars(p) if not np.isnan(p) else "",
            })

    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "category_wilcoxon_table.tsv", sep="\t", index=False)
    print("  Saved category_wilcoxon_table.tsv")

    print("\n  Results:")
    for cat in CAT_ORDER:
        csub = table[table["TM category"] == cat]
        print(f"\n  {cat}:")
        for _, r in csub.iterrows():
            p_str = f"{r['Wilcoxon p']:.1e}" if not np.isnan(r["Wilcoxon p"]) else "n/a"
            print(
                f"    {r['Method']:30s}  n={r['n genes']:5d}  "
                f"med_diff={r['Median diff']:+.4f}  {r['% cat<nonTM']:>4s} lower  "
                f"p={p_str} {r['Sig']}"
            )

    return table


# ---------------------------------------------------------------------------
# Figure 6: Bimodality analysis
# ---------------------------------------------------------------------------

# Channel subfamily prefixes used to annotate genes in the scatter plot.
# Keys are UniProt gene-name prefixes; values are short labels.
_CHANNEL_FAMILY = {
    "SCN": "Nav", "CACNA": "Cav", "KCNA": "Kv1", "KCNB": "Kv2",
    "KCNC": "Kv3", "KCND": "Kv4", "KCNQ": "Kv7/KCNQ", "KCNH": "Kv10-12/EAG",
    "KCNF": "Kv", "KCNG": "Kv", "KCNS": "Kv", "KCNV": "Kv", "HCN": "HCN",
    "CNGA": "CNG", "CNGB": "CNG", "TRPA": "TRP", "TRPC": "TRP",
    "TRPM": "TRP", "TRPV": "TRP", "TRPP": "TRP", "TRPML": "TRP",
    "RYR": "Ryanodine", "ITPR": "IP3R",
}

_FAMILY_COLORS = {
    "Nav": "#d73027", "Cav": "#f46d43", "Kv1": "#4575b4", "Kv2": "#74add1",
    "Kv3": "#abd9e9", "Kv4": "#e0f3f8", "Kv7/KCNQ": "#313695",
    "Kv10-12/EAG": "#fee090", "Kv": "#9ecae1", "HCN": "#41ab5d",
    "CNG": "#addd8e", "TRP": "#9970ab", "Ryanodine": "#bf812d",
    "IP3R": "#dfc27d", "Other": "#aaaaaa",
}


def _gene_family(gene_symbol: str) -> str:
    for prefix, fam in _CHANNEL_FAMILY.items():
        if gene_symbol.upper().startswith(prefix):
            return fam
    return "Other"


def figure_6_bimodality_analysis(cdf, dd_counts, gene_meta):
    """Three-panel investigation of bimodality in constrained TM categories."""
    print("\n--- Figure 6: Bimodality analysis ---")

    # Ion channel genes = those with S4 annotation in cdf
    ion_genes = set(cdf[cdf["category"] == "Voltage sensor (S4)"]["gene"].unique())
    print(f"  {len(ion_genes)} ion channel genes with S4 annotation")

    # Wide table: one row per gene, columns = categories
    target_cats = ["Pore-lining (S5)", "Voltage sensor (S4)", "Non-TM"]
    ion_cdf = cdf[cdf["gene"].isin(ion_genes) & cdf["category"].isin(target_cats)].copy()
    ion_wide = ion_cdf.pivot_table(index="gene", columns="category", values="geomean_oe")

    # Attach DD counts and gene symbols
    ion_wide = ion_wide.join(dd_counts.rename("n_dd"), how="left")
    ion_wide["n_dd"] = ion_wide["n_dd"].fillna(0).astype(int)
    # Use gene_meta for symbols; fall back to uniprot_id if not found
    gene_symbol_map = gene_meta["gene_symbol"].to_dict() if "gene_symbol" in gene_meta.columns else {}
    ion_wide["symbol"] = ion_wide.index.map(lambda u: gene_symbol_map.get(u, u))
    ion_wide["family"] = ion_wide["symbol"].map(_gene_family)

    dd_bins = [-1, 0, 5, 1_000]
    dd_labels = ["0 (none)", "1–5", ">5"]
    ion_wide["dd_group"] = pd.cut(ion_wide["n_dd"], bins=dd_bins, labels=dd_labels)
    dd_colors = {"0 (none)": "#aaaaaa", "1–5": "#f4a582", ">5": "#d6604d"}

    s5_col = "Pore-lining (S5)"
    nt_col = "Non-TM"
    s4_col = "Voltage sensor (S4)"

    # ---- print summary stats ----
    print(f"\n  Breakdown by DD de novo burden:")
    for dd_grp in dd_labels:
        sub = ion_wide[ion_wide["dd_group"] == dd_grp]
        med_s5 = sub[s5_col].median() if s5_col in sub.columns else float("nan")
        med_nt = sub[nt_col].median() if nt_col in sub.columns else float("nan")
        print(f"    {dd_grp:10s}: n={len(sub):3d}, "
              f"median S5={med_s5:.3f}, median nonTM={med_nt:.3f}")

    if s5_col in ion_wide.columns:
        print("\n  Top-10 most constrained by S5 geomean OE:")
        top10 = ion_wide.dropna(subset=[s5_col]).nsmallest(10, s5_col)
        for gene, row in top10.iterrows():
            print(f"    {row['symbol']:12s} ({gene}): S5={row[s5_col]:.3f}, "
                  f"non-TM={row[nt_col]:.3f}, DD={int(row['n_dd'])}")

    # ---- Figure ----
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    fig.suptitle(
        "Bimodality in ion channel TM constraint: DD burden, subfamily, and whole-gene effects\n"
        f"(AIC, no filter; {len(ion_genes)} genes with voltage-sensor S1–S6 segments)",
        fontsize=13, fontweight="bold",
    )

    # --- Panel A: S5 OE histogram split by DD de novo burden ---
    ax = axes[0]
    if s5_col in ion_wide.columns:
        for dd_grp in [">5", "1–5", "0 (none)"]:
            vals = ion_wide[ion_wide["dd_group"] == dd_grp][s5_col].dropna().values
            if len(vals):
                ax.hist(vals, bins=18, alpha=0.65, color=dd_colors[dd_grp],
                        label=f"{dd_grp} DD de novo (n={len(vals)})", density=True)
    ax.set_xlabel("Geomean OE upper — Pore-lining (S5)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title("S5 constraint by DD de novo burden\n(do high-DD genes dominate the low-OE mode?)")
    ax.legend(fontsize=9)
    ax.axvline(0.6, color="red", linestyle="--", alpha=0.4, linewidth=0.8)
    ax.axvline(1.0, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)

    # --- Panel B: S5 OE vs Non-TM OE scatter, colored by channel family ---
    ax = axes[1]
    if s5_col in ion_wide.columns and nt_col in ion_wide.columns:
        scatter_data = ion_wide.dropna(subset=[s5_col, nt_col]).copy()
        families_present = scatter_data["family"].unique()
        for fam in sorted(families_present, key=lambda f: (f == "Other", f)):
            sub = scatter_data[scatter_data["family"] == fam]
            ax.scatter(
                sub[s5_col], sub[nt_col],
                c=_FAMILY_COLORS.get(fam, "#aaaaaa"),
                s=45, alpha=0.75, edgecolors="gray", linewidth=0.4, zorder=3,
                label=fam,
            )
        # Identity line
        lmin = min(scatter_data[s5_col].min(), scatter_data[nt_col].min()) * 0.92
        lmax = max(scatter_data[s5_col].max(), scatter_data[nt_col].max()) * 1.05
        ax.plot([lmin, lmax], [lmin, lmax], "k--", alpha=0.3, linewidth=1)
        # Label most constrained genes
        for gene, row in scatter_data.nsmallest(8, s5_col).iterrows():
            ax.annotate(
                row["symbol"], (row[s5_col], row[nt_col]),
                fontsize=6.5, xytext=(3, 2), textcoords="offset points",
                ha="left", va="bottom",
            )
        ax.set_xlabel("Geomean OE — Pore-lining (S5)", fontsize=11)
        ax.set_ylabel("Geomean OE — Non-TM", fontsize=11)
        ax.set_title(
            "S5 vs. Non-TM OE per gene, by channel subfamily\n"
            "(points above diagonal: non-TM more constrained than S5)"
        )
        ax.legend(fontsize=7, ncol=2, loc="upper left", framealpha=0.85)
        ax.axhline(1.0, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
        ax.axvline(0.6, color="red", linestyle="--", alpha=0.4, linewidth=0.8)

    # --- Panel C: Non-TM OE histogram for the same genes (is non-TM also bimodal?) ---
    ax = axes[2]
    if nt_col in ion_wide.columns:
        for dd_grp in [">5", "1–5", "0 (none)"]:
            vals = ion_wide[ion_wide["dd_group"] == dd_grp][nt_col].dropna().values
            if len(vals):
                ax.hist(vals, bins=18, alpha=0.65, color=dd_colors[dd_grp],
                        label=f"{dd_grp} DD de novo (n={len(vals)})", density=True)
    ax.set_xlabel("Geomean OE upper — Non-TM", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(
        "Non-TM OE for the same ion channel genes\n"
        "(is bimodality whole-gene or segment-specific?)"
    )
    ax.legend(fontsize=9)
    ax.axvline(1.0, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig6_bimodality_analysis.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig6_bimodality_analysis.png")

    return ion_wide


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    df, all_features, dd_counts, gene_meta = load_data()

    pgdf = compute_per_gene_tm_stats(df, dd_counts)
    figure_1_paired_distributions(pgdf)
    cdf = figures_2_3_4_segment_types(df)
    figure_5_plddt_bias(df)
    wilcox_table = compute_category_wilcoxon_table(df)
    figure_6_bimodality_analysis(cdf, dd_counts, gene_meta)

    print("\nAll plots saved to", OUT_DIR)


if __name__ == "__main__":
    main()
