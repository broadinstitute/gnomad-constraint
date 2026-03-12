"""Plot missense o/e alongside matched pLoF o/e by score percentile.

Shows whether missense prediction scores capture constraint signal
beyond what's explained by gene-level pLoF constraint alone.

Two weighting approaches are compared:
  1. Weight by adjusted expected missense count (expected_variants * adj_r)
  2. Weight by possible missense variants
Observed missense is NOT used as a weight (circular).
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
mis_df = pd.read_csv("/tmp/missense_oe_12_26_25.tsv", sep="\t")
plof_df = pd.read_csv("/tmp/matched_plof_all_genes.mane_select.tsv", sep="\t")

# ---------------------------------------------------------------------------
# Score definitions
# ---------------------------------------------------------------------------
score_defs = [
    ("ESM1v", "ESM_1v_neg"),
    ("PopEVE", "popEVE_neg"),
    ("AlphaMissense", "AM"),
]

colors = {
    "ESM1v": "#d62728",       # red
    "PopEVE": "#1f77b4",      # blue
    "AlphaMissense": "#2ca02c",  # green
}

# ---------------------------------------------------------------------------
# Compute o/e ratios
# ---------------------------------------------------------------------------
percentiles = mis_df["percentile"].values  # 1..100

missense_oe = {}
plof_oe_by_exp = {}  # weighted by adjusted expected missense count
plof_oe_by_pos = {}  # weighted by possible missense variants

for display, col_prefix in score_defs:
    missense_oe[display] = (
        mis_df[f"{col_prefix}.obs"].values
        / mis_df[f"{col_prefix}.exp_with_adj_r"].values
    )
    plof_oe_by_exp[display] = (
        plof_df[f"{col_prefix}.plof_obs_x_exp"].values
        / plof_df[f"{col_prefix}.plof_exp_x_exp"].values
    )
    plof_oe_by_pos[display] = (
        plof_df[f"{col_prefix}.plof_obs_x_pos"].values
        / plof_df[f"{col_prefix}.plof_exp_x_pos"].values
    )

# ---------------------------------------------------------------------------
# Reference lines
# ---------------------------------------------------------------------------
SYN_OE = 1.05
PLOF_OE = 0.58

# ---------------------------------------------------------------------------
# Shared style
# ---------------------------------------------------------------------------
OUT = "/Users/jgoodric/PycharmProjects/gnomad-constraint/proemis3d_enrichment_plots"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})


def _style_ax(ax):
    ax.axhline(SYN_OE, color="gray", linestyle="--", linewidth=1.2, zorder=1)
    ax.text(52, SYN_OE + 0.02, "obs/exp synonymous", fontsize=10,
            fontstyle="italic", color="gray", ha="center")
    ax.axhline(PLOF_OE, color="#d62728", linestyle="--", linewidth=1.2, zorder=1)
    ax.text(52, PLOF_OE + 0.02, "obs/exp pLoFs", fontsize=10,
            fontstyle="italic", color="#d62728", ha="center")
    ax.set_xlabel("Score Percentile (high is more pathogenic)", fontsize=12)
    ax.set_ylabel("Observed / Expected Ratio", fontsize=12)
    ax.set_xlim(0, 100)
    ax.set_ylim(0.25, 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ===================================================================
# Plot 1: Unmatched — missense o/e only (replicates the reference image)
# ===================================================================
fig1, ax1 = plt.subplots(figsize=(8, 6))

for display, _ in score_defs:
    ax1.plot(percentiles, missense_oe[display], color=colors[display],
             linewidth=1.8, label=display, zorder=3)

_style_ax(ax1)
legend1 = ax1.legend(loc="lower left", fontsize=10, frameon=True,
                     title="Method", title_fontsize=11)
legend1.get_frame().set_edgecolor("lightgray")

fig1.tight_layout()
fig1.savefig(f"{OUT}/missense_oe_unmatched_plot.png", dpi=200, bbox_inches="tight")
fig1.savefig(f"{OUT}/missense_oe_unmatched_plot.pdf", bbox_inches="tight")
print(f"Saved unmatched plot to {OUT}/missense_oe_unmatched_plot.png")

# ===================================================================
# Plot 2: Matched — weighted by adjusted expected missense count
# ===================================================================
fig2, ax2 = plt.subplots(figsize=(8, 6))

for display, _ in score_defs:
    c = colors[display]
    ax2.plot(percentiles, missense_oe[display], color=c, linewidth=1.8,
             label=f"{display} — missense o/e", zorder=3)
    ax2.plot(percentiles, plof_oe_by_exp[display], color=c, linewidth=1.8,
             linestyle="--", alpha=0.8,
             label=f"{display} — matched pLoF o/e", zorder=3)

_style_ax(ax2)
ax2.set_title("Matched pLoF weighted by expected missense count", fontsize=13)
legend2 = ax2.legend(loc="lower left", fontsize=9, frameon=True,
                     title="Method", title_fontsize=10, ncol=1)
legend2.get_frame().set_edgecolor("lightgray")

fig2.tight_layout()
fig2.savefig(f"{OUT}/matched_plof_oe_by_expected.png", dpi=200, bbox_inches="tight")
fig2.savefig(f"{OUT}/matched_plof_oe_by_expected.pdf", bbox_inches="tight")
print(f"Saved expected-weighted plot to {OUT}/matched_plof_oe_by_expected.png")

# ===================================================================
# Plot 3: Matched — weighted by possible missense variants
# ===================================================================
fig3, ax3 = plt.subplots(figsize=(8, 6))

for display, _ in score_defs:
    c = colors[display]
    ax3.plot(percentiles, missense_oe[display], color=c, linewidth=1.8,
             label=f"{display} — missense o/e", zorder=3)
    ax3.plot(percentiles, plof_oe_by_pos[display], color=c, linewidth=1.8,
             linestyle="--", alpha=0.8,
             label=f"{display} — matched pLoF o/e", zorder=3)

_style_ax(ax3)
ax3.set_title("Matched pLoF weighted by possible missense variants", fontsize=13)
legend3 = ax3.legend(loc="lower left", fontsize=9, frameon=True,
                     title="Method", title_fontsize=10, ncol=1)
legend3.get_frame().set_edgecolor("lightgray")

fig3.tight_layout()
fig3.savefig(f"{OUT}/matched_plof_oe_by_possible.png", dpi=200, bbox_inches="tight")
fig3.savefig(f"{OUT}/matched_plof_oe_by_possible.pdf", bbox_inches="tight")
print(f"Saved possible-weighted plot to {OUT}/matched_plof_oe_by_possible.png")

# ===================================================================
# Plot 4: Side-by-side comparison of both weighting approaches
# ===================================================================
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

for display, _ in score_defs:
    c = colors[display]
    # Left panel: weighted by expected
    ax4a.plot(percentiles, missense_oe[display], color=c, linewidth=1.8,
              label=f"{display} — missense o/e", zorder=3)
    ax4a.plot(percentiles, plof_oe_by_exp[display], color=c, linewidth=1.8,
              linestyle="--", alpha=0.8,
              label=f"{display} — matched pLoF o/e", zorder=3)
    # Right panel: weighted by possible
    ax4b.plot(percentiles, missense_oe[display], color=c, linewidth=1.8,
              label=f"{display} — missense o/e", zorder=3)
    ax4b.plot(percentiles, plof_oe_by_pos[display], color=c, linewidth=1.8,
              linestyle="--", alpha=0.8,
              label=f"{display} — matched pLoF o/e", zorder=3)

_style_ax(ax4a)
_style_ax(ax4b)
ax4a.set_title("Weight: expected missense (adj_r)", fontsize=12)
ax4b.set_title("Weight: possible missense variants", fontsize=12)
ax4b.set_ylabel("")

legend4 = ax4a.legend(loc="lower left", fontsize=8, frameon=True,
                      title="Method", title_fontsize=9, ncol=1)
legend4.get_frame().set_edgecolor("lightgray")

fig4.tight_layout()
fig4.savefig(f"{OUT}/matched_plof_oe_comparison.png", dpi=200, bbox_inches="tight")
fig4.savefig(f"{OUT}/matched_plof_oe_comparison.pdf", bbox_inches="tight")
print(f"Saved comparison plot to {OUT}/matched_plof_oe_comparison.png")
