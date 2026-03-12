"""Visualize Proemis3D constraint regions per gene across forward methods.

Standalone script that reads the union TSV produced by
``compare_forward_methods_report.py --export-tsv`` and generates per-gene
comparison plots in three styles:

1. **tracks** -- genome-browser-style linear tracks (one row per method).
2. **bars** -- compact segmented bars (thinner, tighter spacing).
3. **heatmap** -- residue-level heatmap (x = residue, y = method).

All styles colour constraint regions by OE upper CI using the gnomAD
regional missense constraint palette (dark red = low OE / most
constrained, through orange to cream = high OE / unconstrained).
Null catch-all regions are shown in light grey, and NA / unassigned
residues in white with hatching.

Usage::

    python visualize_gene_methods.py \\
        --tsv forward_combined.tsv \\
        --gene P12345 \\
        [--style all|tracks|bars|heatmap] \\
        [--output-dir .] \\
        [--vmin 0] [--vmax 1.5]
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Font scaling
# ---------------------------------------------------------------------------

_BASE_FIG_WIDTH = 14.0  # Reference width where hardcoded sizes look good.


def _fs(base_size: float, fig_width: float) -> float:
    """Scale a font size proportionally to figure width."""
    return base_size * (fig_width / _BASE_FIG_WIDTH)


def _height_scale(fig_width: float) -> float:
    """Scaling factor for row/figure heights to accommodate larger fonts.

    Uses a gentler square-root scaling so height grows more slowly than
    font size but still keeps text from looking squished.
    """
    ratio = fig_width / _BASE_FIG_WIDTH
    if ratio <= 1.0:
        return 1.0
    return 1.0 + (ratio - 1.0) * 0.5


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _read_and_normalize_tsv(tsv_path: str) -> pd.DataFrame:
    """Read the union TSV and add ``_is_null`` boolean column."""
    print(f"Reading {tsv_path} ...")
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if df.columns.duplicated().any():
        seen: Dict[str, int] = {}
        new_cols: List[str] = []
        for c in df.columns:
            if c in seen:
                seen[c] += 1
                new_cols.append(f"{c}.{seen[c]}")
            else:
                seen[c] = 0
                new_cols.append(c)
        df.columns = new_cols
    if "run_label" not in df.columns:
        raise ValueError("TSV must contain a 'run_label' column.")
    if "is_null" in df.columns:
        df["_is_null"] = (
            df["is_null"].fillna(False).astype(str).str.lower().isin(("true", "1"))
        )
    else:
        df["_is_null"] = False
    return df


# ---------------------------------------------------------------------------
# Shared data preparation
# ---------------------------------------------------------------------------

def prepare_gene_data(
    df: pd.DataFrame,
    uniprot_id: str,
) -> Tuple[Dict[str, pd.DataFrame], int, List[str]]:
    """Filter *df* to a single gene and return per-method DataFrames.

    For each ``run_label`` the returned DataFrame has columns:

    * ``residue_index`` -- 0-based residue position
    * ``region_index`` -- constraint region id (NaN for unassigned)
    * ``oe_upper`` -- OE upper CI for the region (NaN for unassigned / null)
    * ``is_null`` -- bool, True for the null catch-all region

    Returns
    -------
    per_method : dict[str, DataFrame]
        Keyed by *run_label*.
    protein_length : int
        Number of residues (max residue_index + 1).
    run_labels : list[str]
        Ordered list of run labels.
    """
    gene_df = df[df["uniprot_id"] == uniprot_id].copy()
    if gene_df.empty:
        raise ValueError(f"UniProt ID '{uniprot_id}' not found in TSV.")

    # Pick the most common transcript if there are multiple.
    tx_counts = gene_df.groupby("transcript_id")["residue_index"].count()
    best_tx = tx_counts.idxmax()
    gene_df = gene_df[gene_df["transcript_id"] == best_tx]

    protein_length = int(gene_df["residue_index"].max()) + 1
    run_labels = list(gene_df["run_label"].unique())

    per_method: Dict[str, pd.DataFrame] = {}
    for rl in run_labels:
        sub = gene_df[gene_df["run_label"] == rl][
            ["residue_index", "region_index", "oe_upper", "_is_null"]
        ].copy()
        sub = sub.rename(columns={"_is_null": "is_null"})
        sub = sub.sort_values("residue_index").reset_index(drop=True)
        per_method[rl] = sub

    return per_method, protein_length, run_labels


# ---------------------------------------------------------------------------
# Method ordering
# ---------------------------------------------------------------------------

# Canonical sort keys for the three axes of variation in the run labels.
_REGION_METHOD_ORDER = {"AIC": 0, "AIC weight": 1, "LRT": 2}
_PLDDT_ORDER = {"none": 0, "pLDDT remove": 1, "pLDDT exclude": 2}
_PAE_ORDER = {
    "none": 0,
    "PAE truncate_center": 1,
    "PAE filter_center": 2,
    "PAE filter_region": 3,
    "PAE exclude_center": 4,
    "PAE exclude_region": 5,
}


def _parse_method_components(label: str) -> Tuple[str, str, str]:
    """Extract (region_method, plddt, pae) from a run label string."""
    parts = [p.strip() for p in label.split(",")]

    # Region selection method is always the first token.
    region_method = parts[0]
    # Strip the "(no pLDDT/PAE)" suffix if present.
    if "(" in region_method:
        region_method = region_method[: region_method.index("(")].strip()

    plddt = "none"
    pae = "none"
    for p in parts[1:]:
        if p.startswith("pLDDT"):
            plddt = p
        elif p.startswith("PAE"):
            pae = p

    return region_method, plddt, pae


def sort_by_method(run_labels: List[str]) -> List[str]:
    """Sort run labels hierarchically: region method > pLDDT > PAE.

    Methods with the same region-selection approach are grouped together,
    then sub-sorted by pLDDT filter and PAE filter.
    """

    def _sort_key(label: str) -> Tuple[int, int, int, str]:
        rm, pl, pa = _parse_method_components(label)
        return (
            _REGION_METHOD_ORDER.get(rm, 99),
            _PLDDT_ORDER.get(pl, 99),
            _PAE_ORDER.get(pa, 99),
            label,
        )

    return sorted(run_labels, key=_sort_key)


def sort_by_filter_group(run_labels: List[str]) -> List[str]:
    """Sort run labels by filter group first, then by model (AIC before LRT).

    Order: no filter (AIC, LRT) → pLDDT (AIC, LRT) → PAE (AIC, LRT).
    """

    def _sort_key(label: str) -> Tuple[int, int, str]:
        rm, pl, pa = _parse_method_components(label)
        if pl == "none" and pa == "none":
            group = 0
        elif pl != "none":
            group = 1
        else:
            group = 2
        return (group, _REGION_METHOD_ORDER.get(rm, 99), label)

    return sorted(run_labels, key=_sort_key)


def sort_by_cluster(
    per_method: Dict[str, pd.DataFrame],
    protein_length: int,
    run_labels: List[str],
) -> List[str]:
    """Reorder run labels by hierarchical clustering of residue-level patterns.

    Each method is represented as a vector of ``oe_upper`` values per residue
    (NaN for missing / null).  Methods with similar spatial constraint
    patterns end up adjacent.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist, squareform

    n = len(run_labels)
    if n <= 2:
        return list(run_labels)

    # Build residue-level matrix (methods x residues).
    matrix = np.full((n, protein_length), np.nan)
    for idx, rl in enumerate(run_labels):
        sub = per_method[rl]
        for _, row in sub.iterrows():
            ri = int(row["residue_index"])
            if ri < protein_length and not row["is_null"] and pd.notna(row["oe_upper"]):
                matrix[idx, ri] = float(row["oe_upper"])

    # Pairwise correlation distance, treating NaN as 0 after centering.
    # Replace NaN with per-method mean so missing residues don't dominate.
    row_means = np.nanmean(matrix, axis=1, keepdims=True)
    row_means = np.where(np.isnan(row_means), 0, row_means)
    filled = np.where(np.isnan(matrix), row_means, matrix)

    dist = pdist(filled, metric="correlation")
    # Replace any NaN distances (constant rows) with max distance.
    dist = np.where(np.isnan(dist), 1.0, dist)
    Z = linkage(dist, method="average")
    order = leaves_list(Z)

    return [run_labels[i] for i in order]


def _get_region_spans(
    sub: pd.DataFrame,
) -> List[Tuple[int, int, float, str]]:
    """Return contiguous region spans as (start, length, oe_upper, kind).

    *kind* is one of ``"constraint"``, ``"null"``, or ``"na"``.
    """
    spans: List[Tuple[int, int, float, str]] = []
    if sub.empty:
        return spans

    prev_key = None
    span_start = 0
    span_len = 0
    span_oe: float = np.nan
    span_kind = "na"

    for _, row in sub.iterrows():
        ri = int(row["residue_index"])
        oe_val = row["oe_upper"]
        is_null = bool(row["is_null"])
        reg = row["region_index"]

        if is_null:
            kind = "null"
            key = ("null",)
        elif pd.isna(reg) or pd.isna(oe_val):
            kind = "na"
            key = ("na",)
        else:
            kind = "constraint"
            key = ("constraint", int(reg))

        if key != prev_key or ri != span_start + span_len:
            if prev_key is not None:
                spans.append((span_start, span_len, span_oe, span_kind))
            span_start = ri
            span_len = 1
            span_oe = oe_val if kind != "na" else np.nan
            span_kind = kind
            prev_key = key
        else:
            span_len += 1

    if prev_key is not None:
        spans.append((span_start, span_len, span_oe, span_kind))

    return spans


# ---------------------------------------------------------------------------
# Colour helpers — gnomAD regional missense constraint palette
# ---------------------------------------------------------------------------

# Discrete gnomAD RMC bin colours (one flat colour per bin).
_RMC_BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
_RMC_COLORS = [
    "#67000d",   # [0.0, 0.2)
    "#cb181d",   # [0.2, 0.4)
    "#ef6548",   # [0.4, 0.6)
    "#fc8d59",   # [0.6, 0.8)
    "#fdbb84",   # [0.8, 1.0)
    "#fef0d9",   # [1.0, +inf)
]
_NULL_COLOR = "#d9d9d9"
_NA_COLOR = "#f7f7f7"
_EDGE_COLOR = "#333333"

# Discrete colourmap + boundary norm used by all three plot styles.
_RMC_CMAP = mcolors.ListedColormap(_RMC_COLORS)
_RMC_BOUNDARIES = _RMC_BINS + [999.0]  # upper sentinel so last bin catches 1.0+
_RMC_NORM = mcolors.BoundaryNorm(_RMC_BOUNDARIES, _RMC_CMAP.N)


def _oe_color(val: float, cmap, norm) -> Any:
    """Map an OE upper CI value to a discrete bin colour."""
    if np.isnan(val):
        return _NULL_COLOR
    return cmap(norm(val))


def _add_discrete_colorbar(fig, ax, fig_width: float = _BASE_FIG_WIDTH) -> None:
    """Add a discrete colorbar with 6 equal-sized bins matching the RMC palette."""
    from matplotlib.patches import Patch

    # Use evenly-spaced boundaries for display only (no 999 sentinel).
    display_bounds = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    display_norm = mcolors.BoundaryNorm(display_bounds, _RMC_CMAP.N)
    sm = plt.cm.ScalarMappable(cmap=_RMC_CMAP, norm=display_norm)
    sm.set_array([])
    cbar = fig.colorbar(
        sm, ax=ax, fraction=0.02, pad=0.02, spacing="uniform",
    )
    # Tick at each boundary, label the last one as "1.0+"
    cbar.set_ticks(display_bounds)
    cbar.ax.set_yticklabels(
        ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", ""],
        fontsize=_fs(10, fig_width),
    )
    cbar.set_label("OE upper CI", fontsize=_fs(12, fig_width))


def _count_constraint_regions(sub: pd.DataFrame) -> int:
    """Count distinct constraint regions (excludes null and NA).

    Uses the actual ``region_index`` values from the data rather than
    visual spans, so a single region split by NA gaps is counted once.
    """
    mask = (~sub["is_null"]) & sub["region_index"].notna()
    return int(sub.loc[mask, "region_index"].nunique())


def _format_labels(run_labels: List[str]) -> List[str]:
    """Replace commas in run labels with newlines for multi-line y-tick labels."""
    return [rl.replace(", ", "\n") for rl in run_labels]


def _add_null_na_legend(ax, fig_width: float = _BASE_FIG_WIDTH, color_null: bool = False) -> None:
    """Add legend patches for null and NA regions to *ax*.

    When *color_null* is True the catch-all region is coloured by OE so only
    the NA / unassigned patch is added (the catch-all blends into the main
    colour scale).
    """
    from matplotlib.patches import Patch

    legend_elements = []
    if not color_null:
        legend_elements.append(
            Patch(facecolor=_NULL_COLOR, edgecolor=_EDGE_COLOR, linewidth=0.5,
                  label="Catch-all region")
        )
    legend_elements.append(
        Patch(facecolor=_NA_COLOR, edgecolor=_EDGE_COLOR, linewidth=0.5,
              hatch="///", label="NA / unassigned")
    )
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=_fs(10, fig_width),
        framealpha=0.9,
        edgecolor="#cccccc",
    )


# ---------------------------------------------------------------------------
# Plot 1: Linear tracks (genome-browser style)
# ---------------------------------------------------------------------------

def plot_tracks(
    per_method: Dict[str, pd.DataFrame],
    protein_length: int,
    run_labels: List[str],
    title: str,
    path: str,
    vmin: float = 0.0,
    vmax: float = 1.5,
    dpi: int = 150,
) -> None:
    """Genome-browser-style linear tracks, one row per method."""
    cmap = _RMC_CMAP
    norm = _RMC_NORM

    n_methods = len(run_labels)
    display_labels = _format_labels(run_labels)
    fig_width = max(14, protein_length / 70)
    hs = _height_scale(fig_width)
    row_height = 0.7 * hs
    row_gap = 0.3 * hs
    total_height = n_methods * (row_height + row_gap) + 2.0 * hs

    fig, ax = plt.subplots(figsize=(fig_width, total_height), dpi=dpi)

    region_counts: List[int] = []
    for idx, rl in enumerate(run_labels):
        y_center = (n_methods - 1 - idx) * (row_height + row_gap)
        sub = per_method[rl]
        spans = _get_region_spans(sub)
        region_counts.append(_count_constraint_regions(sub))

        for start, length, oe_val, kind in spans:
            if kind == "constraint":
                color = _oe_color(oe_val, cmap, norm)
                ax.broken_barh(
                    [(start, length)],
                    (y_center, row_height),
                    facecolors=color,
                    edgecolors=_EDGE_COLOR,
                    linewidth=0.4,
                )
            elif kind == "null":
                ax.broken_barh(
                    [(start, length)],
                    (y_center, row_height),
                    facecolors=_NULL_COLOR,
                    edgecolors=_EDGE_COLOR,
                    linewidth=0.3,
                )
            else:
                ax.broken_barh(
                    [(start, length)],
                    (y_center, row_height),
                    facecolors=_NA_COLOR,
                    edgecolors="#aaaaaa",
                    linewidth=0,
                    hatch="///",
                )

    ax.set_xlim(-0.5, protein_length + 0.5)
    ax.set_ylim(-row_gap, n_methods * (row_height + row_gap))
    ax.set_yticks(
        [(n_methods - 1 - i) * (row_height + row_gap) + row_height / 2 for i in range(n_methods)]
    )
    ax.set_yticklabels(display_labels, fontsize=_fs(9, fig_width))
    ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
    ax.set_title(title, fontsize=_fs(15, fig_width), pad=12)
    ax.tick_params(axis="x", labelsize=_fs(10, fig_width))

    # Annotate region counts on the right side of each track.
    for idx, (rl, n_reg) in enumerate(zip(run_labels, region_counts)):
        y_text = (n_methods - 1 - idx) * (row_height + row_gap) + row_height / 2
        ax.text(
            protein_length + protein_length * 0.01,
            y_text,
            str(n_reg),
            va="center",
            ha="left",
            fontsize=_fs(11, fig_width),
            color="#555555",
        )
    # Label for the region count column.
    ax.text(
        protein_length + protein_length * 0.01,
        n_methods * (row_height + row_gap) - row_gap * 0.3,
        "# regions",
        va="bottom",
        ha="left",
        fontsize=_fs(11, fig_width),
        fontweight="bold",
        color="#555555",
    )

    _add_discrete_colorbar(fig, ax, fig_width)
    _add_null_na_legend(ax, fig_width)

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved tracks plot: {path}")


# ---------------------------------------------------------------------------
# Plot 2: Segmented bars (compact)
# ---------------------------------------------------------------------------

def plot_bars(
    per_method: Dict[str, pd.DataFrame],
    protein_length: int,
    run_labels: List[str],
    title: str,
    path: str,
    vmin: float = 0.0,
    vmax: float = 1.5,
    dpi: int = 150,
    color_null: bool = False,
) -> None:
    """Compact segmented bars, one thin bar per method.

    When *color_null* is True, catch-all (null) regions are coloured by their
    OE upper CI value using the same palette as constraint regions rather than
    the default flat grey.
    """
    cmap = _RMC_CMAP
    norm = _RMC_NORM

    n_methods = len(run_labels)
    display_labels = _format_labels(run_labels)
    fig_width = max(14, protein_length / 70)
    hs = _height_scale(fig_width)
    bar_height = 0.4 * hs
    bar_gap = 0.15 * hs
    total_height = n_methods * (bar_height + bar_gap) + 2.0 * hs

    fig, ax = plt.subplots(figsize=(fig_width, total_height), dpi=dpi)

    region_counts: List[int] = []
    for idx, rl in enumerate(run_labels):
        y_center = (n_methods - 1 - idx) * (bar_height + bar_gap)
        sub = per_method[rl]
        spans = _get_region_spans(sub)
        region_counts.append(_count_constraint_regions(sub))

        for start, length, oe_val, kind in spans:
            if kind == "constraint":
                color = _oe_color(oe_val, cmap, norm)
                ax.barh(
                    y_center + bar_height / 2,
                    length,
                    left=start,
                    height=bar_height,
                    color=color,
                    edgecolor=_EDGE_COLOR,
                    linewidth=0.3,
                )
            elif kind == "null":
                null_color = _oe_color(oe_val, cmap, norm) if color_null else _NULL_COLOR
                ax.barh(
                    y_center + bar_height / 2,
                    length,
                    left=start,
                    height=bar_height,
                    color=null_color,
                    edgecolor=_EDGE_COLOR,
                    linewidth=0.2,
                )
            else:
                ax.barh(
                    y_center + bar_height / 2,
                    length,
                    left=start,
                    height=bar_height,
                    color=_NA_COLOR,
                    edgecolor="#aaaaaa",
                    linewidth=0,
                    hatch="///",
                )

    ax.set_xlim(-0.5, protein_length + 0.5)
    ax.set_ylim(-bar_gap, n_methods * (bar_height + bar_gap))
    ax.set_yticks(
        [(n_methods - 1 - i) * (bar_height + bar_gap) + bar_height / 2 for i in range(n_methods)]
    )
    ax.set_yticklabels(display_labels, fontsize=_fs(8, fig_width))
    ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
    ax.set_title(title, fontsize=_fs(15, fig_width), pad=12)
    ax.tick_params(axis="x", labelsize=_fs(10, fig_width))

    # Annotate region counts on the right side of each bar.
    for idx, (rl, n_reg) in enumerate(zip(run_labels, region_counts)):
        y_text = (n_methods - 1 - idx) * (bar_height + bar_gap) + bar_height / 2
        ax.text(
            protein_length + protein_length * 0.01,
            y_text,
            str(n_reg),
            va="center",
            ha="left",
            fontsize=_fs(11, fig_width),
            color="#555555",
        )
    ax.text(
        protein_length + protein_length * 0.01,
        n_methods * (bar_height + bar_gap) - bar_gap * 0.3,
        "# regions",
        va="bottom",
        ha="left",
        fontsize=_fs(11, fig_width),
        fontweight="bold",
        color="#555555",
    )

    _add_discrete_colorbar(fig, ax, fig_width)
    _add_null_na_legend(ax, fig_width, color_null=color_null)

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved bars plot: {path}")


# ---------------------------------------------------------------------------
# Plot 3: Residue-level heatmap
# ---------------------------------------------------------------------------

_NULL_SENTINEL = -1.0
_NA_SENTINEL = -2.0


def plot_heatmap(
    per_method: Dict[str, pd.DataFrame],
    protein_length: int,
    run_labels: List[str],
    title: str,
    path: str,
    vmin: float = 0.0,
    vmax: float = 1.5,
    dpi: int = 150,
) -> None:
    """Residue-level heatmap: x = residue, y = method, colour = oe_upper."""
    n_methods = len(run_labels)
    display_labels = _format_labels(run_labels)
    matrix = np.full((n_methods, protein_length), _NA_SENTINEL, dtype=float)

    for idx, rl in enumerate(run_labels):
        sub = per_method[rl]
        for _, row in sub.iterrows():
            ri = int(row["residue_index"])
            if ri >= protein_length:
                continue
            if row["is_null"]:
                matrix[idx, ri] = _NULL_SENTINEL
            elif pd.notna(row["oe_upper"]):
                matrix[idx, ri] = float(row["oe_upper"])
            # else stays _NA_SENTINEL

    # Build a discrete colourmap: NA sentinel, null sentinel, then the 6
    # gnomAD RMC bin colours.
    heatmap_colors = [_NA_COLOR, _NULL_COLOR] + _RMC_COLORS
    heatmap_cmap = mcolors.ListedColormap(heatmap_colors)
    # Boundaries: NA < -1.5 < null < -0.5 < bin edges ...
    heatmap_bounds = [_NA_SENTINEL - 0.5, -1.5, -0.5] + _RMC_BINS[1:] + [999.0]
    heatmap_norm = mcolors.BoundaryNorm(heatmap_bounds, heatmap_cmap.N)

    fig_width = max(14, protein_length / 45)
    hs = _height_scale(fig_width)
    fig_height = max(5 * hs, n_methods * 0.4 * hs + 2.0 * hs)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

    ax.imshow(
        matrix,
        aspect="auto",
        cmap=heatmap_cmap,
        norm=heatmap_norm,
        interpolation="nearest",
        origin="upper",
    )

    # Overlay hatching on NA cells.  For each row, find contiguous NA spans
    # and draw a single hatched rectangle per span for efficiency.
    from matplotlib.patches import Rectangle

    for row_idx in range(n_methods):
        is_na = matrix[row_idx] == _NA_SENTINEL
        if not is_na.any():
            continue
        # Find contiguous True runs.
        diff = np.diff(is_na.astype(int))
        starts = np.where(diff == 1)[0] + 1
        ends = np.where(diff == -1)[0] + 1
        if is_na[0]:
            starts = np.concatenate(([0], starts))
        if is_na[-1]:
            ends = np.concatenate((ends, [protein_length]))
        for s, e in zip(starts, ends):
            ax.add_patch(
                Rectangle(
                    (s - 0.5, row_idx - 0.5),
                    e - s,
                    1,
                    fill=False,
                    hatch="///",
                    edgecolor="#aaaaaa",
                    linewidth=0,
                )
            )

    ax.set_yticks(range(n_methods))
    ax.set_yticklabels(display_labels, fontsize=_fs(9, fig_width))
    ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
    ax.set_title(title, fontsize=_fs(15, fig_width), pad=12)
    ax.tick_params(axis="x", labelsize=_fs(10, fig_width))

    # Colourbar showing only the data bins (not the sentinels).
    _add_discrete_colorbar(fig, ax, fig_width)
    _add_null_na_legend(ax, fig_width)

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved heatmap plot: {path}")


# ---------------------------------------------------------------------------
# Plot 4: AlphaFold structure metrics (pLDDT + PAE)
# ---------------------------------------------------------------------------

_PLDDT_COLORS = [
    (0, "#FF7D45"),    # Very low (0-50): orange
    (50, "#FFDB13"),   # Low (50-70): yellow
    (70, "#65CBF3"),   # Confident (70-90): light blue
    (90, "#0053D6"),   # Very high (90-100): dark blue
]


def _load_plddt(path: str) -> List[float]:
    """Load pLDDT scores from an AlphaFold confidence JSON file."""
    import json

    with open(path) as f:
        data = json.load(f)
    # AlphaFold v4 format: {"confidenceScore": [...]}
    # Older format may have different keys.
    if "confidenceScore" in data:
        return data["confidenceScore"]
    if isinstance(data, list) and "plddt" in data[0]:
        return data[0]["plddt"]
    raise ValueError(f"Cannot parse pLDDT from {path}. Expected 'confidenceScore' key.")


def _load_pae(path: str) -> np.ndarray:
    """Load PAE matrix from an AlphaFold PAE JSON file."""
    import json

    with open(path) as f:
        data = json.load(f)
    # AlphaFold v4 format: [{"predicted_aligned_error": [[...]]}]
    if isinstance(data, list) and "predicted_aligned_error" in data[0]:
        return np.array(data[0]["predicted_aligned_error"], dtype=float)
    # Newer format: {"predicted_aligned_error": [[...]]}
    if "predicted_aligned_error" in data:
        return np.array(data["predicted_aligned_error"], dtype=float)
    raise ValueError(f"Cannot parse PAE from {path}. Expected 'predicted_aligned_error' key.")


def _dist_matrix_from_cif(cif_path: str) -> np.ndarray:
    """Compute Cα distance matrix from an mmCIF file."""
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("af2", cif_path)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True) and "CA" in residue:
                    ca_coords.append(residue["CA"].get_vector().get_array())
    coords = np.array(ca_coords)
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))


def _load_af2_from_gcs(
    gcs_dir: str, uniprot_id: str,
) -> Tuple[Optional[List[float]], Optional[np.ndarray], Optional[np.ndarray]]:
    """Load pLDDT, PAE and distance matrix from gzipped AlphaFold files in GCS.

    Returns ``(plddt, pae, dist_matrix)`` -- any may be ``None`` on failure.
    """
    import gzip
    import json
    import subprocess
    import tempfile

    gcs_dir = gcs_dir.rstrip("/")
    plddt_data: Optional[List[float]] = None
    pae_data: Optional[np.ndarray] = None
    dist_data: Optional[np.ndarray] = None

    for kind, parser in [("confidence", "plddt"), ("predicted_aligned_error", "pae")]:
        gcs_path = f"{gcs_dir}/AF-{uniprot_id}-F1-{kind}_v4.json.gz"
        tmp_path = os.path.join(tempfile.gettempdir(), f"AF-{uniprot_id}-{kind}.json.gz")
        print(f"Downloading {gcs_path} ...")
        result = subprocess.run(
            ["gsutil", "-q", "cp", gcs_path, tmp_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"  Could not download {gcs_path}: {result.stderr.strip()}")
            continue
        with gzip.open(tmp_path, "rt") as f:
            data = json.load(f)
        if parser == "plddt":
            if "confidenceScore" in data:
                plddt_data = data["confidenceScore"]
            elif isinstance(data, list) and len(data) > 0 and "plddt" in data[0]:
                plddt_data = data[0]["plddt"]
            else:
                print(f"  Could not parse pLDDT from {gcs_path}")
        else:
            if isinstance(data, list) and "predicted_aligned_error" in data[0]:
                pae_data = np.array(data[0]["predicted_aligned_error"], dtype=float)
            elif "predicted_aligned_error" in data:
                pae_data = np.array(data["predicted_aligned_error"], dtype=float)
            else:
                print(f"  Could not parse PAE from {gcs_path}")

    # Distance matrix from CIF.
    cif_gcs = f"{gcs_dir}/AF-{uniprot_id}-F1-model_v4.cif.gz"
    cif_tmp = os.path.join(tempfile.gettempdir(), f"AF-{uniprot_id}-model.cif.gz")
    print(f"Downloading {cif_gcs} ...")
    result = subprocess.run(
        ["gsutil", "-q", "cp", cif_gcs, cif_tmp],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        import gzip as _gzip
        cif_plain = cif_tmp.replace(".gz", "")
        with _gzip.open(cif_tmp, "rb") as fin, open(cif_plain, "wb") as fout:
            fout.write(fin.read())
        try:
            dist_data = _dist_matrix_from_cif(cif_plain)
        except Exception as exc:
            print(f"  Could not compute distance matrix: {exc}")
    else:
        print(f"  Could not download {cif_gcs}: {result.stderr.strip()}")

    return plddt_data, pae_data, dist_data


def _fetch_alphafold_json(uniprot_id: str, kind: str) -> Optional[str]:
    """Fetch AlphaFold DB JSON via the prediction API.

    Uses ``https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`` to
    discover the correct URLs (version may vary).

    Returns file content as a string, or None on failure.
    """
    import json
    import urllib.request
    import urllib.error

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        print(f"Querying AlphaFold API for {uniprot_id} ...")
        with urllib.request.urlopen(api_url, timeout=15) as resp:
            entries = json.loads(resp.read())
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        print(f"  AlphaFold API query failed: {exc}")
        return None

    if not entries:
        print(f"  No AlphaFold entry for {uniprot_id}.")
        return None

    # Use the first (canonical) entry.
    entry = entries[0]
    if kind == "pae":
        url = entry.get("paeDocUrl")
    elif kind == "plddt":
        url = entry.get("plddtDocUrl")
    else:
        return None

    if not url:
        print(f"  No {kind} URL in AlphaFold entry for {uniprot_id}.")
        return None

    try:
        print(f"Fetching {url} ...")
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read().decode()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        print(f"  Could not fetch {kind}: {exc}")
        return None


# ---------------------------------------------------------------------------
# UniProt feature loading & plotting
# ---------------------------------------------------------------------------

# Feature types to display and their colours (order = top to bottom).
_FEATURE_DISPLAY_ORDER = [
    "domain",
    "topological domain",
    "transmembrane region",
    "intramembrane region",
    "region of interest",
    "short sequence motif",
    "coiled-coil region",
    "repeat",
    "signal peptide",
    "propeptide",
    "disulfide bond",
    "zinc finger region",
    "DNA-binding region",
    "active site",
    "metal ion-binding site",
    "binding site",
    "glycosylation site",
]

# Merge old-style UPPERCASE keys into their lowercase equivalents.
_FEATURE_KEY_ALIASES: Dict[str, str] = {
    "DOMAIN": "domain",
    "CHAIN": "chain",
    "chain": "chain",
    "REGION": "region of interest",
    "TRANSMEM": "transmembrane region",
    "TOPO_DOM": "topological domain",
    "COMPBIAS": "compositionally biased region",
    "COILED": "coiled-coil region",
    "SIGNAL": "signal peptide",
    "DISULFID": "disulfide bond",
    "ZN_FING": "zinc finger region",
    "DNA_BIND": "DNA-binding region",
    "dna-binding region": "DNA-binding region",
    "REPEAT": "repeat",
    "METAL": "metal ion-binding site",
    "BINDING": "binding site",
    "CARBOHYD": "glycosylation site",
    "ACT_SITE": "active site",
    "NP_BIND": "binding site",
    "CA_BIND": "metal ion-binding site",
    "MOTIF": "region of interest",
    "SITE": "active site",
}

_FEATURE_COLORS: Dict[str, str] = {
    "domain": "#2166ac",
    "topological domain": "#4393c3",
    "transmembrane region": "#d6604d",
    "intramembrane region": "#f4a442",
    "region of interest": "#92c5de",
    "short sequence motif": "#a6611a",
    "coiled-coil region": "#fddbc7",
    "compositionally biased region": "#bababa",
    "repeat": "#c2a5cf",
    "signal peptide": "#1b7837",
    "propeptide": "#7fbf7b",
    "disulfide bond": "#e7298a",
    "zinc finger region": "#66a61e",
    "DNA-binding region": "#e6ab02",
    "active site": "#b2182b",
    "metal ion-binding site": "#f4a582",
    "binding site": "#d95f02",
    "glycosylation site": "#7570b3",
}


def _load_features_json(
    path: str, uniprot_id: str,
) -> Optional[Dict[str, List[dict]]]:
    """Load UniProt features for *uniprot_id* from a local JSON file.

    The JSON is expected to be ``{uniprot_id: {feature_type: [{start, end, note, ...}]}}``.
    """
    import json

    with open(path) as fh:
        data = json.load(fh)
    gene_data = data.get(uniprot_id)
    if gene_data is None:
        print(f"  {uniprot_id} not found in features JSON {path}")
        return None
    return gene_data


def _export_features_from_ht(
    ht_path: str, uniprot_ids: List[str], out_path: str,
) -> str:
    """Export UniProt features from a Hail Table to a local JSON file.

    Requires the Hail *batch* backend.  Returns *out_path*.
    """
    import json

    import hail as hl

    hl.init(backend="batch", tmp_dir="gs://gnomad-tmp-4day", idempotent=True)
    ht = hl.read_table(ht_path)

    uid_set = hl.literal(set(uniprot_ids))
    ht = ht.filter(uid_set.contains(ht.uniprot_id))
    ht_exp = ht.explode("uniprot_2021_04_features")

    rows = ht_exp.select(
        ft=ht_exp.uniprot_2021_04_features.feature_type,
        note=ht_exp.uniprot_2021_04_features.note,
        start=ht_exp.uniprot_2021_04_features.interval.start,
        end=ht_exp.uniprot_2021_04_features.interval.end,
        fid=ht_exp.uniprot_2021_04_features.feature_id,
    ).collect()

    result: Dict[str, Dict[str, list]] = {}
    seen: Dict[str, set] = {}
    for r in rows:
        uid = r.uniprot_id
        ft = r.ft
        key = (uid, ft, r.start, r.end, r.note)
        if key in seen.get(uid, set()):
            continue
        seen.setdefault(uid, set()).add(key)
        result.setdefault(uid, {}).setdefault(ft, []).append(
            {"start": r.start, "end": r.end, "note": r.note, "feature_id": r.fid}
        )

    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"Exported features for {len(result)} gene(s) to {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Variant data loading
# ---------------------------------------------------------------------------

# Variant track definitions: (display_name, colour, marker).
_VARIANT_TRACKS = [
    ("Autism (case)", "#e41a1c", "^"),
    ("Autism (control)", "#377eb8", "v"),
    ("DD de novo", "#ff7f00", "s"),
    ("gnomAD v4 de novo", "#984ea3", "o"),
]


def _load_variants(
    tsv_path: str, uniprot_id: str,
) -> Dict[str, List[int]]:
    """Load variant residue positions from the autism/DD TSV for a gene.

    Returns ``{track_name: [residue_index, ...]}`` for every track that
    has at least one variant.
    """
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    sub = df[df["uniprot_id"] == uniprot_id].copy()
    if sub.empty:
        return {}

    result: Dict[str, List[int]] = {}

    # Autism case (ASD).
    mask_asd = sub["variant_level_annotations.autism.role"].str.contains(
        "ASD", na=False,
    )
    pos = sub.loc[mask_asd, "residue_index"].dropna().astype(int).tolist()
    if pos:
        result["Autism (case)"] = sorted(set(pos))

    # Autism control.
    mask_ctrl = sub["variant_level_annotations.autism.role"].str.contains(
        "control", na=False,
    )
    pos = sub.loc[mask_ctrl, "residue_index"].dropna().astype(int).tolist()
    if pos:
        result["Autism (control)"] = sorted(set(pos))

    # DD de novo.
    mask_dd = sub[
        "variant_level_annotations.dd_denovo_no_transcript_match.case_control"
    ].str.contains("DD", na=False)
    pos = sub.loc[mask_dd, "residue_index"].dropna().astype(int).tolist()
    if pos:
        result["DD de novo"] = sorted(set(pos))

    # gnomAD v4 de novo.
    mask_gdn = sub[
        "variant_level_annotations.gnomad_de_novo.de_novo_AC.AC_all"
    ].notna()
    pos = sub.loc[mask_gdn, "residue_index"].dropna().astype(int).tolist()
    if pos:
        result["gnomAD v4 de novo"] = sorted(set(pos))

    return result


# ---------------------------------------------------------------------------
# Plot 5: UniProt features + variant tracks
# ---------------------------------------------------------------------------

def plot_features(
    features: Dict[str, List[dict]],
    protein_length: int,
    title: str,
    path: str,
    variants: Optional[Dict[str, List[int]]] = None,
    dpi: int = 150,
) -> None:
    """Plot UniProt protein features and variant tracks.

    One row per feature type / variant source.  UniProt spans are coloured
    bars; point features and variants are drawn as markers.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    if variants is None:
        variants = {}

    # Normalise feature keys: merge old-style UPPERCASE names into canonical
    # lowercase equivalents so both formats are displayed.
    normalised: Dict[str, List[dict]] = {}
    for key, entries in features.items():
        canonical = _FEATURE_KEY_ALIASES.get(key, key)
        normalised.setdefault(canonical, []).extend(entries)
    features = normalised

    # Feature types that have data.
    display_types = [
        ft for ft in _FEATURE_DISPLAY_ORDER if ft in features and features[ft]
    ]

    # Variant tracks that have data (keep canonical order).
    variant_track_names = [
        name for name, _, _ in _VARIANT_TRACKS if name in variants
    ]
    variant_meta = {name: (col, mkr) for name, col, mkr in _VARIANT_TRACKS}

    n_feature_rows = len(display_types)
    n_variant_rows = len(variant_track_names)
    n_rows = n_feature_rows + n_variant_rows

    if n_rows == 0:
        print("No displayable features or variants; skipping feature plot.")
        return

    fig_width = max(14, protein_length / 45)
    hs = _height_scale(fig_width)
    row_height = 0.7 * hs
    row_gap = 0.3 * hs
    row_step = row_height + row_gap
    fig_height = n_rows * row_step + 2.5 * hs

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

    # All row labels and y-centers (bottom to top): variant tracks first, then features.
    all_labels: List[str] = []
    y_centers: List[float] = []

    # --- Variant tracks (bottom rows) ---
    for row_idx, vt_name in enumerate(reversed(variant_track_names)):
        y_center = row_idx * row_step
        colour, marker = variant_meta[vt_name]
        positions = variants[vt_name]
        for pos in positions:
            ax.plot(
                pos, y_center, marker=marker,
                markersize=_fs(6, fig_width),
                color=colour, markeredgecolor="black", markeredgewidth=0.4,
                zorder=3, linestyle="none",
            )
        all_labels.append(vt_name)
        y_centers.append(y_center)

    # --- UniProt feature tracks (top rows) ---
    for feat_idx, ft in enumerate(reversed(display_types)):
        abs_idx = n_variant_rows + feat_idx
        y_center = abs_idx * row_step
        colour = _FEATURE_COLORS.get(ft, "#999999")
        entries = features[ft]

        for entry in entries:
            s, e = entry["start"], entry["end"]
            if s == e:
                ax.plot(
                    s, y_center, marker="D", markersize=_fs(6, fig_width),
                    color=colour, markeredgecolor="black", markeredgewidth=0.4,
                    zorder=3,
                )
            else:
                ax.broken_barh(
                    [(s, e - s + 1)],
                    (y_center - row_height / 2, row_height),
                    facecolors=colour,
                    edgecolors="black",
                    linewidth=0.3,
                    zorder=2,
                )
                span_width = e - s + 1
                if span_width > protein_length * 0.04:
                    label = entry.get("note", "")
                    if len(label) > 25:
                        label = label[:22] + "..."
                    ax.text(
                        s + span_width / 2, y_center, label,
                        ha="center", va="center", fontsize=_fs(7, fig_width),
                        color="white", fontweight="bold", clip_on=True,
                        zorder=4,
                    )
        all_labels.append(ft.title())
        y_centers.append(y_center)

    # Draw a thin separator line between features and variants.
    if n_variant_rows > 0 and n_feature_rows > 0:
        sep_y = (n_variant_rows - 0.5) * row_step
        ax.axhline(sep_y, color="#bbbbbb", linewidth=0.8, linestyle="--")

    ax.set_xlim(-0.5, protein_length - 0.5)
    ax.set_ylim(-row_step * 0.6, (n_rows - 0.4) * row_step)
    ax.set_yticks(y_centers)
    ax.set_yticklabels(all_labels, fontsize=_fs(9, fig_width))
    ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
    ax.set_title(
        f"{title}  —  UniProt features & variants",
        fontsize=_fs(15, fig_width), pad=10,
    )
    ax.tick_params(axis="x", labelsize=_fs(10, fig_width))

    # Combined legend: feature patches + variant markers.
    legend_handles = [
        Patch(facecolor=_FEATURE_COLORS.get(ft, "#999"), edgecolor="black",
              linewidth=0.5, label=ft.title())
        for ft in display_types
    ]
    for vt_name in variant_track_names:
        colour, marker = variant_meta[vt_name]
        legend_handles.append(
            Line2D(
                [0], [0], marker=marker, color="w", markerfacecolor=colour,
                markeredgecolor="black", markeredgewidth=0.4,
                markersize=_fs(7, fig_width), label=vt_name,
            )
        )
    n_legend = len(legend_handles)
    ncol = min(n_legend, 4) if n_legend > 6 else min(n_legend, 3)
    ax.legend(
        handles=legend_handles, loc="upper center",
        bbox_to_anchor=(0.5, -0.15), fontsize=_fs(8, fig_width),
        framealpha=0.9, ncol=ncol,
        borderaxespad=0, columnspacing=1.5, handletextpad=0.5,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved features plot: {path}")


def plot_structure(
    plddt: Optional[List[float]],
    pae: Optional[np.ndarray],
    dist_matrix: Optional[np.ndarray],
    protein_length: int,
    title: str,
    path: str,
    plddt_cutoff: float = 70.0,
    pae_cutoff: float = 15.0,
    dpi: int = 150,
) -> None:
    """Plot AlphaFold pLDDT, PAE and distance matrix for a protein.

    Creates a figure with up to three panels:
    - Top: pLDDT per residue (coloured by confidence band).
    - Middle: PAE binary map (red = above cutoff / filtered).
    - Bottom: Cα distance matrix heatmap.
    """
    has_plddt = plddt is not None
    has_pae = pae is not None
    has_dist = dist_matrix is not None
    n_panels = int(has_plddt) + int(has_pae) + int(has_dist)
    if n_panels == 0:
        return

    height_ratios = []
    if has_plddt:
        height_ratios.append(1)
    if has_pae:
        height_ratios.append(3)
    if has_dist:
        height_ratios.append(3)

    fig_width = max(14, protein_length / 45)
    hs = _height_scale(fig_width)
    fig_height = sum(height_ratios) * 2.5 * hs + 1.5 * hs
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(fig_width, fig_height),
        dpi=dpi,
        gridspec_kw={"height_ratios": height_ratios},
        squeeze=False,
    )
    axes = axes.ravel()
    panel_idx = 0

    # --- pLDDT panel ---
    if has_plddt:
        ax = axes[panel_idx]
        panel_idx += 1
        x = np.arange(len(plddt))
        scores = np.array(plddt[:protein_length])

        # Colour each bar by the AlphaFold confidence band.
        bar_colors = []
        for s in scores:
            if s >= 90:
                bar_colors.append("#0053D6")
            elif s >= 70:
                bar_colors.append("#65CBF3")
            elif s >= 50:
                bar_colors.append("#FFDB13")
            else:
                bar_colors.append("#FF7D45")

        ax.bar(x[:protein_length], scores, width=1.0, color=bar_colors, edgecolor="none")
        ax.axhline(plddt_cutoff, color="#cb181d", linewidth=1.2, linestyle="--", alpha=0.8)
        ax.text(
            protein_length * 0.99, plddt_cutoff + 1.5,
            f"pLDDT cutoff = {plddt_cutoff:g}",
            ha="right", va="bottom", fontsize=_fs(10, fig_width),
            color="#cb181d", fontweight="bold",
        )
        ax.set_xlim(-0.5, protein_length - 0.5)
        ax.set_ylim(0, 100)
        ax.set_ylabel("pLDDT", fontsize=_fs(12, fig_width))
        ax.set_title(f"{title}  —  AlphaFold structure metrics",
                     fontsize=_fs(15, fig_width), pad=10)
        ax.tick_params(labelsize=_fs(10, fig_width))

        # Band legend.
        from matplotlib.patches import Patch

        legend_patches = [
            Patch(color="#0053D6", label="Very high (≥90)"),
            Patch(color="#65CBF3", label="Confident (70–90)"),
            Patch(color="#FFDB13", label="Low (50–70)"),
            Patch(color="#FF7D45", label="Very low (<50)"),
        ]
        ax.legend(
            handles=legend_patches, loc="lower left",
            fontsize=_fs(9, fig_width), framealpha=0.9, ncol=4,
        )

    # --- PAE panel (binary: filtered vs kept) ---
    if has_pae:
        ax = axes[panel_idx]
        panel_idx += 1
        pae_sub = pae[:protein_length, :protein_length]

        # Binary mask: 1 = above cutoff (filtered / red), 0 = at or below (kept / white).
        filtered = (pae_sub > pae_cutoff).astype(float)
        binary_cmap = mcolors.ListedColormap(["white", "#cb181d"])
        ax.imshow(
            filtered,
            aspect="auto",
            cmap=binary_cmap,
            vmin=0,
            vmax=1,
            interpolation="nearest",
            origin="upper",
        )
        ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
        ax.set_ylabel("Residue index", fontsize=_fs(13, fig_width))
        pae_label = f"PAE cutoff = {pae_cutoff} Å"
        if not has_plddt:
            ax.set_title(f"{title}  —  {pae_label}",
                         fontsize=_fs(15, fig_width), pad=10)
        else:
            ax.set_title(pae_label, fontsize=_fs(13, fig_width), pad=6)
        ax.tick_params(labelsize=_fs(10, fig_width))

        # Legend instead of colorbar.
        from matplotlib.patches import Patch as _Patch

        ax.legend(
            handles=[
                _Patch(facecolor="#cb181d", edgecolor="#333", linewidth=0.5,
                       label=f"PAE > {pae_cutoff} Å (filtered)"),
                _Patch(facecolor="white", edgecolor="#333", linewidth=0.5,
                       label=f"PAE ≤ {pae_cutoff} Å (kept)"),
            ],
            loc="upper right", fontsize=_fs(10, fig_width), framealpha=0.9,
        )

    # --- Distance matrix panel (bottom) ---
    if has_dist:
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        ax = axes[panel_idx]
        panel_idx += 1
        dist_sub = dist_matrix[:protein_length, :protein_length]
        im = ax.imshow(
            dist_sub,
            aspect="auto",
            cmap="viridis_r",
            interpolation="nearest",
            origin="upper",
        )
        ax.set_ylabel("Residue index", fontsize=_fs(13, fig_width))
        ax.set_xlabel("Residue index", fontsize=_fs(13, fig_width))
        dist_title = "Cα distance matrix (Å)"
        if not has_plddt and not has_pae:
            ax.set_title(f"{title}  —  {dist_title}",
                         fontsize=_fs(15, fig_width), pad=10)
        else:
            ax.set_title(dist_title, fontsize=_fs(13, fig_width), pad=6)
        ax.tick_params(labelsize=_fs(10, fig_width))

        # Use a divider-based colorbar so the axes width stays aligned.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.08)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("Distance (Å)", fontsize=_fs(11, fig_width))
        cbar.ax.tick_params(labelsize=_fs(9, fig_width))

    # Match widths: add invisible colorbars to non-dist panels so they share
    # the same effective width as the distance-matrix panel.
    if has_dist:
        from mpl_toolkits.axes_grid1 import make_axes_locatable as _mal

        for i in range(panel_idx):
            target_ax = axes[i]
            if target_ax is ax:
                continue
            div = _mal(target_ax)
            spacer = div.append_axes("right", size="2%", pad=0.08)
            spacer.set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved structure plot: {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize Proemis3D constraint regions for a gene across methods."
    )
    parser.add_argument(
        "--tsv",
        required=True,
        help="Path to the union TSV produced by compare_forward_methods_report.py.",
    )
    parser.add_argument(
        "--gene",
        required=True,
        nargs="+",
        help="One or more UniProt IDs to visualize (e.g. P12345 P13637).",
    )
    parser.add_argument(
        "--style",
        choices=["all", "tracks", "bars", "heatmap"],
        default="all",
        help="Which plot style(s) to generate. Default: all.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for output PNGs. Default: current directory.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=0.0,
        help="Minimum OE upper CI for the colourmap. Default: 0.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=1.5,
        help="Maximum OE upper CI for the colourmap. Default: 1.5.",
    )
    parser.add_argument(
        "--sort",
        choices=["original", "method", "filter", "cluster"],
        default="filter",
        help=(
            "How to order methods on the y-axis. "
            "'filter' groups by filter type (no filter → pLDDT → PAE) with "
            "AIC before LRT within each group (default); "
            "'method' groups by region-selection method (AIC/LRT) then filter; "
            "'cluster' reorders by hierarchical clustering of residue-level "
            "OE patterns; 'original' keeps the TSV order."
        ),
    )
    parser.add_argument(
        "--af2-gcs-dir",
        default=None,
        help=(
            "GCS directory containing gzipped AlphaFold JSON files "
            "(e.g. gs://gnomad-julia/alphafold2). Files are expected as "
            "AF-{gene}-F1-confidence_v4.json.gz and "
            "AF-{gene}-F1-predicted_aligned_error_v4.json.gz."
        ),
    )
    parser.add_argument(
        "--plddt-json",
        default=None,
        help="Path to a local AlphaFold pLDDT confidence JSON file.",
    )
    parser.add_argument(
        "--pae-json",
        default=None,
        help="Path to a local AlphaFold PAE JSON file.",
    )
    parser.add_argument(
        "--plddt-cutoff",
        type=float,
        default=70.0,
        help="pLDDT cutoff line on the structure plot. Default: 70.",
    )
    parser.add_argument(
        "--pae-cutoff",
        type=float,
        default=15.0,
        help="PAE cutoff in Å for the structure plot. Pairs above this are coloured red. Default: 15.",
    )
    parser.add_argument(
        "--no-structure",
        action="store_true",
        help="Skip the structure (pLDDT + PAE + distance matrix) plot entirely.",
    )
    parser.add_argument(
        "--features-json",
        default=None,
        help=(
            "Path to a local JSON file with UniProt features. "
            "Format: {uniprot_id: {feature_type: [{start, end, note, ...}]}}."
        ),
    )
    parser.add_argument(
        "--features-ht",
        default=None,
        help=(
            "GCS path to a Hail Table with UniProt features keyed by "
            "(uniprot_id, residue_index). Features will be exported to a "
            "local JSON automatically (requires Hail batch backend)."
        ),
    )
    parser.add_argument(
        "--variants-tsv",
        default=None,
        help=(
            "Path to a TSV with autism, DD de novo and gnomAD v4 de novo "
            "variant annotations. Variants are shown as marker tracks in "
            "the features plot."
        ),
    )
    parser.add_argument(
        "--no-features",
        action="store_true",
        help="Skip the UniProt features plot entirely.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="DPI for output images. Default: 150.",
    )
    args = parser.parse_args()

    df = _read_and_normalize_tsv(args.tsv)
    os.makedirs(args.output_dir, exist_ok=True)

    styles = (
        ["tracks", "bars", "heatmap"] if args.style == "all" else [args.style]
    )
    plot_funcs = {
        "tracks": plot_tracks,
        "bars": plot_bars,
        "heatmap": plot_heatmap,
    }

    # --- Pre-load / export UniProt features ---
    all_features: Dict[str, Dict[str, List[dict]]] = {}
    if not args.no_features:
        if args.features_json:
            import json

            with open(args.features_json) as _fh:
                all_features = json.load(_fh)
        elif args.features_ht:
            feat_json_path = os.path.join(args.output_dir, "uniprot_features.json")
            _export_features_from_ht(args.features_ht, args.gene, feat_json_path)
            import json

            with open(feat_json_path) as _fh:
                all_features = json.load(_fh)

    for gene in args.gene:
        print(f"\n{'='*60}\nProcessing gene: {gene}\n{'='*60}")
        try:
            per_method, protein_length, run_labels = prepare_gene_data(df, gene)
        except Exception as exc:
            print(f"Skipping {gene}: {exc}")
            continue

        if args.sort == "filter":
            run_labels = sort_by_filter_group(run_labels)
        elif args.sort == "method":
            run_labels = sort_by_method(run_labels)
        elif args.sort == "cluster":
            run_labels = sort_by_cluster(per_method, protein_length, run_labels)

        gene_safe = gene.replace("/", "_")
        title = f"{gene}  ({protein_length} residues, {len(run_labels)} methods)"

        for style in styles:
            if style == "bars":
                # Generate grey catch-all version.
                out_path = os.path.join(args.output_dir, f"{gene_safe}_bars.png")
                plot_bars(
                    per_method=per_method,
                    protein_length=protein_length,
                    run_labels=run_labels,
                    title=title,
                    path=out_path,
                    vmin=args.vmin,
                    vmax=args.vmax,
                    dpi=args.dpi,
                    color_null=False,
                )
                # Generate OE-colored catch-all version.
                out_path_colored = os.path.join(
                    args.output_dir, f"{gene_safe}_bars_colored.png"
                )
                plot_bars(
                    per_method=per_method,
                    protein_length=protein_length,
                    run_labels=run_labels,
                    title=title,
                    path=out_path_colored,
                    vmin=args.vmin,
                    vmax=args.vmax,
                    dpi=args.dpi,
                    color_null=True,
                )
            else:
                out_path = os.path.join(args.output_dir, f"{gene_safe}_{style}.png")
                plot_funcs[style](
                    per_method=per_method,
                    protein_length=protein_length,
                    run_labels=run_labels,
                    title=title,
                    path=out_path,
                    vmin=args.vmin,
                    vmax=args.vmax,
                    dpi=args.dpi,
                )

        # --- Structure plot (pLDDT + PAE + distance matrix) ---
        if not args.no_structure:
            import tempfile

            plddt_data: Optional[List[float]] = None
            pae_data: Optional[np.ndarray] = None
            dist_data: Optional[np.ndarray] = None

            # Priority 1: local JSON files (only usable for single-gene runs)
            if args.plddt_json:
                plddt_data = _load_plddt(args.plddt_json)
            if args.pae_json:
                pae_data = _load_pae(args.pae_json)

            # Priority 2: GCS bucket
            if (plddt_data is None or pae_data is None or dist_data is None) and args.af2_gcs_dir:
                gcs_plddt, gcs_pae, gcs_dist = _load_af2_from_gcs(
                    args.af2_gcs_dir, gene,
                )
                if plddt_data is None:
                    plddt_data = gcs_plddt
                if pae_data is None:
                    pae_data = gcs_pae
                if dist_data is None:
                    dist_data = gcs_dist

            # Priority 3: AlphaFold DB API (pLDDT and PAE only; no dist)
            if plddt_data is None or pae_data is None:
                for kind in ["plddt", "pae"]:
                    if kind == "plddt" and plddt_data is not None:
                        continue
                    if kind == "pae" and pae_data is not None:
                        continue
                    raw = _fetch_alphafold_json(gene, kind)
                    if raw:
                        tmp = os.path.join(
                            tempfile.gettempdir(), f"{gene_safe}_{kind}.json"
                        )
                        with open(tmp, "w") as fh:
                            fh.write(raw)
                        if kind == "plddt":
                            plddt_data = _load_plddt(tmp)
                        else:
                            pae_data = _load_pae(tmp)

            if plddt_data is not None or pae_data is not None or dist_data is not None:
                out_path = os.path.join(args.output_dir, f"{gene_safe}_structure.png")
                plot_structure(
                    plddt=plddt_data,
                    pae=pae_data,
                    dist_matrix=dist_data,
                    protein_length=protein_length,
                    title=gene,
                    path=out_path,
                    plddt_cutoff=args.plddt_cutoff,
                    pae_cutoff=args.pae_cutoff,
                    dpi=args.dpi,
                )
            else:
                print(f"No structure data available for {gene}; skipping structure plot.")

        # --- UniProt features plot ---
        if not args.no_features and (gene in all_features or args.variants_tsv):
            gene_features = all_features.get(gene, {})
            gene_variants: Dict[str, List[int]] = {}
            if args.variants_tsv:
                gene_variants = _load_variants(args.variants_tsv, gene)
                if gene_variants:
                    print(f"Loaded variants: {', '.join(f'{k} ({len(v)})' for k, v in gene_variants.items())}")
            out_path = os.path.join(args.output_dir, f"{gene_safe}_features.png")
            plot_features(
                features=gene_features,
                protein_length=protein_length,
                title=gene,
                path=out_path,
                variants=gene_variants,
                dpi=args.dpi,
            )


if __name__ == "__main__":
    main()
