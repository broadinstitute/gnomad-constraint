#!/usr/bin/env python3
"""
Compare Proemis3D forward outputs across filtering methods and model comparison (AIC, AIC weight, LRT).

Two modes:
1. Export mode: read all HTs from config, add run_label + method/cutoff columns, union, export to TSV.
2. Report from TSV: read combined TSV with pandas, aggregate by run_label, write Markdown + CSV report.

Usage:
  # Create combined TSV from config, then write report (one Hail job + pandas).
  # With Hail batch backend, --output-tsv must be a remote URI (e.g. gs://bucket/forward_combined.tsv).
  python compare_forward_methods_report.py --config runs.yaml --output-tsv gs://your-bucket/forward_combined.tsv --output report.md --csv summary.csv

  # Report from existing TSV (no Hail; fast)
  python compare_forward_methods_report.py --tsv combined.tsv --output report.md --csv summary.csv

  # Limit runs when building TSV (quick test)
  python compare_forward_methods_report.py --config runs.yaml --output-tsv combined.tsv --output report.md --limit 2

  # Run export on Dataproc (--backend spark is default)
  gcloud dataproc jobs submit pyspark compare_forward_methods_report.py --cluster jg3 --region us-central1 -- \\
    --config forward_runs_144_genes.yaml --output-tsv gs://your-bucket/forward_combined.tsv --output report.md --csv summary.csv

Config YAML format (runs.yaml):
  runs:
    - path: "gs://bucket/forward.aic.ht"
      label: "AIC (no pLDDT/PAE)"

Comparing debug (union) TSV vs single-run export:
  - The union TSV has one row per (uniprot_id, transcript_id, residue_index) *per run_label*.
  - To compare to a single-run export (e.g. one gene's .tsv from a 158-gene run): filter the union
    to the run_label that matches the single run (same pLDDT/PAE method and cutoffs), then filter
    to the same uniprot_id and transcript_id. Row counts and values should match only if both
    exports came from the same underlying forward HT (same pipeline run and inputs).
  - If numbers differ: the two TSVs likely came from different pipeline runs (e.g. 144-gene vs
    158-gene), so inputs (obs/exp, AF2, PAE, pLDDT) may differ; or you are comparing different
    run_labels (e.g. different PAE/pLDDT method).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


def _cluster_heatmap_order(
    matrix: np.ndarray,
) -> Tuple[np.ndarray, List[int], List[int], Any, Any]:
    """Apply hierarchical clustering to rows and columns.
    Returns (reordered_matrix, row_order, col_order, link_row, link_col) for dendrograms.
    link_row/link_col are None when there is only one row/column (no clustering in that dimension).
    """
    link_row, link_col = None, None
    fill = float(np.nanmean(matrix)) if np.any(np.isfinite(matrix)) else 0.0
    X = np.nan_to_num(matrix, nan=fill, posinf=fill, neginf=fill)
    Xcol = np.nan_to_num(matrix.T, nan=fill, posinf=fill, neginf=fill)
    row_order = list(range(matrix.shape[0]))
    col_order = list(range(matrix.shape[1]))
    if matrix.shape[0] > 1:
        link_row = linkage(pdist(X, metric="euclidean"), method="average")
        d = dendrogram(link_row, no_plot=True)
        row_order = d["leaves"]
    if matrix.shape[1] > 1:
        link_col = linkage(pdist(Xcol, metric="euclidean"), method="average")
        d = dendrogram(link_col, no_plot=True)
        col_order = d["leaves"]
    reordered = matrix[np.ix_(row_order, col_order)]
    return reordered, row_order, col_order, link_row, link_col


def _plot_clustered_heatmap(
    matrix: np.ndarray,
    row_labels: List[str],
    col_labels: List[str],
    path: str,
    title: str,
    cbar_label: str,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Tuple[float, float] = (12, 10),
    dpi: int = 150,
) -> None:
    """Draw a heatmap with hierarchical clustering and row/column dendrograms.

    Dendrograms are aligned with the heatmap by explicitly setting axis limits:
    scipy ``dendrogram`` places leaves at x = 5, 15, 25, ... (spacing 10), so for
    *n* leaves the full span is [0, n*10].  ``imshow`` places pixel centres at
    0, 1, 2, ... so the extent is [-0.5, n-0.5].  Both axes share the same
    physical GridSpec width/height, so we set matching limits.

    The colorbar is given its own GridSpec column so it does not steal width from
    the heatmap (which would misalign it with the column dendrogram above).
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    reordered, row_order, col_order, link_row, link_col = _cluster_heatmap_order(matrix)
    ordered_row_labels = [row_labels[i] for i in row_order]
    ordered_col_labels = [col_labels[i] for i in col_order]
    n_rows, n_cols = reordered.shape
    dendro_w = 1.5
    dendro_h = 1.2
    cbar_w = 0.3
    has_dendro = link_row is not None and link_col is not None
    if has_dendro:
        fig = plt.figure(
            figsize=(figsize[0] + dendro_w + cbar_w, figsize[1] + dendro_h)
        )
        gs = GridSpec(
            2,
            3,
            figure=fig,
            width_ratios=[dendro_w, figsize[0], cbar_w],
            height_ratios=[dendro_h, figsize[1]],
            hspace=0.02,
            wspace=0.02,
        )
        ax_row = fig.add_subplot(gs[1, 0])
        ax_col = fig.add_subplot(gs[0, 1])
        ax_heat = fig.add_subplot(gs[1, 1])
        ax_cbar = fig.add_subplot(gs[1, 2])

        # Row dendrogram (left).  Leaf positions: y = 5, 15, ... so span [0, n_rows*10].
        # imshow origin='upper' puts row 0 at top, so invert the dendrogram y-axis.
        dendrogram(
            link_row,
            ax=ax_row,
            orientation="left",
            leaf_rotation=0,
            color_threshold=0,
            above_threshold_color="0.3",
        )
        ax_row.set_ylim(0, n_rows * 10)
        ax_row.invert_yaxis()
        ax_row.set_axis_off()

        # Column dendrogram (top).  Leaf positions: x = 5, 15, ... so span [0,
        # n_cols*10].
        dendrogram(
            link_col,
            ax=ax_col,
            orientation="top",
            leaf_rotation=90,
            color_threshold=0,
            above_threshold_color="0.3",
        )
        ax_col.set_xlim(0, n_cols * 10)
        ax_col.set_axis_off()

        # Hide the unused top-left and top-right cells.
        for pos in [(0, 0), (0, 2)]:
            ax_empty = fig.add_subplot(gs[pos[0], pos[1]])
            ax_empty.set_axis_off()
    else:
        fig = plt.figure(figsize=(figsize[0] + cbar_w, figsize[1]))
        gs = GridSpec(1, 2, figure=fig, width_ratios=[figsize[0], cbar_w], wspace=0.02)
        ax_heat = fig.add_subplot(gs[0, 0])
        ax_cbar = fig.add_subplot(gs[0, 1])

    kwargs = {"aspect": "auto", "cmap": cmap, "interpolation": "nearest"}
    if vmin is not None:
        kwargs["vmin"] = vmin
    if vmax is not None:
        kwargs["vmax"] = vmax
    im = ax_heat.imshow(reordered, **kwargs)
    ax_heat.set_xticks(np.arange(n_cols))
    ax_heat.set_xticklabels(
        ordered_col_labels, rotation=90, ha="right", fontsize=max(4, 10 - n_cols // 20)
    )
    ax_heat.set_yticks(np.arange(n_rows))
    ax_heat.set_yticklabels(ordered_row_labels, fontsize=max(4, 9 - n_rows // 30))
    ax_heat.set_xlabel("Run")
    ax_heat.set_ylabel("Gene (Uniprot ID)")
    ax_heat.set_title(title)
    fig.colorbar(im, cax=ax_cbar, label=cbar_label)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


# Path parser: no Hail/pandas required
def _parse_forward_path(path: str) -> Dict[str, Optional[str]]:
    """Parse forward HT path into model_comparison, pae/plddt method and cutoff."""
    name = Path(path).stem  # no .ht
    out = {
        "model_comparison": None,
        "pae_cutoff": None,
        "pae_method": None,
        "plddt_cutoff": None,
        "plddt_method": None,
    }
    m = re.search(r"\.(aic|aic_weight|lrt)(?:\.|$)", name)
    if m:
        out["model_comparison"] = m.group(1)
    m = re.search(r"\.pae_cutoff_([0-9.]+)", name)
    if m:
        out["pae_cutoff"] = m.group(1)
    m = re.search(r"\.pae_method_([^.]+)", name)
    if m:
        out["pae_method"] = m.group(1)
    m = re.search(r"\.plddt_cutoff_([0-9.]+)", name)
    if m:
        out["plddt_cutoff"] = m.group(1)
    m = re.search(r"\.plddt_method_([^.]+)", name)
    if m:
        out["plddt_method"] = m.group(1)
    return out


def _model_plddt_pae_display(
    model_comparison: Any,
    plddt_method: Any,
    pae_method: Any,
) -> Tuple[str, str, str]:
    """Map parsed method names to short display for summary table. Returns (model, plddt, pae)."""

    def _ok(x: Any) -> bool:
        if x is None:
            return False
        try:
            import pandas as pd

            if pd.isna(x):
                return False
        except Exception:
            pass
        s = str(x).strip().lower()
        return bool(s) and s not in ("none", "nan", "<na>")

    model = str(model_comparison).strip() if _ok(model_comparison) else None
    plddt = str(plddt_method).strip() if _ok(plddt_method) else None
    pae = str(pae_method).strip() if _ok(pae_method) else None

    model_display = (
        "AIC"
        if model == "aic"
        else (
            "AIC weight"
            if model == "aic_weight"
            else "LRT" if model == "lrt" else (model or "–")
        )
    )
    plddt_display = (
        (
            "remove"
            if plddt == "remove_low_plddt_residues"
            else "exclude" if plddt == "exclude_low_plddt_from_stats" else "–"
        )
        if plddt
        else "–"
    )
    # PAE method: truncate_on_pairwise_pae_with_center -> truncate_center, etc.
    pae_display = "–"
    if pae:
        if (
            "truncate_on_pairwise_pae_with_center" in pae
            or pae == "truncate_on_pairwise_pae_with_center"
        ):
            pae_display = "truncate_center"
        elif (
            "filter_on_pairwise_pae_with_center" in pae
            or pae == "filter_on_pairwise_pae_with_center"
        ):
            pae_display = "filter_center"
        elif (
            "filter_on_pairwise_pae_in_region" in pae
            or pae == "filter_on_pairwise_pae_in_region"
        ):
            pae_display = "filter_region"
        elif (
            "exclude_on_pairwise_pae_with_center" in pae
            or pae == "exclude_on_pairwise_pae_with_center"
        ):
            pae_display = "exclude_center"
        elif (
            "exclude_on_pairwise_pae_in_region" in pae
            or pae == "exclude_on_pairwise_pae_in_region"
        ):
            pae_display = "exclude_region"
        else:
            pae_display = pae
    return model_display, plddt_display, pae_display


def load_config(config_path: str) -> List[Tuple[str, str]]:
    """Load (path, label) pairs from YAML or JSON config."""
    with open(config_path) as f:
        raw = f.read()
    if config_path.endswith(".json"):
        data = json.loads(raw)
    else:
        try:
            import yaml

            data = yaml.safe_load(raw)
        except ImportError:
            raise ImportError("PyYAML required for YAML config. pip install pyyaml")
    runs = data.get("runs", data)
    return [(r["path"], r["label"]) for r in runs]


# ---------- Hail: load HTs, add columns, union, export to TSV ----------
def _ensure_hl(backend: str = "spark") -> None:
    import hail as hl

    try:
        if backend == "batch":
            hl.init(backend="batch", tmp_dir="gs://gnomad-tmp-4day")
        else:
            # Spark backend (default); use when running on Dataproc or locally with
            # Spark
            hl.init(tmp_dir="gs://gnomad-tmp-4day")
    except Exception:
        pass


def export_union_ht_to_tsv(
    runs: List[Tuple[str, str]],
    tsv_path: str,
    limit: Optional[int] = None,
    backend: str = "spark",
) -> None:
    """Load each HT, add run_label and method/cutoff columns, union, export to TSV."""
    import hail as hl

    # Hail batch backend only supports remote paths (e.g. GCS)
    if backend == "batch" and not (
        tsv_path.startswith("gs://") or tsv_path.startswith("hdfs://")
    ):
        raise ValueError(
            "The Hail batch backend requires a remote path for --output-tsv. "
            "Use a GCS URI instead, e.g.: --output-tsv gs://your-bucket/forward_combined.tsv"
        )

    _ensure_hl(backend=backend)
    if limit is not None:
        runs = runs[:limit]
        print(f"Limited to first {len(runs)} run(s) (--limit {limit}).")

    def _opt_str(x: Optional[str]):
        return hl.literal(x) if x is not None else hl.literal("None")

    tables = []
    for i, (path, label) in enumerate(runs):
        print(f"Loading {i+1}/{len(runs)}: {path}")
        ht = hl.read_table(path).naive_coalesce(1).cache()
        parsed = _parse_forward_path(path)
        ht = ht.annotate(
            run_label=hl.literal(label),
            model_comparison=_opt_str(parsed["model_comparison"]),
            pae_cutoff=_opt_str(parsed["pae_cutoff"]),
            pae_method=_opt_str(parsed["pae_method"]),
            plddt_cutoff=_opt_str(parsed["plddt_cutoff"]),
            plddt_method=_opt_str(parsed["plddt_method"]),
        ).cache()
        tables.append(ht)

    print("Unioning tables...")
    union_ht = tables[0]
    for t in tables[1:]:
        union_ht = union_ht.union(t)
    union_ht = union_ht.cache()

    print(f"Exporting to {tsv_path}...")
    union_ht.export(tsv_path, header=True)
    print(f"Wrote {tsv_path}")


# ---------- Pandas: read TSV, aggregate by run_label, produce stats ----------
def _apply_gene_mapping(df: Any, mapping_tsv_path: str) -> Any:
    """Add gene_id to df by merging with a TSV that has transcript_id and gene_id. Returns a DataFrame with unique columns and gene_id added (for one-per-gene table)."""
    import pandas as pd

    df = _df_unique_columns(df)
    if "transcript_id" not in df.columns:
        return df
    mapping = pd.read_csv(mapping_tsv_path, sep="\t", low_memory=False)
    if "gene_id" not in mapping.columns or "transcript_id" not in mapping.columns:
        return df
    mapping = mapping[["transcript_id", "gene_id"]].drop_duplicates(
        subset=["transcript_id"]
    )
    return df.merge(mapping, on="transcript_id", how="left")


def _read_and_normalize_tsv(tsv_path: str):
    """Read combined TSV and normalize columns (e.g. _is_null). Returns DataFrame."""
    import pandas as pd

    print(f"Reading {tsv_path}...")
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    # Ensure unique column names (TSV may have duplicate headers; pandas then
    # raises "column label is not unique")
    if df.columns.duplicated().any():
        seen = {}
        new_cols = []
        for c in df.columns:
            if c in seen:
                seen[c] += 1
                new_cols.append(f"{c}.{seen[c]}")
            else:
                seen[c] = 0
                new_cols.append(c)
        df.columns = new_cols
    if "run_label" not in df.columns:
        raise ValueError("TSV must have run_label column (from --output-tsv export).")
    if "is_null" in df.columns:
        df["_is_null"] = (
            df["is_null"].fillna(False).astype(str).str.lower().isin(("true", "1"))
        )
    else:
        df["_is_null"] = False
    return df


def _first_column_index(df: Any, name: str) -> int:
    """Return the first column index with the given name. Raises KeyError if not found."""
    for i, c in enumerate(df.columns):
        if c == name:
            return i
    raise KeyError(name)


def _df_unique_columns(df: Any) -> Any:
    """Return a copy of the DataFrame with only the first occurrence of each column name (avoids 'column label is not unique' errors)."""
    import pandas as pd

    # Detect duplicates by count; duplicated() can be unreliable with some Index types
    if len(df.columns) == len(set(df.columns)):
        return df
    keep_idx = []
    seen = set()
    for i, c in enumerate(df.columns):
        if c not in seen:
            keep_idx.append(i)
            seen.add(c)
    return df.iloc[:, keep_idx].copy()


def restrict_to_one_transcript_per_gene(df: Any) -> Any:
    """Restrict DataFrame to one (uniprot_id, transcript_id) per gene. Uses gene_id when present (yields one transcript per gene, e.g. 144 genes); otherwise one per uniprot_id. Deterministic: first transcript when sorted."""
    import pandas as pd

    df = _df_unique_columns(df)
    try:
        uniprot_idx = _first_column_index(df, "uniprot_id")
        transcript_idx = _first_column_index(df, "transcript_id")
    except KeyError:
        return df
    # Prefer gene_id so we get one transcript per gene (e.g. 144 genes when
    # multiple transcripts share a gene)
    try:
        gene_idx = _first_column_index(df, "gene_id")
        use_gene_id = True
    except KeyError:
        gene_idx = None
        use_gene_id = False
    if use_gene_id:
        sub = df.iloc[:, [gene_idx, uniprot_idx, transcript_idx]].copy()
        sub.columns = ["gene_id", "uniprot_id", "transcript_id"]
        group_col = "gene_id"
    else:
        sub = df.iloc[:, [uniprot_idx, transcript_idx]].copy()
        sub.columns = ["uniprot_id", "transcript_id"]
        group_col = "uniprot_id"
    one_per_gene = (
        sub.drop_duplicates()
        .sort_values(["uniprot_id", "transcript_id"])
        .groupby(group_col, as_index=False)
        .first()
    )
    # df was deduplicated above; merge on the key columns
    return df.merge(
        one_per_gene[["uniprot_id", "transcript_id"]],
        on=["uniprot_id", "transcript_id"],
        how="inner",
    )


def compute_stats_from_tsv(tsv_path_or_df: Any) -> List[Tuple[str, Dict[str, Any]]]:
    """Read combined TSV (or use provided DataFrame) and compute per-run stats. Returns list of (label, stats_dict)."""
    import pandas as pd

    if isinstance(tsv_path_or_df, str):
        df = _read_and_normalize_tsv(tsv_path_or_df)
    else:
        df = tsv_path_or_df

    # Treat "None" (string from Hail export) and empty/NaN as no filter when
    # detecting baseline
    def _method_empty(val) -> bool:
        if pd.isna(val) or val == "":
            return True
        return isinstance(val, str) and val.strip() == "None"

    # Baseline total residues from unfiltered runs (no pLDDT/PAE) per model_comparison
    baseline_total_by_model = {}
    for run_label, grp in df.groupby("run_label", sort=False):
        first = grp.iloc[0]
        model = first.get("model_comparison")
        if (
            pd.isna(model)
            or model == ""
            or (isinstance(model, str) and model.strip() == "None")
        ):
            continue
        is_unfiltered = _method_empty(first.get("pae_method")) and _method_empty(
            first.get("plddt_method")
        )
        if is_unfiltered:
            baseline_total_by_model[model] = len(grp)

    results = []
    for run_label, grp in df.groupby("run_label", sort=False):
        # Per (uniprot_id, transcript_id) stats
        per_tx = (
            grp.groupby(["uniprot_id", "transcript_id"], dropna=False)
            .agg(
                n_residues=("residue_index", "count"),
                n_residues_missing_obs=("obs", lambda s: s.isna().sum()),
                n_residues_missing_exp=("exp", lambda s: s.isna().sum()),
                n_residues_in_null=("_is_null", "sum"),
            )
            .reset_index()
        )
        n_regions = (
            grp[~grp["_is_null"]]
            .groupby(["uniprot_id", "transcript_id"])["region_index"]
            .nunique()
        )
        per_tx = per_tx.merge(
            n_regions.rename("n_constraint_regions"),
            left_on=["uniprot_id", "transcript_id"],
            right_index=True,
            how="left",
        )
        per_tx["n_constraint_regions"] = (
            per_tx["n_constraint_regions"].fillna(0).astype(int)
        )

        n_genes = per_tx["uniprot_id"].nunique()
        n_transcripts = len(per_tx)
        total_residues = per_tx["n_residues"].sum()
        total_missing_obs = per_tx["n_residues_missing_obs"].sum()
        total_missing_exp = per_tx["n_residues_missing_exp"].sum()
        n_tx_all_null = (per_tx["n_constraint_regions"] == 0).sum()

        mean_regions = per_tx["n_constraint_regions"].mean()
        median_regions = per_tx["n_constraint_regions"].median()
        min_regions = per_tx["n_constraint_regions"].min()
        max_regions = per_tx["n_constraint_regions"].max()
        mean_null_frac = (
            per_tx["n_residues_in_null"] / per_tx["n_residues"].replace(0, pd.NA)
        ).mean()
        mean_residues_per_tx = per_tx["n_residues"].mean()

        # Region lengths (constraint regions only)
        region_sizes = (
            grp[~grp["_is_null"]]
            .groupby(["uniprot_id", "transcript_id", "region_index"])
            .size()
            .rename("region_length")
        )
        if len(region_sizes) > 0:
            mean_region_length = region_sizes.mean()
            median_region_length = region_sizes.median()
            min_region_length = region_sizes.min()
            max_region_length = region_sizes.max()
        else:
            mean_region_length = None
            median_region_length = None
            min_region_length = None
            max_region_length = None

        total_residues_val = int(total_residues)
        n_transcripts_val = int(n_transcripts)
        # Filtered out % = residues not in this run vs unfiltered baseline (same model)
        model = grp.iloc[0].get("model_comparison")
        baseline_total = (
            baseline_total_by_model.get(model) if pd.notna(model) and model else None
        )
        if baseline_total and baseline_total > 0:
            frac_filtered_out = (baseline_total - total_residues_val) / baseline_total
        else:
            frac_filtered_out = 0.0
            baseline_total = None

        # Per-run distribution data for distribution plots (not written to CSV)
        null_frac_per_tx = (
            per_tx["n_residues_in_null"] / per_tx["n_residues"].replace(0, pd.NA)
        ).dropna()
        distributions = {
            "n_regions_per_tx": per_tx["n_constraint_regions"].values,
            "region_lengths": (
                region_sizes.values
                if len(region_sizes) > 0
                else np.array([], dtype=float)
            ),
            "null_frac_per_tx": null_frac_per_tx.values,
        }

        first = grp.iloc[0]
        model_display, plddt_display, pae_display = _model_plddt_pae_display(
            first.get("model_comparison"),
            first.get("plddt_method"),
            first.get("pae_method"),
        )

        results.append(
            (
                run_label,
                {
                    "model_display": model_display,
                    "plddt_display": plddt_display,
                    "pae_display": pae_display,
                    "n_genes": int(n_genes),
                    "n_transcripts": n_transcripts_val,
                    "total_residues": total_residues_val,
                    "baseline_total_residues": baseline_total,
                    "frac_residues_filtered_out": frac_filtered_out,
                    "total_residues_missing_obs": int(total_missing_obs),
                    "total_residues_missing_exp": int(total_missing_exp),
                    "frac_residues_missing_obs": (
                        total_missing_obs / total_residues if total_residues else 0
                    ),
                    "frac_residues_missing_exp": (
                        total_missing_exp / total_residues if total_residues else 0
                    ),
                    "n_constraint_regions_mean": float(mean_regions),
                    "n_constraint_regions_median": float(median_regions),
                    "n_constraint_regions_min": int(min_regions),
                    "n_constraint_regions_max": int(max_regions),
                    "region_length_mean": (
                        float(mean_region_length)
                        if mean_region_length is not None
                        else None
                    ),
                    "region_length_median": (
                        float(median_region_length)
                        if median_region_length is not None
                        else None
                    ),
                    "region_length_min": (
                        int(min_region_length)
                        if min_region_length is not None
                        else None
                    ),
                    "region_length_max": (
                        int(max_region_length)
                        if max_region_length is not None
                        else None
                    ),
                    "frac_residues_in_null_mean": float(mean_null_frac),
                    "mean_residues_per_transcript": float(mean_residues_per_tx),
                    "n_transcripts_all_null": int(n_tx_all_null),
                    "frac_transcripts_all_null": (
                        n_tx_all_null / n_transcripts_val if n_transcripts_val else 0
                    ),
                    "_distributions": distributions,
                },
            )
        )
    return results


def compute_per_gene_run_matrices(
    tsv_path_or_df: Any,
) -> Tuple[Any, Any, List[str], List[str]]:
    """Build per-gene, per-run matrices for heatmaps.
    Returns (matrix_regions, matrix_null_frac, run_labels, gene_ids).
    matrix_regions[i, j] = mean constraint regions per transcript for gene i in run j.
    matrix_null_frac[i, j] = mean null fraction for gene i in run j.
    """
    import pandas as pd

    if isinstance(tsv_path_or_df, str):
        df = _read_and_normalize_tsv(tsv_path_or_df)
    else:
        df = tsv_path_or_df

    rows = []
    for run_label, grp in df.groupby("run_label", sort=False):
        per_tx = (
            grp.groupby(["uniprot_id", "transcript_id"], dropna=False)
            .agg(
                n_residues=("residue_index", "count"),
                n_residues_in_null=("_is_null", "sum"),
            )
            .reset_index()
        )
        n_regions = (
            grp[~grp["_is_null"]]
            .groupby(["uniprot_id", "transcript_id"])["region_index"]
            .nunique()
        )
        per_tx = per_tx.merge(
            n_regions.rename("n_constraint_regions"),
            left_on=["uniprot_id", "transcript_id"],
            right_index=True,
            how="left",
        )
        per_tx["n_constraint_regions"] = (
            per_tx["n_constraint_regions"].fillna(0).astype(int)
        )
        per_tx["null_frac"] = per_tx["n_residues_in_null"] / per_tx[
            "n_residues"
        ].replace(0, pd.NA)
        per_gene = (
            per_tx.groupby("uniprot_id", dropna=False)
            .agg(
                mean_regions=("n_constraint_regions", "mean"),
                mean_null_frac=("null_frac", "mean"),
            )
            .reset_index()
        )
        per_gene["run_label"] = run_label
        rows.append(per_gene)

    gene_run = pd.concat(rows, ignore_index=True)
    run_labels = gene_run["run_label"].unique().tolist()
    gene_ids = gene_run["uniprot_id"].unique().tolist()

    matrix_regions = np.full((len(gene_ids), len(run_labels)), np.nan)
    matrix_null_frac = np.full((len(gene_ids), len(run_labels)), np.nan)
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    run_to_idx = {r: i for i, r in enumerate(run_labels)}

    for _, row in gene_run.iterrows():
        i = gene_to_idx[row["uniprot_id"]]
        j = run_to_idx[row["run_label"]]
        matrix_regions[i, j] = row["mean_regions"]
        matrix_null_frac[i, j] = (
            row["mean_null_frac"] if pd.notna(row["mean_null_frac"]) else np.nan
        )

    return matrix_regions, matrix_null_frac, run_labels, gene_ids


def compute_gene_examples(
    tsv_path_or_df: Any,
    n_stable: int = 1,
    n_outliers: int = 2,
) -> Optional[Dict[str, List[Dict[str, Any]]]]:
    """Identify genes that are stable across methods vs outliers for different reasons.
    Returns dict with keys: stable, outlier_region_variance, outlier_often_all_null, outlier_null_frac_range.
    Each value is a list of dicts with gene (uniprot_id), reason summary, and optional stats.
    """
    import pandas as pd

    if isinstance(tsv_path_or_df, str):
        df = _read_and_normalize_tsv(tsv_path_or_df)
    else:
        df = tsv_path_or_df

    # Per (run_label, uniprot_id): mean n_constraint_regions, mean null_frac,
    # any transcript all-null
    rows = []
    for run_label, grp in df.groupby("run_label", sort=False):
        per_tx = (
            grp.groupby(["uniprot_id", "transcript_id"], dropna=False)
            .agg(
                n_residues=("residue_index", "count"),
                n_residues_in_null=("_is_null", "sum"),
            )
            .reset_index()
        )
        n_regions = (
            grp[~grp["_is_null"]]
            .groupby(["uniprot_id", "transcript_id"])["region_index"]
            .nunique()
        )
        per_tx = per_tx.merge(
            n_regions.rename("n_constraint_regions"),
            left_on=["uniprot_id", "transcript_id"],
            right_index=True,
            how="left",
        )
        per_tx["n_constraint_regions"] = (
            per_tx["n_constraint_regions"].fillna(0).astype(int)
        )
        per_tx["null_frac"] = per_tx["n_residues_in_null"] / per_tx[
            "n_residues"
        ].replace(0, pd.NA)

        per_gene = (
            per_tx.groupby("uniprot_id", dropna=False)
            .agg(
                mean_regions=("n_constraint_regions", "mean"),
                mean_null_frac=("null_frac", "mean"),
                any_all_null=("n_constraint_regions", lambda s: (s == 0).any()),
            )
            .reset_index()
        )
        per_gene["run_label"] = run_label
        rows.append(per_gene)

    gene_run = pd.concat(rows, ignore_index=True)

    # Per-gene cross-run stats
    agg = (
        gene_run.groupby("uniprot_id", dropna=False)
        .agg(
            std_mean_regions=("mean_regions", "std"),
            mean_regions_median=("mean_regions", "median"),
            frac_runs_any_all_null=("any_all_null", "mean"),
            range_null_frac=(
                "mean_null_frac",
                lambda s: s.max() - s.min() if s.notna().any() else 0,
            ),
            n_runs=("run_label", "count"),
        )
        .reset_index()
    )

    # Require gene present in all runs (no missing runs)
    n_runs = gene_run["run_label"].nunique()
    agg = agg[agg["n_runs"] == n_runs].copy()
    agg["std_mean_regions"] = agg["std_mean_regions"].fillna(0)

    out: Dict[str, List[Dict[str, Any]]] = {
        "stable": [],
        "outlier_region_variance": [],
        "outlier_often_all_null": [],
        "outlier_null_frac_range": [],
    }

    # Stable: low std(mean_regions), never/seldom all-null
    stable_candidates = agg[agg["frac_runs_any_all_null"] < 0.1].copy()
    stable_candidates = stable_candidates.sort_values(
        "std_mean_regions", ascending=True
    )
    for _, row in stable_candidates.head(n_stable).iterrows():
        out["stable"].append(
            {
                "gene": row["uniprot_id"],
                "std_mean_regions": float(row["std_mean_regions"]),
                "frac_runs_any_all_null": float(row["frac_runs_any_all_null"]),
            }
        )

    # Outlier: region count varies most
    by_std = agg.sort_values("std_mean_regions", ascending=False)
    for _, row in by_std.head(n_outliers).iterrows():
        out["outlier_region_variance"].append(
            {
                "gene": row["uniprot_id"],
                "std_mean_regions": float(row["std_mean_regions"]),
                "mean_regions_median": float(row["mean_regions_median"]),
            }
        )

    # Outlier: often all-null across runs
    by_all_null = agg.sort_values("frac_runs_any_all_null", ascending=False)
    for _, row in by_all_null.head(n_outliers).iterrows():
        if row["frac_runs_any_all_null"] > 0:
            out["outlier_often_all_null"].append(
                {
                    "gene": row["uniprot_id"],
                    "frac_runs_any_all_null": float(row["frac_runs_any_all_null"]),
                }
            )

    # Outlier: null fraction range across runs
    by_range = agg.sort_values("range_null_frac", ascending=False)
    for _, row in by_range.head(n_outliers).iterrows():
        if row["range_null_frac"] > 0:
            out["outlier_null_frac_range"].append(
                {
                    "gene": row["uniprot_id"],
                    "range_null_frac": float(row["range_null_frac"]),
                }
            )

    return out


def get_gene_example_details(
    tsv_path_or_df: Any,
    gene_ids: List[str],
) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """
    For each gene_id, get per-run, per-transcript region details: n_regions, region lengths and
    residue ranges (start, end) per region, null span, null count.
    Returns: { gene_id: { run_label: [ { transcript_id, n_regions, region_lengths, region_ranges,
              null_range, n_null, n_residues }, ... ], ... }, ... }
    """
    import pandas as pd

    if isinstance(tsv_path_or_df, str):
        df = _read_and_normalize_tsv(tsv_path_or_df)
    else:
        df = tsv_path_or_df

    df_genes = df[df["uniprot_id"].isin(gene_ids)]
    out: Dict[str, Dict[str, List[Dict[str, Any]]]] = {g: {} for g in gene_ids}

    for run_label, grp in df_genes.groupby("run_label", sort=False):
        for (uniprot_id, transcript_id), tx_grp in grp.groupby(
            ["uniprot_id", "transcript_id"], dropna=False
        ):
            n_residues = len(tx_grp)
            n_null = int(tx_grp["_is_null"].sum())
            constraint = tx_grp[~tx_grp["_is_null"]]
            null_rows = tx_grp[tx_grp["_is_null"]]

            if len(constraint) == 0:
                region_lengths: List[int] = []
                region_ranges: List[Dict[str, int]] = []
                n_regions = 0
            else:
                region_agg = constraint.groupby("region_index")["residue_index"].agg(
                    ["min", "max", "size"]
                )
                region_lengths = region_agg["size"].tolist()
                region_ranges = [
                    {
                        "start": int(row["min"]),
                        "end": int(row["max"]),
                        "length": int(row["size"]),
                    }
                    for _, row in region_agg.iterrows()
                ]
                n_regions = len(region_ranges)

            null_range: Optional[Dict[str, int]] = None
            if len(null_rows) > 0:
                null_range = {
                    "start": int(null_rows["residue_index"].min()),
                    "end": int(null_rows["residue_index"].max()),
                    "length": int(len(null_rows)),
                }

            rec = {
                "transcript_id": transcript_id,
                "n_regions": n_regions,
                "region_lengths": region_lengths,
                "region_ranges": region_ranges,
                "null_range": null_range,
                "n_null": n_null,
                "n_residues": n_residues,
            }
            if uniprot_id not in out:
                out[uniprot_id] = {}
            if run_label not in out[uniprot_id]:
                out[uniprot_id][run_label] = []
            out[uniprot_id][run_label].append(rec)

    return out


def _region_ranges_text(t: Dict[str, Any]) -> str:
    """Format region_ranges and null_range for one transcript as readable text."""
    parts = []
    ranges = t.get("region_ranges") or []
    for i, rr in enumerate(ranges, 1):
        parts.append(
            f"Region {i}: residues {rr['start']}–{rr['end']} ({rr['length']} residues)"
        )
    nr = t.get("null_range")
    if nr:
        parts.append(
            f"Null (catch-all): residues {nr['start']}–{nr['end']} ({nr['length']} residues)"
        )
    elif t.get("n_null", 0) > 0:
        parts.append(
            f"Null: {t['n_null']} residues (span not contiguous or not stored)"
        )
    return "; ".join(parts) if parts else "—"


def _format_gene_detail_section(
    gene_id: str,
    category: str,
    reason: str,
    details: Dict[str, List[Dict[str, Any]]],
    max_runs_show: int = 12,
    runs_for_ranges: int = 3,
) -> List[str]:
    """Format one gene's detailed subsection: summary, region definitions (residue ranges), and explicit differences."""
    lines: List[str] = []
    lines.append(f"#### **{gene_id}** ({category})")
    lines.append("")
    lines.append(reason)
    lines.append("")
    if not details:
        lines.append("*(No per-run details available.)*")
        lines.append("")
        return lines

    run_labels = sorted(details.keys())
    n_regions_per_run = [sum(t["n_regions"] for t in details[r]) for r in run_labels]
    mean_regions_per_tx = [
        sum(t["n_regions"] for t in details[r]) / max(1, len(details[r]))
        for r in run_labels
    ]
    null_frac_per_run = []
    for r in run_labels:
        total_res = sum(t["n_residues"] for t in details[r])
        total_null = sum(t["n_null"] for t in details[r])
        null_frac_per_run.append(100.0 * total_null / total_res if total_res else 0)

    lines.append(f"**Summary across {len(run_labels)} runs:**")
    lines.append("")
    lines.append(
        f"- Number of constraint regions (total over transcripts) ranges from **{min(n_regions_per_run)}** to **{max(n_regions_per_run)}** across runs."
    )
    lines.append(
        f"- Mean constraint regions per transcript ranges from **{min(mean_regions_per_tx):.1f}** to **{max(mean_regions_per_tx):.1f}**."
    )
    lines.append(
        f"- Null fraction (residues in catch-all region) ranges from **{min(null_frac_per_run):.1f}%** to **{max(null_frac_per_run):.1f}%**."
    )
    lines.append("")

    # Pick runs that best show the spread: fewest regions, most regions, and
    # one with mid null frac
    idx_min_regions = min(range(len(run_labels)), key=lambda i: n_regions_per_run[i])
    idx_max_regions = max(range(len(run_labels)), key=lambda i: n_regions_per_run[i])
    runs_to_describe = [run_labels[idx_min_regions], run_labels[idx_max_regions]]
    if len(run_labels) > 2:
        mid_idx = len(run_labels) // 2
        by_null = sorted(range(len(run_labels)), key=lambda i: null_frac_per_run[i])
        idx_mid_null = by_null[mid_idx]
        r_mid = run_labels[idx_mid_null]
        if r_mid not in runs_to_describe:
            runs_to_describe.append(r_mid)
    runs_to_describe = runs_to_describe[:runs_for_ranges]

    lines.append("**Constraint regions (residue ranges) — contrasting runs:**")
    lines.append("")
    lines.append(
        "Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region."
    )
    lines.append("")
    for r in runs_to_describe:
        tx_list = details[r]
        total_regions = sum(t["n_regions"] for t in tx_list)
        total_null = sum(t["n_null"] for t in tx_list)
        total_res = sum(t["n_residues"] for t in tx_list)
        null_pct = 100.0 * total_null / total_res if total_res else 0
        lines.append(
            f"- **Run: {r}** — {total_regions} constraint region(s) total, {null_pct:.1f}% in null."
        )
        for t in tx_list:
            lines.append(
                f"  - **Transcript {t['transcript_id']}**: {_region_ranges_text(t)}"
            )
        lines.append("")

    lines.append("**What is different across methods:**")
    lines.append("")
    min_regions_run = run_labels[n_regions_per_run.index(min(n_regions_per_run))]
    max_regions_run = run_labels[n_regions_per_run.index(max(n_regions_per_run))]
    min_null_run = run_labels[null_frac_per_run.index(min(null_frac_per_run))]
    max_null_run = run_labels[null_frac_per_run.index(max(null_frac_per_run))]
    lines.append(
        f"- **Region count:** **{min_regions_run}** has the fewest constraint regions ({min(n_regions_per_run)}); "
        f"**{max_regions_run}** has the most ({max(n_regions_per_run)}). "
        f"Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null."
    )
    lines.append(
        f"- **Null fraction:** Lowest in **{min_null_run}** ({min(null_frac_per_run):.1f}%); highest in **{max_null_run}** ({max(null_frac_per_run):.1f}%). "
        f"Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all."
    )
    all_null_runs = [
        r for r in run_labels if sum(t["n_regions"] for t in details[r]) == 0
    ]
    if all_null_runs:
        example = all_null_runs[0][:80] + ("..." if len(all_null_runs[0]) > 80 else "")
        lines.append(
            f"- **All-null runs:** In **{len(all_null_runs)}** run(s) this gene has *no* constraint regions (all residues in null). Example: {example}"
        )
    lines.append("")
    return lines


# ---------- Plots ----------
def write_plots(
    results: List[Tuple[str, Dict[str, Any]]],
    plot_dir: str,
    tsv_path_or_df: Any = None,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]], List[Tuple[str, str]]]:
    """Generate summary plots, distribution plots, and optional per-gene heatmaps.
    Returns (summary_plot_paths, distribution_plot_paths, per_gene_plot_paths).
    If tsv_path_or_df is provided, per-gene heatmaps (genes x runs) are generated.
    Requires matplotlib. Returns ([], [], []) if matplotlib not available or results have errors.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return [], [], []

    valid = [(l, s) for l, s in results if "error" not in s]
    if not valid:
        return [], [], []

    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    summary_paths: List[Tuple[str, str]] = []
    distribution_paths: List[Tuple[str, str]] = []
    per_gene_paths: List[Tuple[str, str]] = []

    # ----- Distribution plots (full distributions: violin plots) -----
    dist_data = [
        (l, s.get("_distributions")) for l, s in valid if s.get("_distributions")
    ]
    if dist_data:
        n_runs = len(dist_data)
        d_labels = [l for l, _ in dist_data]
        width = 6  # x-axis: metric (distribution)
        height = max(6, n_runs * 0.35)  # y-axis: run labels
        y_pos = np.arange(n_runs)
        violin_width = min(0.8, max(0.2, 14 / n_runs))

        def _violin_plot(
            ax, data_per_run, positions, labels, xlabel, title, color="steelblue"
        ):
            # Horizontal violins: y = run labels, x = metric (vert=False).
            parts = ax.violinplot(
                data_per_run,
                positions=positions,
                widths=violin_width * 0.9,
                showmeans=False,
                showmedians=False,
                showextrema=False,
                vert=False,
            )
            for pc in parts["bodies"]:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
            ax.set_yticks(positions)
            ax.set_yticklabels(labels, fontsize=6)
            ax.set_xlabel(xlabel)
            ax.set_title(title)

        # 1. Distribution: constraint regions per transcript (violin per run)
        fig, ax = plt.subplots(figsize=(width, height))
        _violin_plot(
            ax,
            [d["n_regions_per_tx"] for _, d in dist_data],
            y_pos,
            d_labels,
            "Constraint regions per transcript",
            "Distribution of constraint regions per transcript",
            color="steelblue",
        )
        plt.tight_layout()
        p = str(Path(plot_dir) / "dist_regions_per_tx.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        plt.close()
        distribution_paths.append(
            (
                "dist_regions_per_tx.png",
                "Distribution: constraint regions per transcript",
            )
        )

        # 2. Distribution: region length (residues per constraint region)
        dist_with_regions = [
            (l, d) for l, d in dist_data if len(d["region_lengths"]) > 0
        ]
        if dist_with_regions:
            n_runs_r = len(dist_with_regions)
            width_r = 6
            height_r = max(6, n_runs_r * 0.35)
            y_pos_r = np.arange(n_runs_r)
            violin_width_r = min(0.8, max(0.2, 14 / n_runs_r))
            fig, ax = plt.subplots(figsize=(width_r, height_r))
            parts = ax.violinplot(
                [d["region_lengths"] for _, d in dist_with_regions],
                positions=y_pos_r,
                widths=violin_width_r * 0.9,
                showmeans=False,
                showmedians=False,
                showextrema=False,
                vert=False,
            )
            for pc in parts["bodies"]:
                pc.set_facecolor("seagreen")
                pc.set_alpha(0.7)
            ax.set_yticks(y_pos_r)
            ax.set_yticklabels([l for l, _ in dist_with_regions], fontsize=6)
            ax.set_xlabel("Residues per constraint region")
            ax.set_title("Distribution of constraint region length")
            plt.tight_layout()
            p = str(Path(plot_dir) / "dist_region_length.png")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            plt.close()
            distribution_paths.append(
                (
                    "dist_region_length.png",
                    "Distribution: region length (residues per region)",
                )
            )

        # 3. Distribution: null fraction per transcript (violin per run)
        fig, ax = plt.subplots(figsize=(width, height))
        _violin_plot(
            ax,
            [d["null_frac_per_tx"] * 100 for _, d in dist_data],
            y_pos,
            d_labels,
            "Null fraction (%)",
            "Distribution of null fraction per transcript",
            color="coral",
        )
        plt.tight_layout()
        p = str(Path(plot_dir) / "dist_null_frac_per_tx.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        plt.close()
        distribution_paths.append(
            ("dist_null_frac_per_tx.png", "Distribution: null fraction per transcript")
        )

    # ----- Per-gene heatmaps (genes x runs) -----
    if tsv_path_or_df is not None:
        try:
            matrix_regions, matrix_null_frac, run_labels, gene_ids = (
                compute_per_gene_run_matrices(tsv_path_or_df)
            )
            n_genes = len(gene_ids)
            n_runs_h = len(run_labels)
            if n_genes > 0 and n_runs_h > 0:
                # Row height and col width for readability (cap figure size)
                row_h = max(0.08, min(0.25, 400 / n_genes))
                col_w = max(0.15, min(0.5, 300 / n_runs_h))
                fig_h = min(50, max(8, n_genes * row_h))
                fig_w = min(40, max(8, n_runs_h * col_w))
                figsize = (fig_w, fig_h)

                # Heatmap 1: mean constraint regions (hierarchical clustering on genes
                # and runs)
                _plot_clustered_heatmap(
                    matrix_regions,
                    gene_ids,
                    run_labels,
                    str(Path(plot_dir) / "heatmap_regions_per_gene.png"),
                    title="Mean constraint regions per transcript (per gene, per run)",
                    cbar_label="Mean regions",
                    cmap="viridis",
                    figsize=figsize,
                )
                per_gene_paths.append(
                    (
                        "heatmap_regions_per_gene.png",
                        "Per-gene heatmap: mean constraint regions per transcript (genes × runs, hierarchically clustered)",
                    )
                )

                # Heatmap 2: mean null fraction (hierarchical clustering)
                _plot_clustered_heatmap(
                    matrix_null_frac * 100,
                    gene_ids,
                    run_labels,
                    str(Path(plot_dir) / "heatmap_null_frac_per_gene.png"),
                    title="Mean null fraction % (per gene, per run)",
                    cbar_label="Null %",
                    cmap="Oranges",
                    vmin=0,
                    vmax=100,
                    figsize=figsize,
                )
                per_gene_paths.append(
                    (
                        "heatmap_null_frac_per_gene.png",
                        "Per-gene heatmap: mean null fraction % (genes × runs, hierarchically clustered)",
                    )
                )
        except Exception as e:
            import warnings

            warnings.warn(f"Per-gene heatmaps skipped: {e}", UserWarning, stacklevel=1)

    return summary_paths, distribution_paths, per_gene_paths


# ---------- Report writing ----------
def write_markdown_report(
    results: List[Tuple[str, Dict[str, Any]]],
    output_path: str,
    title: str = "Proemis3D forward method comparison (158 test genes)",
    results_one_per_gene: Optional[List[Tuple[str, Dict[str, Any]]]] = None,
    summary_plot_paths: Optional[List[Tuple[str, str]]] = None,
    distribution_plot_paths: Optional[List[Tuple[str, str]]] = None,
    per_gene_plot_paths: Optional[List[Tuple[str, str]]] = None,
    gene_examples: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    gene_example_details: Optional[Dict[str, Dict[str, List[Dict[str, Any]]]]] = None,
) -> None:
    """Write a Markdown report with tables and short summary. Optionally embed plots and gene examples (with full region details if gene_example_details provided). If results_one_per_gene is provided, a second summary table (one protein/transcript per gene) is included."""
    # Sort by model, then pLDDT, then PAE for readable grouping
    _order_model = {"AIC": 0, "AIC weight": 1, "LRT": 2}
    _order_plddt = {"–": 0, "exclude": 1, "remove": 2}
    _order_pae = {
        "–": 0,
        "truncate_center": 1,
        "filter_center": 2,
        "filter_region": 3,
        "exclude_center": 4,
        "exclude_region": 5,
    }

    def _sort_key(item):
        _, stats = item
        if "error" in stats:
            return (99, 99, 99, "")
        m = stats.get("model_display") or "–"
        p = stats.get("plddt_display") or "–"
        a = stats.get("pae_display") or "–"
        return (
            _order_model.get(m, 99),
            _order_plddt.get(p, 99),
            _order_pae.get(a, 99),
            str(m) + str(p) + str(a),
        )

    sorted_results = sorted(results, key=_sort_key)

    lines = [
        f"# {title}",
        "",
        "## Summary table (all runs)",
        "",
        "| Model comparison method | pLDDT filter | PAE filter | % Assigned missing | "
        "Const. regions (mean / median / min / max) | Region length (mean / median / min / max) | "
        "Null frac (mean) | Tx all-null % |",
        "|-------------------------|--------------|------------|--------------------|"
        "-----------------------------------------|------------------------------------------|"
        "------------------|----------------|",
    ]

    for label, stats in sorted_results:
        if "error" in stats:
            lines.append(f"| — | — | — | ERROR: {stats['error']} | — | — | — | — |")
            continue
        model = stats.get("model_display") or "–"
        plddt = stats.get("plddt_display") or "–"
        pae = stats.get("pae_display") or "–"
        # When every residue has a row (NA for filtered), row count equals
        # baseline so use fraction with NA
        baseline_total = stats.get("baseline_total_residues")
        total_residues = stats.get("total_residues")
        use_na_frac = baseline_total is None or total_residues == baseline_total
        missing_pct = (
            stats.get("frac_residues_missing_obs", 0) * 100
            if use_na_frac
            else stats.get("frac_residues_filtered_out", 0) * 100
        )
        regions_str = (
            f"{stats['n_constraint_regions_mean']:.1f} / "
            f"{stats['n_constraint_regions_median']:.0f} / "
            f"{stats['n_constraint_regions_min']} / {stats['n_constraint_regions_max']}"
        )
        rl_mean, rl_med = stats.get("region_length_mean"), stats.get(
            "region_length_median"
        )
        rl_min, rl_max = stats.get("region_length_min"), stats.get("region_length_max")
        len_str = (
            f"{rl_mean:.1f} / {rl_med:.0f} / {rl_min} / {rl_max}"
            if rl_mean is not None and rl_med is not None
            else "—"
        )
        null_pct = stats["frac_transcripts_all_null"] * 100
        lines.append(
            f"| {model} | {plddt} | {pae} | {missing_pct:.2f}% | "
            f"{regions_str} | {len_str} | "
            f"{stats['frac_residues_in_null_mean']*100:.1f}% | {null_pct:.1f}% |"
        )

    # Second summary table: one transcript per gene (e.g. 144 genes)
    if results_one_per_gene:
        sorted_one_per_gene = sorted(results_one_per_gene, key=_sort_key)
        n_genes_one = None
        for _, st in sorted_one_per_gene:
            if "error" not in st and st.get("n_transcripts") is not None:
                n_genes_one = st["n_transcripts"]
                break
        subheading = f" (one transcript per gene"
        if n_genes_one is not None:
            subheading += f", {n_genes_one} genes"
        subheading += ")"
        lines.extend(
            [
                "",
                f"## Summary table{subheading}",
                "",
                "| Model comparison method | pLDDT filter | PAE filter | % Assigned missing | "
                "Const. regions (mean / median / min / max) | Region length (mean / median / min / max) | "
                "Null frac (mean) | Tx all-null % |",
                "|-------------------------|--------------|------------|--------------------|"
                "-----------------------------------------|------------------------------------------|"
                "------------------|----------------|",
            ]
        )
        for label, stats in sorted_one_per_gene:
            if "error" in stats:
                lines.append(f"| — | — | — | ERROR: {stats['error']} | — | — | — | — |")
                continue
            model = stats.get("model_display") or "–"
            plddt = stats.get("plddt_display") or "–"
            pae = stats.get("pae_display") or "–"
            baseline_total = stats.get("baseline_total_residues")
            total_residues = stats.get("total_residues")
            use_na_frac = baseline_total is None or total_residues == baseline_total
            missing_pct = (
                stats.get("frac_residues_missing_obs", 0) * 100
                if use_na_frac
                else stats.get("frac_residues_filtered_out", 0) * 100
            )
            regions_str = (
                f"{stats['n_constraint_regions_mean']:.1f} / "
                f"{stats['n_constraint_regions_median']:.0f} / "
                f"{stats['n_constraint_regions_min']} / {stats['n_constraint_regions_max']}"
            )
            rl_mean, rl_med = stats.get("region_length_mean"), stats.get(
                "region_length_median"
            )
            rl_min, rl_max = stats.get("region_length_min"), stats.get(
                "region_length_max"
            )
            len_str = (
                f"{rl_mean:.1f} / {rl_med:.0f} / {rl_min} / {rl_max}"
                if rl_mean is not None and rl_med is not None
                else "—"
            )
            null_pct = stats["frac_transcripts_all_null"] * 100
            lines.append(
                f"| {model} | {plddt} | {pae} | {missing_pct:.2f}% | "
                f"{regions_str} | {len_str} | "
                f"{stats['frac_residues_in_null_mean']*100:.1f}% | {null_pct:.1f}% |"
            )

    lines.extend(
        [
            "",
            "## Metrics explained",
            "",
            "- **% Assigned missing**: When outputs have one row per residue (with NA for filtered): fraction of residues with NA for obs/exp in this run. When only assigned residues are emitted: fraction of residues absent from this run vs the unfiltered run (same model).",
            "- **Const. regions**: Number of constraint regions (non-null) per transcript: mean, median, min, max.",
            "- **Region length**: Residues per constraint region: mean, median, min, max.",
            "- **Null frac (mean)**: Mean over transcripts of (residues in null region / total residues).",
            "- **Tx all-null %**: Fraction of transcripts with zero constraint regions (all residues in null).",
            "",
            "## Expectation for missing values (NA) by filtering method",
            "",
            "**% Assigned missing** is the fraction of residues with **NA** for obs/exp in the per-residue output.",
            "",
            "The forward algorithm builds a set of **valid_residues**: the union of all residues that appear in at least one candidate region. Only valid_residues get a region assignment (either a constraint region or the null catch-all). The final output left-joins all residues (from the OE array) against the forward results, so residues **not** in valid_residues get **NA** for region_index, obs, exp, oe, etc.",
            "",
            "Whether a residue ends up in valid_residues depends on whether it was **hard-filtered** (physically removed from the distance matrix) or **soft-excluded** (kept in the distance matrix but marked `exclude_from_stats`).",
            "",
            "### pLDDT filter",
            "",
            "- **–** (none): All residues are in valid_residues. **0% NA**.",
            "- **exclude** (`exclude_low_plddt_from_stats`): Low-pLDDT residues stay in the distance matrix and appear in candidate regions → they are **in valid_residues** and get a region assignment (constraint or null). Their obs/exp are excluded from region-level nLL/OE calculations via `excluded_residues`, but the per-residue export row is **not NA**. **0% NA** from pLDDT exclude alone.",
            "- **remove** (`remove_low_plddt_residues`): Low-pLDDT residues are physically removed from the distance matrix → they never appear in any candidate region → **not in valid_residues** → **NA** in the export. Typically **~40–47% NA**.",
            "",
            "### PAE filter",
            "",
            "- **–** (none): No residues removed by PAE. **0% NA** from PAE.",
            "- **truncate_center**: All residues **after** the first with PAE(center, neighbor) &gt; cutoff are removed from the distance matrix (sequential cutoff). Residues that are truncated from **every** center's candidate region are not in valid_residues → **NA**.",
            "- **filter_center**: Only residues with PAE(center, neighbor) &gt; cutoff are removed from the distance matrix. Residues filtered from **every** center's candidate region are not in valid_residues → **NA**.",
            "- **filter_region**: Residues whose maximum pairwise PAE to any residue in the region &gt; cutoff are removed from the distance matrix. Residues filtered from **every** center are not in valid_residues → **NA**.",
            "- **exclude_center** (`exclude_on_pairwise_pae_with_center`): High-PAE residues stay in the distance matrix and appear in candidate regions → **in valid_residues**. They are marked `exclude_from_stats` (excluded from region nLL/OE), but their export row is **not NA**. Very few residues end up NA (only those that happen to not appear in any candidate region for other reasons).",
            "- **exclude_region** (`exclude_on_pairwise_pae_in_region`): Same as exclude_center but based on pairwise PAE to any residue in the region. Residues stay in valid_residues → export row is **not NA**. Slightly more residues may be excluded from stats than exclude_center, but still very low NA.",
            "",
            "### Combined",
            "",
            "Runs that use both pLDDT and PAE filters can have more NA residues. For example, **remove** + **truncate_center** removes low-pLDDT residues and truncates high-PAE residues, so both sets are not in valid_residues → NA. Conversely, **exclude** + **exclude_center** keeps all residues in valid_residues (0% NA from either) but excludes them from region statistics.",
            "",
            "The **% Assigned missing** column summarizes the fraction of residues with NA in each run.",
            "",
            "## What to look for when choosing a method",
            "",
            "- **More constraint regions** (higher mean/median) → more granular; **fewer** → more conservative.",
            "- **Larger null fraction** → more residues in catch-all; smaller → more residues assigned to constraint regions.",
            "- **High Tx all-null %** → many genes end up with no constraint regions (may be too strict).",
            "- **% Assigned missing** → residues with NA (or absent when only assigned residues are emitted); baseline = unfiltered run for that model.",
            "- **Region length**: very short regions may be noise; very long may be overmerged.",
            "- **LRT** is usually most conservative (fewest regions); **AIC** most permissive; **AIC weight** tunable.",
            "",
        ]
    )

    if gene_examples:
        lines.append("## Example genes to characterize")
        lines.append("")
        lines.append(
            "Genes below are candidates for full characterization: one stable across methods, and several outliers for different reasons. For each gene we show summary stats, **constraint regions as residue ranges** (for contrasting runs), and **what differs** across methods."
        )
        lines.append("")
        details_map = gene_example_details or {}
        if gene_examples.get("stable"):
            lines.append("### Stable (low variation in region count across methods)")
            lines.append("")
            for g in gene_examples["stable"]:
                det = details_map.get(g["gene"])
                if det:
                    lines.extend(
                        _format_gene_detail_section(
                            g["gene"],
                            "stable",
                            f"Low variation: std(mean regions across runs) = {g['std_mean_regions']:.2f}; rarely all-null.",
                            det,
                        )
                    )
                else:
                    lines.append(
                        f"- **{g['gene']}** — std(mean regions across runs) = {g['std_mean_regions']:.2f}; rarely all-null."
                    )
            lines.append("")
        if gene_examples.get("outlier_region_variance"):
            lines.append("### Outlier: region count varies most across methods")
            lines.append("")
            for g in gene_examples["outlier_region_variance"]:
                det = details_map.get(g["gene"])
                if det:
                    lines.extend(
                        _format_gene_detail_section(
                            g["gene"],
                            "region count variance",
                            f"Region count varies a lot: std(mean regions) = {g['std_mean_regions']:.2f}, median = {g['mean_regions_median']:.1f}.",
                            det,
                        )
                    )
                else:
                    lines.append(
                        f"- **{g['gene']}** — std(mean regions) = {g['std_mean_regions']:.2f}, median mean regions = {g['mean_regions_median']:.1f}."
                    )
            lines.append("")
        if gene_examples.get("outlier_often_all_null"):
            lines.append(
                "### Outlier: often all-null (no constraint regions) across runs"
            )
            lines.append("")
            for g in gene_examples["outlier_often_all_null"]:
                pct = g["frac_runs_any_all_null"] * 100
                det = details_map.get(g["gene"])
                if det:
                    lines.extend(
                        _format_gene_detail_section(
                            g["gene"],
                            "often all-null",
                            f"In {pct:.0f}% of runs this gene has no constraint regions (all residues in null).",
                            det,
                        )
                    )
                else:
                    lines.append(f"- **{g['gene']}** — all-null in {pct:.0f}% of runs.")
            lines.append("")
        if gene_examples.get("outlier_null_frac_range"):
            lines.append("### Outlier: null fraction varies most across runs")
            lines.append("")
            for g in gene_examples["outlier_null_frac_range"]:
                det = details_map.get(g["gene"])
                if det:
                    lines.extend(
                        _format_gene_detail_section(
                            g["gene"],
                            "null fraction range",
                            f"Null fraction across runs has range {g['range_null_frac']:.2f} (full spread).",
                            det,
                        )
                    )
                else:
                    lines.append(
                        f"- **{g['gene']}** — range of mean null fraction across runs = {g['range_null_frac']:.2f}."
                    )
            lines.append("")
        if gene_examples.get("user_specified"):
            lines.append("### User-specified genes")
            lines.append("")
            for g in gene_examples["user_specified"]:
                det = details_map.get(g["gene"])
                if det:
                    lines.extend(
                        _format_gene_detail_section(
                            g["gene"],
                            "user-specified",
                            "User-requested gene for characterization.",
                            det,
                        )
                    )
                else:
                    lines.append(
                        f"- **{g['gene']}** — user-specified (no per-run details available)."
                    )
            lines.append("")

    if distribution_plot_paths or per_gene_plot_paths:
        lines.append("## Plots")
        lines.append("")
    if distribution_plot_paths:
        lines.append(
            "Violin plots show the full distribution of each metric across transcripts or regions (one violin per run)."
        )
        lines.append("")
        for rel_path, caption in distribution_plot_paths:
            lines.append(f"### {caption}")
            lines.append("")
            lines.append(f"![{caption}]({rel_path})")
            lines.append("")
    if per_gene_plot_paths:
        lines.append(
            "Per-gene heatmaps show each metric (mean constraint regions or null fraction) for every gene (rows) and run (columns)."
        )
        lines.append("")
        for rel_path, caption in per_gene_plot_paths:
            lines.append(f"### {caption}")
            lines.append("")
            lines.append(f"![{caption}]({rel_path})")
            lines.append("")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def write_csv_summary(
    results: List[Tuple[str, Dict[str, Any]]], output_path: str
) -> None:
    """Write a flat CSV of key metrics for each run (for plotting in R/Python)."""
    import csv

    rows = []
    for label, stats in results:
        if "error" in stats:
            rows.append({"label": label, "error": stats["error"]})
            continue
        row = {
            "label": label,
            **{
                k: ("" if v is None else v)
                for k, v in stats.items()
                if not isinstance(v, (dict, list, np.ndarray))
            },
        }
        rows.append(row)

    if not rows:
        return
    keys = list(rows[0].keys())
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Compare Proemis3D forward outputs: export HTs to one TSV, then report from TSV (pandas)."
    )
    parser.add_argument(
        "--config",
        type=str,
        help="Path to YAML/JSON config with 'runs': [{path, label}, ...]",
    )
    parser.add_argument(
        "--paths",
        nargs="+",
        metavar="PATH:LABEL",
        help="Alternatives to config: list of 'path.ht:Label' entries",
    )
    parser.add_argument(
        "--output-tsv",
        type=str,
        metavar="PATH",
        help="Export combined HT (union of all runs with run_label + method columns) to this TSV. Requires --config or --paths.",
    )
    parser.add_argument(
        "--tsv",
        type=str,
        metavar="PATH",
        help="Read from this TSV to generate report (no Hail). Use when TSV already exists.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="forward_methods_report.md",
        help="Output Markdown report path",
    )
    parser.add_argument(
        "--csv",
        type=str,
        default=None,
        help="Optional CSV path for summary metrics",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Proemis3D forward method comparison (158 test genes)",
        help="Report title",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Limit to first N runs when building TSV (for quick testing)",
    )
    parser.add_argument(
        "--backend",
        type=str,
        choices=("spark", "batch"),
        default="spark",
        help="Hail backend: 'spark' for Dataproc/local Spark (default), 'batch' for Hail Batch Service",
    )
    parser.add_argument(
        "--plot-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="Generate plots and save to DIR; embed in report (requires matplotlib). Use e.g. 'report_plots' next to the report.",
    )
    parser.add_argument(
        "--gene-examples",
        type=int,
        default=0,
        metavar="N",
        help="Add a section with example genes to characterize: 1 stable (low variation) and N outliers per category (region variance, often all-null, null frac range). 0 = off.",
    )
    parser.add_argument(
        "--extra-genes",
        nargs="+",
        metavar="UNIPROT_ID",
        help="Extra UniProt IDs to include in the 'Example genes to characterize' section (e.g. P13637 for ATP1A3). These are added under a 'User-specified' heading with the same detail output as auto-selected genes. Implies --gene-examples 2 if --gene-examples is 0.",
    )
    parser.add_argument(
        "--gene-mapping",
        type=str,
        default=None,
        metavar="TSV",
        help="Optional TSV with columns transcript_id and gene_id. When provided, merged with the main TSV so the 'one transcript per gene' table groups by gene_id (e.g. 144 genes). Without this, grouping uses uniprot_id (one per protein).",
    )
    args = parser.parse_args()

    # Determine TSV path for report
    tsv_path = args.tsv
    if tsv_path is None and args.output_tsv:
        tsv_path = args.output_tsv

    if tsv_path is None:
        print(
            "Provide either --tsv (report from existing TSV) or --output-tsv (build TSV from config/paths).",
            file=sys.stderr,
        )
        sys.exit(1)

    # Build TSV from Hail if requested
    if args.output_tsv:
        if not args.config and not args.paths:
            print("For --output-tsv provide --config or --paths.", file=sys.stderr)
            sys.exit(1)
        if args.config:
            runs = load_config(args.config)
        else:
            runs = []
            for s in args.paths:
                if ":" in s:
                    path, label = s.rsplit(":", 1)
                    runs.append((path.strip(), label.strip()))
                else:
                    runs.append((s, Path(s).stem))
        if not runs:
            print("No runs to compare.", file=sys.stderr)
            sys.exit(1)
        export_union_ht_to_tsv(
            runs, args.output_tsv, limit=args.limit, backend=args.backend
        )

    # Report from TSV (pandas); read once for stats and optional gene examples
    # / one-per-gene table
    df = _read_and_normalize_tsv(tsv_path)
    if args.gene_mapping:
        df = _apply_gene_mapping(df, args.gene_mapping)
    results = compute_stats_from_tsv(df)
    results_one_per_gene = compute_stats_from_tsv(
        restrict_to_one_transcript_per_gene(df)
    )
    # If --extra-genes is provided, ensure gene examples section is enabled.
    effective_gene_examples = args.gene_examples or 0
    if args.extra_genes and effective_gene_examples == 0:
        effective_gene_examples = 2

    if effective_gene_examples > 0:
        gene_examples = compute_gene_examples(
            df, n_stable=1, n_outliers=min(effective_gene_examples, 5)
        )
        # Add user-specified extra genes under 'user_specified' category.
        if args.extra_genes:
            available_uniprots = set(df["uniprot_id"].dropna().unique())
            user_specified = []
            for uid in args.extra_genes:
                if uid in available_uniprots:
                    user_specified.append({"gene": uid})
                else:
                    print(
                        f"Warning: UniProt ID '{uid}' not found in TSV data, skipping.",
                        file=sys.stderr,
                    )
            if user_specified:
                gene_examples["user_specified"] = user_specified

        all_gene_ids = []
        for key in (
            "stable",
            "outlier_region_variance",
            "outlier_often_all_null",
            "outlier_null_frac_range",
            "user_specified",
        ):
            for g in gene_examples.get(key) or []:
                all_gene_ids.append(g["gene"])
        gene_example_details = (
            get_gene_example_details(df, all_gene_ids) if all_gene_ids else None
        )
    else:
        gene_examples = None
        gene_example_details = None

    summary_plot_paths: List[Tuple[str, str]] = []
    distribution_plot_paths: List[Tuple[str, str]] = []
    per_gene_plot_paths: List[Tuple[str, str]] = []
    if args.plot_dir:
        summary_files, distribution_files, per_gene_files = write_plots(
            results, args.plot_dir, tsv_path_or_df=tsv_path
        )
        report_dir = Path(args.output).resolve().parent
        plot_dir_resolved = Path(args.plot_dir).resolve()
        try:
            rel_dir = plot_dir_resolved.relative_to(report_dir)
            summary_plot_paths = [(str(rel_dir / f), cap) for f, cap in summary_files]
            distribution_plot_paths = [
                (str(rel_dir / f), cap) for f, cap in distribution_files
            ]
            per_gene_plot_paths = [(str(rel_dir / f), cap) for f, cap in per_gene_files]
        except ValueError:
            summary_plot_paths = [
                (f"{args.plot_dir}/{f}", cap) for f, cap in summary_files
            ]
            distribution_plot_paths = [
                (f"{args.plot_dir}/{f}", cap) for f, cap in distribution_files
            ]
            per_gene_plot_paths = [
                (f"{args.plot_dir}/{f}", cap) for f, cap in per_gene_files
            ]
        n_plots = len(summary_files) + len(distribution_files) + len(per_gene_files)
        if n_plots:
            print(f"Wrote {n_plots} plots to {args.plot_dir}")
        else:
            print(
                "No plots generated (install matplotlib for --plot-dir)",
                file=sys.stderr,
            )
    write_markdown_report(
        results,
        args.output,
        title=args.title,
        results_one_per_gene=results_one_per_gene,
        summary_plot_paths=summary_plot_paths or None,
        distribution_plot_paths=distribution_plot_paths or None,
        per_gene_plot_paths=per_gene_plot_paths or None,
        gene_examples=gene_examples,
        gene_example_details=gene_example_details,
    )
    print(f"Wrote report to {args.output}")

    if args.csv:
        write_csv_summary(results, args.csv)
        print(f"Wrote CSV to {args.csv}")


if __name__ == "__main__":
    main()
