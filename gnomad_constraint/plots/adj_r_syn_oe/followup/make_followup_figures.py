"""Draw every figure in rr_syn_dnm_followup.md.

Some performance values are copied from Ruchit Panchal's rr_syn_dnm_report.v2.html
(sections 1, 3, 5, 8, 9 and 10). The other panels read the tables written by the
analysis scripts in this directory (remap_rr_to_shet.py, rr_shrinkage_experiments.py,
residual_noise.py, outlier_drivers.py, decile_explained_removed.py,
rr_background_experiments.py, rr_logistic_fits.py, histone_checks.py and
histone_gene_conversion.py).
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

NONE_C = "#898781"
ADJ_R_C = "#1baf7a"
LEGACY_C = "#2a78d6"
RUCHIT_C = "#eb6834"
REC_C = "#2b2b2b"
"""Ridge-penalized rr: near-black so it stands out from every series hue."""
ROULETTE_C = "#6250d6"
"""Per-variant rates learned from gnomAD rare variants (Roulette and home-built)."""

SAME_FILL, SAME_EDGE = "#f0efec", "#c3c2b7"
SYN_FILL, SYN_EDGE = "#cde2fb", "#2a78d6"
NEW_FILL, NEW_EDGE = "#fbe1d6", "#eb6834"
RR_FILL, RR_EDGE = "#fdf0cc", "#eda100"
REC_FILL, REC_EDGE = "#ffffff", REC_C
RATE_FILL, RATE_EDGE = "#e7e3fa", ROULETTE_C

CLASS_C = {"A/T": "#1a9e6e", "C non-CpG": "#8f84e0", "CpG": "#d9622b"}

THEORY = dict(color=INK2, lw=1.3, ls=(0, (1, 2)))
CHANCE_FLOOR = 0.70
"""Upper bound on mean z² from sampling alone for the core set: each site is
observed or not, so a transcript's variance is at most sum(n p (1 - p)) per
(window, context) cell, not its expected count."""

plt.rcParams.update(
    {
        "font.family": ["Arial", "DejaVu Sans"],
        "font.size": 15,
        "axes.labelsize": 15,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 13.5,
        "text.color": INK,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": INK2,
        "axes.titlesize": 16,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    }
)


def box(ax, x, y, w, h, fill, edge, lines, size=8.8, dash=False, center=False, lw=1.3):
    """Draw a rounded box with a bold first line and regular following lines."""
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0,rounding_size=6",
            facecolor=fill,
            edgecolor=edge,
            linewidth=lw,
            linestyle=(0, (4, 3)) if dash else "-",
        )
    )
    step = size * 1.75
    y0 = y + h / 2 - step * (len(lines) - 1) / 2
    for i, line in enumerate(lines):
        ax.text(
            x + (w / 2 if center else 10),
            y0 + i * step,
            line,
            ha="center" if center else "left",
            va="center",
            fontsize=size if i == 0 else size - 0.4,
            fontweight="bold" if i == 0 else "normal",
            color=INK if i == 0 else INK2,
        )


def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=1.2):
    """Draw a straight arrow."""
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=0, shrinkB=0),
    )


def draw_diagram():
    """Draw the shared pipeline and how each rr version builds its table."""
    W, H = 1820, 1395
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    # A. Shared pipeline.
    ax.text(
        10,
        22,
        "A. Constraint pipeline: identical in every version",
        fontsize=15,
        fontweight="bold",
        va="center",
    )
    ax.text(
        10,
        48,
        "No correction = steps 1 → 4 → 6. Every rr version adds step 5 and changes"
        " nothing else.",
        fontsize=12.5,
        color=INK2,
        va="center",
    )
    top = [
        ["1  Context table", "all possible SNVs;", "context, methylation,", "coverage"],
        [
            "2  Mutation rate μ",
            "per context and",
            "methylation level,",
            "gnomAD genomes",
        ],
        ["3  Plateau + coverage", "models, calibrated on", "exome synonymous", "sites"],
        ["4  Per-variant", "expected count"],
        [
            "5  Expected × rr",
            "per variant, before summing;",
            "rate versions also × a",
            "per-variant multiplier",
        ],
        ["6  Sum per transcript", "obs/exp, LOEUF,", "z-scores, flags"],
    ]
    bw, gap, by, bh = 278, 30, 72, 112
    for i, lines in enumerate(top):
        x = 10 + i * (bw + gap)
        is_rr = i == 4
        box(
            ax,
            x,
            by,
            bw,
            bh,
            RR_FILL if is_rr else SAME_FILL,
            RR_EDGE if is_rr else SAME_EDGE,
            lines,
            size=12,
            dash=is_rr,
        )
        if i < len(top) - 1:
            arrow(ax, x + bw + 2, by + bh / 2, x + bw + gap - 2, by + bh / 2)

    # B. rr recipes, three lanes.
    lane_x0, lane_w, lane_gap = 190, 316, 10
    lane_x = [lane_x0 + i * (lane_w + lane_gap) for i in range(5)]
    lane_cx = [x + lane_w / 2 for x in lane_x]
    rr_cx = 10 + 4 * (bw + gap) + bw / 2
    yb = 214
    ax.text(10, yb, "B. How the rr table", fontsize=15, fontweight="bold", va="center")
    ax.text(
        10,
        yb + 24,
        "for step 5 is built",
        fontsize=15,
        fontweight="bold",
        va="center",
    )
    ax.plot([lane_cx[0], lane_cx[-1]], [yb, yb], color=RR_EDGE, lw=1.4)
    for cx in lane_cx:
        ax.plot([cx, cx], [yb, yb + 14], color=RR_EDGE, lw=1.4)
    arrow(ax, rr_cx, yb, rr_cx, by + bh + 3, color=RR_EDGE, lw=1.6)

    headers = [
        (ADJ_R_C, "Gnocchi adj_r", "non-coding constraint (Chen et al. 2023)"),
        (LEGACY_C, "Legacy syn_adj_r", "Ruchit's notebook refit (2025)"),
        (RUCHIT_C, "Ruchit's pipeline now", "branch rp/add_mu_adjustment"),
        (REC_C, "Ridge fit", "Ruchit's pipeline with a ridge penalty"),
        (ROULETTE_C, "Per-variant rates + ridge", "a quick first try, untuned"),
    ]
    hy, hh = 232, 60
    for x, (c, name, sub) in zip(lane_x, headers):
        ax.add_patch(
            FancyBboxPatch(
                (x, hy),
                lane_w,
                hh,
                boxstyle="round,pad=0,rounding_size=6",
                facecolor="white",
                edgecolor=c,
                linewidth=1.8,
            )
        )
        ax.scatter([x + 18], [hy + 19], s=90, color=c, zorder=3)
        ax.text(x + 32, hy + 19, name, fontsize=13.5, fontweight="bold", va="center")
        ax.text(x + 12, hy + 43, sub[:44], fontsize=11.5, color=INK2, va="center")

    S, B, N, R, T = "same", "syn", "new", "rec", "rate"
    rows = [
        (
            "Mutated sites\n(cases)",
            [
                (
                    S,
                    [
                        "WGS trio DNMs",
                        "deCODE + PsychENCODE;",
                        "genome-wide, mostly non-coding",
                    ],
                ),
                (
                    B,
                    [
                        "Synonymous DNMs, 7 cohorts",
                        "17,533 rows (duplicates kept);",
                        "chrX dropped",
                    ],
                ),
                (
                    N,
                    [
                        "Synonymous DNMs, deduplicated",
                        "17k set: 17,059 (202 on chrX), or",
                        "Satterstrom 2026 set: 23,336",
                    ],
                ),
                (N, ["Satterstrom 2026 DNMs", "23,336; no union with the 17k set"]),
                (N, ["Same as the ridge fit"]),
            ],
        ),
        (
            "Background\n(controls)",
            [
                (S, ["Non-mutated genome sites", "about 10× the DNM count"]),
                (
                    B,
                    [
                        "Synonymous sites",
                        "not a DNM, absent from gnomAD v4.1",
                        "exomes; 5% sample (459,612)",
                    ],
                ),
                (B, ["Same design", "one site per locus (479,021)"]),
                (B, ["Same as Ruchit's"]),
                (B, ["Same as Ruchit's"]),
            ],
        ),
        (
            "Regional\nfeatures",
            [
                (
                    S,
                    [
                        "13 features × 4 scales",
                        "1 kb to 1 Mb windows;",
                        "published matrix, chr1–22 only",
                    ],
                ),
                (S, ["Same"]),
                (N, ["Same + rebuilt chrX/Y", "7 of 13 features have X/Y data"]),
                (R, ["Same 13 features × 4 scales", "chrX/Y features not used"]),
                (R, ["Same as the ridge fit"]),
            ],
        ),
        (
            "Per-variant\nrates",
            [
                (S, ["None", "trinucleotide μ only"]),
                (S, ["Same: none"]),
                (S, ["Same: none"]),
                (S, ["Same: none"]),
                (
                    T,
                    [
                        "Roulette (Seplyarskiy 2023), or",
                        "home-built from gnomAD v4.1",
                        "genome rare variants; CpG = 1",
                    ],
                ),
            ],
        ),
        (
            "Feature\nselection",
            [
                (
                    S,
                    [
                        "Bonferroni screen",
                        "0.05/4/13 (CpG: 0.05/4/8);",
                        "about all 32 contexts modeled",
                    ],
                ),
                (
                    B,
                    [
                        "Same rule, fewer DNMs",
                        "12 of 32 contexts modeled;",
                        "the other 20 get rr = 1",
                    ],
                ),
                (
                    N,
                    [
                        "Tunable cutoff, Bonferroni default",
                        "13 of 32 contexts (17k DNMs);",
                        "16 of 32 (Satterstrom)",
                    ],
                ),
                (
                    R,
                    [
                        "No p-value screen",
                        "every feature goes in; the ridge",
                        "penalty decides how much each counts",
                    ],
                ),
                (R, ["Same as the ridge fit"]),
            ],
        ),
        (
            "Per-context\nmodel",
            [
                (
                    S,
                    [
                        "z-score → PCA → logistic",
                        "methylation-linked features",
                        "dropped for CpG contexts",
                    ],
                ),
                (S, ["Same"]),
                (S, ["Same", "exact port (matches legacy to 5e-13)"]),
                (
                    R,
                    [
                        "z-score → logistic + ridge",
                        "penalty chosen on held-out DNMs",
                        "(no gnomAD data used)",
                    ],
                ),
                (
                    T,
                    [
                        "Same, with log(rate) as a fixed",
                        "offset: rr learns only what",
                        "the rates miss",
                    ],
                ),
            ],
        ),
        (
            "rr scaling",
            [
                (
                    S,
                    [
                        "rr = P(window) / P(average window)",
                        "averages 1 over the genome",
                    ],
                ),
                (S, ["Same"]),
                (
                    N,
                    [
                        "Same, then ÷ exome-syn mean",
                        "per context (open question);",
                        "strength α tested, not in code",
                    ],
                ),
                (N, ["Normalized per context", "no α needed"]),
                (T, ["Normalized per context", "with the multiplier included"]),
            ],
        ),
        (
            "Applied to\nexpected",
            [
                (
                    S,
                    [
                        "× per-variant expected",
                        "your analysis scripts;",
                        "summed in PR #55",
                    ],
                ),
                (S, ["Same"]),
                (
                    N,
                    [
                        "In --aggregate-per-variant-expected",
                        "on by default;",
                        "unadjusted values kept as unadj_*",
                    ],
                ),
                (N, ["Same as Ruchit's"]),
                (T, ["Expected × multiplier × rr", "--variant-rate-ht-path"]),
            ],
        ),
        (
            "chrX / chrY",
            [
                (S, ["No rr", "feature matrix is chr1–22"]),
                (S, ["Same: no rr"]),
                (
                    N,
                    [
                        "Partial rr (open question)",
                        "unmasked: 3 of 13 contexts;",
                        "masked (flag): all, but worse",
                    ],
                ),
                (R, ["No rr", "the gain would be under 1%"]),
                (R, ["Same as the ridge fit"]),
            ],
        ),
        (
            "Code",
            [
                (S, ["atgu/gnomad_nc_constraint", "plus released model pickles"]),
                (B, ["Sandbox notebooks", "not in a repo"]),
                (N, ["regional_rate_pipeline.py", "in gnomad-constraint"]),
                (N, ["Ruchit's pipeline", "ridge in fit_rr_context_model"]),
                (
                    T,
                    [
                        "regional_rate_pipeline.py",
                        "--use-variant-rates;",
                        "build_homebuilt_rates.py",
                    ],
                ),
            ],
        ),
    ]
    styles = {
        S: (SAME_FILL, SAME_EDGE),
        B: (SYN_FILL, SYN_EDGE),
        N: (NEW_FILL, NEW_EDGE),
        R: (REC_FILL, REC_EDGE),
        T: (RATE_FILL, RATE_EDGE),
    }
    ry, rh, rgap = hy + hh + 16, 82, 12
    for r, (label, cells) in enumerate(rows):
        y = ry + r * (rh + rgap)
        ax.text(
            176,
            y + rh / 2,
            label,
            ha="right",
            va="center",
            fontsize=12.5,
            fontweight="bold",
            color=INK2,
            linespacing=1.3,
        )
        for x, (kind, lines) in zip(lane_x, cells):
            fill, edge = styles[kind]
            lw = 2.2 if kind == R else 1.3
            if len(lines) == 1:
                box(
                    ax,
                    x,
                    y,
                    lane_w,
                    rh,
                    fill,
                    edge,
                    lines,
                    size=11.8,
                    center=True,
                    lw=lw,
                )
            else:
                box(ax, x, y, lane_w, rh, fill, edge, lines, size=11.8, lw=lw)
            if r < len(rows) - 1:
                arrow(ax, x + lane_w / 2, y + rh + 1, x + lane_w / 2, y + rh + rgap - 1)
        if r == 0:
            y_first = y
    y_end = ry + len(rows) * (rh + rgap)
    for x in lane_x:
        arrow(ax, x + lane_w / 2, hy + hh + 1, x + lane_w / 2, y_first - 1)

    # Legend.
    ly = y_end + 16
    items = [
        (SAME_FILL, SAME_EDGE, "Same as Gnocchi's original"),
        (SYN_FILL, SYN_EDGE, "Changed for syn_adj_r (Ruchit keeps it)"),
        (NEW_FILL, NEW_EDGE, "New in Ruchit's pipeline"),
        (REC_FILL, REC_EDGE, "Changed for the ridge fit"),
        (RATE_FILL, RATE_EDGE, "Added for the per-variant rates"),
        (RR_FILL, RR_EDGE, "Where rr enters"),
    ]
    for n, (fill, edge, text) in enumerate(items):
        lx = 190 + (n % 3) * 540
        ly_n = ly + (n // 3) * 34
        ax.add_patch(
            FancyBboxPatch(
                (lx, ly_n - 10),
                26,
                20,
                boxstyle="round,pad=0,rounding_size=3",
                facecolor=fill,
                edgecolor=edge,
                linewidth=1.3,
            )
        )
        ax.text(lx + 34, ly_n, text, fontsize=12.5, va="center", color=INK2)
    ax.text(
        190,
        ly + 76,
        "Option tried on several versions (the current method, the ridge fit): leave"
        " genes with a QC flag out of the DNMs and background (--exclude-qc-flagged-genes).",
        fontsize=12.5,
        va="center",
        color=INK2,
    )

    fig.savefig(
        os.path.join(OUT_DIR, "rr_pipelines_diagram.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def read_outliers():
    """Read synonymous outlier counts written by remap_rr_to_shet.py --outliers."""
    return pd.read_csv(
        os.path.join(OUT_DIR, "syn_outliers_by_flag.rr_eval_export.tsv"), sep="\t"
    )


def out_count(out, corr, direction, flag=None):
    """Count outliers of one correction and direction, optionally of one flag group."""
    d = out[(out.correction == corr) & (out.direction == direction)]
    if flag is not None:
        d = d[d.flag == flag]
    return int(d.n.sum())


def style_axis(ax, grid_axis="x"):
    """Apply recessive axes: hairline grid, no top/right spines."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


PERF_VERSIONS = [
    ("none", "No correction", NONE_C, None, "-"),
    ("adj_r", "Gnocchi adj_r", ADJ_R_C, None, "-"),
    ("syn_adj_r_norm", "Legacy syn_adj_r", LEGACY_C, None, "-"),
    (
        "bonferroni_full",
        "Ruchit's pipeline\n(current, Bonferroni)",
        RUCHIT_C,
        None,
        "-",
    ),
    (
        "1e-2_a0.625_full",
        "Ruchit's latest pick\n(1e-2, α 0.625)",
        RUCHIT_C,
        "///",
        (0, (4, 2.5)),
    ),
    ("ridge_qcexcl", "Ridge (recommended)", REC_C, None, "-"),
    ("roulette_ridge", "Roulette + ridge", ROULETTE_C, None, "-"),
    (
        "homebuilt_ridge",
        "Home-built rates\n+ ridge",
        ROULETTE_C,
        "///",
        (0, (4, 2.5)),
    ),
]
"""The versions compared on the cell-level scoring: key, label, colour, hatch, line."""


def read_cell_scores(table):
    """Read a cell-level score table of the full run and of the per-variant rates."""
    return pd.concat(
        pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.{table}.tsv"), sep="\t")
        for prefix in ("full_run", "rr_rates")
    ).drop_duplicates(["version"] + (["alpha"] if table == "strength" else []))


def draw_performance():
    """Draw the rr performance panels from the cell-level scores of every version."""
    met = read_cell_scores("metrics").set_index("version")
    strength = read_cell_scores("strength")
    out = pd.concat(
        pd.read_csv(os.path.join(OUT_DIR, f"{p}.outliers_by_flag.tsv"), sep="\t")
        for p in ("full_run", "rr_rates")
    ).drop_duplicates(["version", "direction", "category"])
    keys = [v[0] for v in PERF_VERSIONS]
    names = [v[1] for v in PERF_VERSIONS]
    style = {v[0]: v[2:] for v in PERF_VERSIONS}
    short = {
        "adj_r": "Gnocchi adj_r",
        "syn_adj_r_norm": "Legacy syn_adj_r",
        "bonferroni_full": "Ruchit's pipeline",
        "1e-2_a0.625_full": "Ruchit's latest pick",
        "ridge_qcexcl": "Ridge (recommended)",
        "roulette_ridge": "Roulette + ridge",
        "homebuilt_ridge": "Home-built + ridge",
    }
    names = [short.get(k, "No correction") for k in keys]
    ypos = list(range(len(keys)))[::-1]

    def curve(k):
        d = strength[strength.version == k].sort_values("alpha")
        return [0.0, *d.alpha], [met.loc["none", "core_noise"], *d.core_noise]

    fig, axes = plt.subplots(3, 2, figsize=(13, 16))
    fig.subplots_adjust(
        left=0.2, right=0.97, top=0.9, bottom=0.05, hspace=0.45, wspace=0.6
    )
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [],
            [],
            color=style[k][0],
            ls=style[k][2],
            lw=2.6,
            marker="D" if k == "ridge_qcexcl" else "o",
            ms=9,
            mfc=SURFACE if style[k][1] else style[k][0],
            mec=style[k][0],
            label=short.get(k, "No correction"),
        )
        for k in keys
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.995),
        ncol=4,
        frameon=False,
        fontsize=14,
        handlelength=2.6,
        columnspacing=1.6,
    )

    # (a) Per-gene noise.
    ax = axes[0, 0]
    for y, k in zip(ypos, keys):
        c, h, _ = style[k]
        v = met.loc[k, "core_noise"]
        ax.barh(y, v, height=0.6, color=c, hatch=h, hatchcolor="white", linewidth=0)
        ax.text(v + 0.06, y, f"{v:.2f}", va="center", fontsize=14, fontweight="bold")
    ax.axvline(met.loc["none", "core_noise"], color=MUTED, lw=1, ls=(0, (3, 3)))
    ax.axvline(CHANCE_FLOOR, zorder=0, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.08,
        len(keys) - 0.4,
        f"chance ≤ {CHANCE_FLOOR:.2f}",
        fontsize=13,
        color=INK2,
        va="center",
    )
    ax.set_yticks(ypos, names)
    ax.set_xlim(0, 6.2)
    ax.set_ylim(-0.6, len(keys) - 0.05)
    ax.set_title("a  Per-gene synonymous noise")
    ax.set_xlabel("mean z², overall level removed (lower = better)")
    style_axis(ax)

    # (b) Outliers in genes with and without a QC flag.
    ax = axes[0, 1]
    qc = {"low_exome_mapping_quality", "low_exome_coverage"}
    for y, k in zip(ypos, keys):
        c, h, _ = style[k]
        d = out[out.version == k]
        flagged = int(d[d.category.isin(qc)].n.sum())
        unflagged = int(d.n.sum()) - flagged
        ax.barh(y, flagged, height=0.6, color=AXIS, edgecolor=SURFACE, linewidth=2)
        ax.text(flagged / 2, y, f"{flagged}", ha="center", va="center", fontsize=14)
        ax.barh(
            y,
            unflagged,
            left=flagged,
            height=0.6,
            color=c,
            hatch=h,
            hatchcolor="white",
            edgecolor=SURFACE,
            linewidth=2,
        )
        ax.text(flagged + unflagged + 8, y, f"{unflagged}", va="center", fontsize=14)
    ax.set_yticks(ypos, [""] * len(keys))
    ax.set_xlim(0, 520)
    ax.set_ylim(-0.6, len(keys) - 0.05)
    ax.set_title("b  Synonymous outliers (|z| > 8)")
    ax.set_xlabel(
        "transcripts; gray = genes with a QC flag,\ncolour = genes without one"
    )
    style_axis(ax)

    # (c) Calibration slope.
    ax = axes[1, 0]
    skeys = keys[1:]
    sy = list(range(len(skeys)))[::-1]
    for y, k in zip(sy, skeys):
        c, h, _ = style[k]
        v = met.loc[k, "slope"]
        ax.barh(y, v, height=0.6, color=c, hatch=h, hatchcolor="white", linewidth=0)
        ax.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=14)
    ax.axvline(1.0, **THEORY)
    ax.text(
        0.98,
        len(skeys) - 0.45,
        "1 = right strength",
        ha="right",
        fontsize=13,
        color=INK2,
        va="bottom",
    )
    ax.set_yticks(sy, names[1:])
    ax.set_xlim(0, 1.25)
    ax.set_ylim(-0.6, len(skeys) + 0.3)
    ax.set_title("c  Calibration slope (below 1 = rr too strong)")
    ax.set_xlabel("log-log slope of uncorrected obs/exp on gene rr")
    style_axis(ax)

    # (d) Strength.
    ax = axes[1, 1]
    ax.axvline(1.0, color=AXIS, lw=1)
    a0, n0 = curve("ridge_qcexcl")
    best, start = min(n0), met.loc["none", "core_noise"]
    grid = np.linspace(0, 1.5, 151)
    ax.plot(
        grid,
        best + (start - best) * (grid - 1) ** 2,
        label="Perfectly scaled (theory)",
        **THEORY,
    )
    for k in keys[2:]:
        c, h, ls = style[k]
        a, n = curve(k)
        rec = k == "ridge_qcexcl"
        ax.plot(a, n, color=c, ls=ls, lw=2.8 if rec else 2, label=short[k])
        i = int(np.argmin(n))
        ax.scatter(
            [a[i]],
            [n[i]],
            s=74 if rec else 61,
            marker="D" if rec else "o",
            color=c,
            edgecolor=SURFACE,
            lw=1.5,
            zorder=4,
        )
        ax.scatter(
            [met.loc[k, "slope"]],
            [1.86],
            marker="^",
            s=80,
            color=SURFACE if h else c,
            edgecolor=c,
            lw=1.5,
            zorder=3,
            clip_on=False,
        )
    ax.text(0.02, 1.88, "▲ = slope from panel c", fontsize=13, color=INK2, va="bottom")
    ax.set_xlim(-0.03, 1.53)
    ax.set_ylim(1.85, 3.0)
    ax.set_title("d  Noise at each strength")
    ax.set_xlabel("strength α (rr^α); dots = best α")
    ax.set_ylabel("per-gene noise (mean z²)")
    style_axis(ax, grid_axis="y")

    # (e) Best strength vs calibration slope.
    ax = axes[2, 0]
    ax.plot([0, 1.3], [0, 1.3], **THEORY)
    ax.text(
        0.45,
        0.5,
        "best α = slope",
        rotation=45,
        ha="center",
        va="center",
        fontsize=13,
        color=INK2,
    )
    pts = []
    for k in keys[1:]:
        a, n = curve(k)
        pts.append((k, met.loc[k, "slope"], a[int(np.argmin(n))]))
    # Versions with the same best α (the α grid is coarse) are spread vertically.
    for alpha in sorted({p[2] for p in pts}):
        same = sorted([p for p in pts if p[2] == alpha], key=lambda p: p[1])
        for m, (k, x, y) in enumerate(same):
            dy = (m - (len(same) - 1) / 2) * 0.06
            c, h, _ = style[k]
            ax.scatter(
                [x],
                [y + dy],
                s=112,
                marker="D" if k == "ridge_qcexcl" else "o",
                color=SURFACE if h else c,
                edgecolor=c,
                lw=1.8,
                zorder=3,
            )
    ax.set_xlim(0, 1.3)
    ax.set_ylim(0, 1.3)
    ax.set_aspect("equal", adjustable="box")
    ax.text(
        0.03,
        1.25,
        "versions with the same best α\nare spread vertically",
        fontsize=12.5,
        color=INK2,
        va="top",
    )
    ax.set_title("e  Panels c and d agree: best α ≈ slope")
    ax.set_xlabel("calibration slope (panel c)")
    ax.set_ylabel("best strength α (panel d)")
    style_axis(ax, grid_axis="both")

    # (f) chrX/Y, from Ruchit's v2 report: no rr there costs almost nothing.
    ax = axes[2, 1]
    xnames = [
        "No correction",
        "Legacy syn_adj_r\n(no rr on X/Y)",
        "Ruchit 17k, unmasked\n(3 of 13 contexts)",
        "Ruchit 17k, masked",
        "Ridge and the rate\nversions (no rr on X/Y)",
    ]
    xvals = [4.364, 4.364, 4.332, 4.435, met.loc["ridge_qcexcl", "xy_noise"]]
    xcols = [NONE_C, LEGACY_C, RUCHIT_C, RUCHIT_C, REC_C]
    xy = list(range(len(xnames)))[::-1]
    ax.axvline(4.364, color=MUTED, lw=1, ls=(0, (3, 3)))
    full = 4.364 * 4.447 / 5.081
    ax.axvline(full, **THEORY)
    ax.text(
        full + 0.01,
        len(xnames) - 0.7,
        "if rr helped X/Y as much\nas autosomes" " (−12%)",
        fontsize=13,
        color=INK2,
        va="center",
    )
    for y, v, c in zip(xy, xvals, xcols):
        ax.scatter([v], [y], s=112, color=c, edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.text(v + 0.006, y + 0.22, f"{v:.3f}", fontsize=14)
    ax.set_yticks(xy, xnames)
    ax.set_xlim(3.78, 4.48)
    ax.set_ylim(-0.6, len(xnames) - 0.4)
    ax.set_title("f  chrX/Y noise (769 transcripts)")
    ax.set_xlabel("mean z², overall level removed (lower = better)")
    style_axis(ax)

    fig.savefig(
        os.path.join(OUT_DIR, "rr_performance.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_deciles(
    summary="syn_oe_by_decile.full_run.tsv", out="rr_by_decile.png", note=None
):
    """
    Draw syn obs/exp by LOEUF and s_het decile from Ruchit's evaluation export.

    :param summary: Decile summary TSV (decile_explained_removed.py, with the full
        run's and the per-variant rates' per-transcript tables).
    :param out: Output PNG name.
    :param note: Optional line drawn above the panels.
    """
    df = pd.read_csv(os.path.join(OUT_DIR, summary), sep="\t")
    n = int(df[(df.binning == "shet") & (df.correction == "none")].n.sum())

    def oe(binning, correction):
        d = df[(df.binning == binning) & (df.correction == correction)]
        return d.sort_values("decile").oe.tolist()

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    fig.subplots_adjust(
        left=0.08, right=0.97, top=0.94, bottom=0.07, hspace=0.75, wspace=0.28
    )
    dec = list(range(10))
    short = {
        "1e-2_a0.625_full": "Ruchit's latest pick",
        "bonferroni_full": "Ruchit's pipeline",
        "homebuilt_ridge": "Home-built + ridge",
    }
    series = [
        (short.get(k, label.replace("\n", " ")), c, ls, k)
        for k, label, c, _, ls in PERF_VERSIONS
    ]

    for ax, binning, title, xlabel in (
        (
            axes[0, 0],
            "loeuf",
            "g  Synonymous obs/exp by LOEUF decile",
            "LOEUF decile (v4.1.2); 0 = most constrained",
        ),
        (
            axes[0, 1],
            "shet",
            "h  Synonymous obs/exp by s_het decile",
            "gnomAD v4 s_het decile; 0 = most constrained",
        ),
    ):
        ax.axhline(1.0, **THEORY)
        ax.text(9.35, 0.997, "expected", fontsize=13, color=INK2, va="top")
        for label, c, ls, corr in series:
            v = oe(binning, corr)
            ax.plot(
                dec,
                v,
                color=c,
                lw=2.8 if corr == "ridge_qcexcl" else 2,
                ls=ls,
                label=label,
            )
            ax.scatter(
                [9], [v[-1]], s=58, color=c, edgecolor=SURFACE, linewidth=1.5, zorder=3
            )
        ax.set_xticks(dec)
        ax.set_xlim(-0.3, 10.6)
        ax.set_ylim(0.84, 1.14)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        style_axis(ax, grid_axis="y")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.555),
        ncol=4,
        frameon=False,
        fontsize=13,
    )

    for ax, corrs, title in (
        (
            axes[1, 0],
            [("none", NONE_C, "")],
            "i  No correction: LOEUF vs s_het binning",
        ),
        (
            axes[1, 1],
            [
                ("ridge_qcexcl", REC_C, "Ridge"),
                ("roulette_ridge", ROULETTE_C, "Roulette + ridge"),
            ],
            "j  With rr: LOEUF vs s_het binning",
        ),
    ):
        ax.axhline(1.0, **THEORY)
        ax.text(0, 0.996, "expected = 1", fontsize=13, color=INK2, va="top")
        for corr, c, name in corrs:
            lw = 2.8 if corr == "ridge_qcexcl" else 2
            for binning, ls, what in (
                ("loeuf", "-", "by LOEUF"),
                ("shet", (0, (4, 2.5)), "by s_het"),
            ):
                v = oe(binning, corr)
                ax.plot(
                    dec,
                    v,
                    color=c,
                    lw=lw,
                    ls=ls,
                    label=f"{name}, {what}" if name else what,
                )
        if len(corrs) == 1:
            loeuf, shet = oe("loeuf", "none"), oe("shet", "none")
            ax.text(9.2, loeuf[-1], "by LOEUF", va="center", fontsize=14, color=INK2)
            ax.text(
                9.2, shet[-1] - 0.004, "by s_het", va="center", fontsize=14, color=INK2
            )
        else:
            ax.legend(
                loc="upper left",
                frameon=False,
                fontsize=12.5,
                bbox_to_anchor=(0.0, 0.97),
            )
        ax.set_xticks(dec)
        ax.set_xlim(-0.3, 11.4)
        ax.set_ylim(0.97, 1.13)
        ax.set_title(title)
        ax.set_xlabel(f"decile; 0 = most constrained (same {n:,} transcripts)")
        style_axis(ax, grid_axis="y")

    if note:
        fig.text(0.08, 0.985, note, fontsize=15, color=INK2, va="top")
    fig.savefig(
        os.path.join(OUT_DIR, out), dpi=170, bbox_inches="tight", pad_inches=0.3
    )
    plt.close(fig)


FULL_RUN_VERSIONS = [
    ("none", "No correction", NONE_C, None, "-"),
    ("bonferroni_full", "Current method", RUCHIT_C, None, "-"),
    (
        "bonferroni_qcexcl",
        "Current method, QC out",
        RUCHIT_C,
        "///",
        (0, (4, 2.5)),
    ),
    ("ridge_full", "Ridge", REC_C, "///", (0, (4, 2.5))),
    (
        "ridge_qcexcl",
        "Ridge, QC out (recommended)",
        REC_C,
        None,
        "-",
    ),
    ("roulette_ridge", "Roulette + ridge", ROULETTE_C, None, "-"),
    ("homebuilt_ridge", "Home-built + ridge", ROULETTE_C, "///", (0, (4, 2.5))),
]
"""The full pipeline runs: key, label, colour, bar hatch, line style."""

FULL_RUN_LINES = [
    "none",
    "bonferroni_full",
    "ridge_qcexcl",
    "roulette_ridge",
    "homebuilt_ridge",
]
"""The versions drawn as lines (the QC-exclusion twins are within 0.01 of these)."""


def full_run_style(key):
    """Get a full-run version's label, colour, hatch and line style."""
    return {v[0]: v[1:] for v in FULL_RUN_VERSIONS}[key]


def hbar(ax, y, v, left, c, h):
    """Draw one horizontal bar, hatched in white if ``h``."""
    ax.barh(
        y,
        v,
        left=left,
        height=0.6,
        color=c,
        hatch=h,
        hatchcolor="white",
        edgecolor=SURFACE,
        linewidth=1.5,
    )


def draw_new_metrics():
    """Draw the checks added since Ruchit's reports, for every version."""
    met = read_cell_scores("metrics").set_index("version")
    hot = pd.read_csv(os.path.join(OUT_DIR, "downstream.hotspots.tsv"), sep="\t")
    keys = [v[0] for v in PERF_VERSIONS]
    short = {
        "none": "No correction",
        "adj_r": "Gnocchi adj_r",
        "syn_adj_r_norm": "Legacy syn_adj_r",
        "bonferroni_full": "Ruchit's pipeline",
        "1e-2_a0.625_full": "Ruchit's latest pick",
        "ridge_qcexcl": "Ridge (recommended)",
        "roulette_ridge": "Roulette + ridge",
        "homebuilt_ridge": "Home-built + ridge",
    }
    names = [short[v[0]] for v in PERF_VERSIONS]
    style = {v[0]: v[2:] for v in PERF_VERSIONS}
    fig, axes = plt.subplots(2, 2, figsize=(14, 12.5))
    fig.subplots_adjust(
        left=0.18, right=0.97, top=0.89, bottom=0.08, hspace=0.5, wspace=0.5
    )
    fig.suptitle(
        "The checks added since Ruchit's reports",
        x=0.03,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    fig.text(
        0.03,
        0.94,
        "a–c: cell-level scoring of every version, as in the figure above."
        "  d: the constraint pipeline's own counts.",
        fontsize=13.5,
        color=INK2,
    )

    # (a) Noise on held-out chromosomes, and (b) on genes without a QC flag.
    for ax, col, title, xlim in (
        (axes[0, 0], "oof_noise", "a  Noise on held-out chromosomes", (1.8, 3.0)),
        (
            axes[0, 1],
            "outcome_free_noise",
            "b  Noise on genes without a QC flag",
            (1.8, 6.0),
        ),
    ):
        ks = [
            k
            for k in keys
            if col != "oof_noise" or pd.notna(met.loc[k, col]) or k == "none"
        ]
        yy = list(range(len(ks)))[::-1]
        for y, k in zip(yy, ks):
            c, h, _ = style[k]
            v = (
                met.loc[k, col]
                if pd.notna(met.loc[k, col])
                else met.loc[k, "core_noise"]
            )
            allc = met.loc[k, "core_noise"] if col == "oof_noise" else None
            if allc is not None and k != "none":
                ax.plot([allc, v], [y, y], color=GRID, lw=2, zorder=1)
                ax.scatter([allc], [y], s=150, marker="|", color=INK2, lw=2, zorder=3)
            ax.scatter(
                [v],
                [y],
                s=112,
                color=SURFACE if h else c,
                edgecolor=c,
                lw=1.8,
                zorder=4,
            )
            ax.text(
                v + 0.03 * (xlim[1] - xlim[0]), y, f"{v:.2f}", va="center", fontsize=13
            )
        ax.axvline(
            met.loc["none", col if col != "oof_noise" else "core_noise"],
            color=MUTED,
            lw=1,
            ls=(0, (3, 3)),
            zorder=0,
        )
        ax.set_yticks(yy, [names[keys.index(k)] for k in ks])
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.6, len(ks) - 0.4)
        ax.set_title(title)
        ax.set_xlabel("mean z², level removed (lower = better)")
        style_axis(ax)
    axes[0, 0].set_xlabel("mean z², level removed; ● held out, | all chromosomes")

    # (c) Calibration slope by context class.
    ax = axes[1, 0]
    ks = keys[1:]
    yy = list(range(len(ks)))[::-1]
    ax.axvline(1.0, **THEORY)
    for j, (cls, c) in enumerate(CLASS_C.items()):
        for y, k in zip(yy, ks):
            v = met.loc[k, f"slope_{cls}"]
            if pd.notna(v):
                ax.scatter(
                    [v],
                    [y + (1 - j) * 0.26],
                    s=95,
                    color=c,
                    edgecolor=SURFACE,
                    lw=1.1,
                    zorder=3,
                    label=cls if y == yy[0] else None,
                )
    ax.set_xscale("log")
    ticks = [0.1, 0.2, 0.5, 1, 2, 4]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0.08, 4.5)
    ax.set_yticks(yy, names[1:])
    ax.set_ylim(-0.6, len(ks) - 0.2)
    ax.set_title("c  Slope by context class")
    ax.set_xlabel("← too strong     slope (log)     too weak →")
    ax.legend(loc="lower right", frameon=False, fontsize=14, handletextpad=0.2)
    style_axis(ax)

    # (d) Maternal hotspot regions.
    ax = axes[1, 1]
    regions = list(dict.fromkeys(hot.region))
    show = {
        "none": ("No correction", NONE_C, "o", False),
        "ridge_qcexcl": ("Ridge (recommended)", REC_C, "D", False),
        "roulette_ridge": ("Roulette + ridge", ROULETTE_C, "s", False),
        "homebuilt_ridge": ("Home-built + ridge", ROULETTE_C, "^", True),
    }
    for i, region in enumerate(regions):
        y = len(regions) - 1 - i
        for j, (k, (label, c, m, hollow)) in enumerate(show.items()):
            r = hot[
                (hot.region == region)
                & (hot.set == "control genes")
                & (hot.version == k)
            ]
            ax.scatter(
                r.oe,
                [y + (1.5 - j) * 0.2] * len(r),
                s=95,
                marker=m,
                color=SURFACE if hollow else c,
                edgecolor=c,
                lw=1.6,
                zorder=3,
                label=label if i == 0 else None,
            )
    ax.axvline(1.0, **THEORY)
    ax.set_yticks(range(len(regions))[::-1], regions)
    ax.set_xlim(0.9, 1.36)
    ax.set_title("d  Maternal hotspot regions")
    ax.set_xlabel(
        "synonymous obs/exp of the region's genes\n(no QC flag, not"
        " outliers), relative to all genes"
    )
    ax.legend(loc="lower right", frameon=False, fontsize=13, handletextpad=0.2)
    style_axis(ax)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_new_metrics.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_full_run_release():
    """Draw what the full pipeline runs do to LoF and missense metrics."""
    sc = pd.read_csv(os.path.join(OUT_DIR, "release.scores.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(OUT_DIR, "release.deciles.tsv"), sep="\t")
    sc = sc.set_index("version")
    fig, axes = plt.subplots(2, 2, figsize=(15, 11.5))
    fig.subplots_adjust(
        left=0.2, right=0.97, top=0.86, bottom=0.08, hspace=0.45, wspace=0.42
    )
    fig.suptitle(
        "What each version does to LoF and missense metrics",
        x=0.03,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    fig.text(
        0.03,
        0.925,
        "17,641 autosomal MANE Select transcripts. LOEUF still separates the 314 ClinGen"
        " haploinsufficient genes\nequally well in every version (AUC 0.867–0.870).",
        fontsize=13.5,
        color=INK2,
    )
    base = dec[dec.version == "none"].set_index(["csq", "decile_none"]).oe
    for ax, csq, title in (
        (axes[0, 0], "lof", "a  LoF obs/exp, relative to no rr"),
        (axes[0, 1], "mis", "b  Missense obs/exp, relative to no rr"),
    ):
        ax.axhline(1.0, **THEORY)
        for key in FULL_RUN_LINES[1:]:
            label, c, _, ls = full_run_style(key)
            v = dec[(dec.version == key) & (dec.csq == csq)].sort_values("decile_none")
            ax.plot(
                v.decile_none,
                v.oe.to_numpy() / base.loc[csq].reindex(v.decile_none).to_numpy(),
                color=c,
                ls=ls,
                lw=2.8 if key == "ridge_qcexcl" else 2,
                label=label.replace("\n", " "),
            )
        ax.set_xticks(range(10))
        ax.set_xlim(-0.3, 9.3)
        ax.set_ylim(0.9, 1.1)
        ax.set_title(title)
        ax.set_xlabel("LOEUF decile without rr; 0 = most constrained")
        style_axis(ax, grid_axis="y")
    axes[0, 0].set_ylabel("obs/exp ÷ obs/exp without rr\n(above 1 = expected fell)")
    axes[0, 0].legend(loc="lower left", frameon=False, fontsize=12.5)

    keys = [v[0] for v in FULL_RUN_VERSIONS]
    names = [full_run_style(k)[0] for k in keys]
    ypos = list(range(len(keys)))[::-1]

    # (c) Genes crossing LOEUF 0.6.
    ax = axes[1, 0]
    for y, k in zip(ypos, keys):
        if k == "none":
            continue
        _, c, h, _ = full_run_style(k)
        out, into = int(sc.loc[k, "out_of_loeuf_lt_0.6"]), int(
            sc.loc[k, "into_loeuf_lt_0.6"]
        )
        hbar(ax, y, -out, 0, c, h)
        hbar(ax, y, into, 0, c, h)
        ax.text(-out - 12, y, f"−{out}", ha="right", va="center", fontsize=13)
        ax.text(into + 12, y, f"+{into}", va="center", fontsize=13)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(ypos, names)
    ax.set_xlim(-720, 300)
    ax.set_ylim(-0.6, len(keys) - 0.1)
    ax.set_title("c  Genes crossing LOEUF 0.6")
    ax.set_xlabel("← no longer below 0.6      newly below 0.6 →")
    style_axis(ax)

    # (d) Missense-constrained genes and the z-score scale.
    ax = axes[1, 1]
    for y, k in zip(ypos, keys):
        _, c, h, _ = full_run_style(k)
        n = int(sc.loc[k, "mis_z_ge_3.09"])
        hbar(ax, y, n, 0, c, h)
        ax.text(
            n + 60,
            y,
            f"{n:,}  (÷ {sc.loc[k, 'mis_z_scale']:.2f})",
            va="center",
            fontsize=13,
        )
    ax.axvline(sc.loc["none", "mis_z_ge_3.09"], color=MUTED, lw=1, ls=(0, (3, 3)))
    ax.set_yticks(ypos, [""] * len(keys))
    ax.set_xlim(0, 6200)
    ax.set_ylim(-0.6, len(keys) - 0.1)
    ax.set_title("d  Missense z ≥ 3.09")
    ax.set_xlabel("genes (÷ the spread raw z is divided by)")
    style_axis(ax)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_full_run_release.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


BG_PIPELINE_RUNS = [
    ("absent_locus", "Absent from gnomAD,\none row per locus (current)"),
    ("absent_variant", "Absent from gnomAD,\none row per variant"),
    ("all_locus", "All sites,\none row per locus"),
    ("all_variant", "All sites,\none row per variant"),
    ("ctx_proportional", "Absent, sampled by context\nin proportion to all sites"),
    ("ctx_even", "Absent, sampled\nevenly by context"),
]
"""The pipeline background runs (run_backgrounds.sh), in figure order."""


def draw_background_pipeline():
    """Draw the pipeline background runs' noise and class slopes."""
    from matplotlib.lines import Line2D

    m = pd.read_csv(os.path.join(OUT_DIR, "rr_bg_pipeline.metrics.tsv"), sep="\t")
    m = m.set_index("version")
    fig, axes = plt.subplots(1, 3, figsize=(18, 8), sharey=True)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.76, bottom=0.12, wspace=0.14)
    fig.suptitle(
        "The training backgrounds, rerun through the pipeline",
        x=0.02,
        ha="left",
        fontsize=19,
        weight="bold",
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                ls="",
                marker="o",
                ms=12,
                color=MUTED,
                mec=SURFACE,
                label="default penalty (0.3 per site)",
            ),
            Line2D(
                [],
                [],
                ls="",
                marker="o",
                ms=12,
                mfc=SURFACE,
                mec=MUTED,
                mew=1.8,
                label="penalty per DNM (10; same shrinkage for any background size)",
            ),
            Line2D(
                [],
                [],
                ls=(0, (3, 3)),
                color=REC_C,
                label="current background, default penalty",
            ),
        ],
        loc="upper left",
        bbox_to_anchor=(0.015, 0.935),
        ncol=3,
        frameon=False,
        fontsize=14,
    )
    ypos = list(range(len(BG_PIPELINE_RUNS)))[::-1]
    names = [n for _, n in BG_PIPELINE_RUNS]
    for ax, col, title, xlim in (
        (axes[0], "oof_noise", "a  Noise, held-out chromosomes", (2.42, 2.66)),
        (
            axes[1],
            "outcome_free_noise",
            "b  Noise, genes without a QC flag",
            (2.74, 3.1),
        ),
    ):
        ref = m.loc["absent_locus_ridge", col]
        ax.axvline(ref, color=REC_C, lw=1.2, ls=(0, (3, 3)), zorder=0)
        for y, (run, _) in zip(ypos, BG_PIPELINE_RUNS):
            a, b = m.loc[f"{run}_ridge", col], m.loc[f"{run}_ridge_dnm_10", col]
            c = REC_C if run.startswith("absent_locus") else MUTED
            ax.plot([a, b], [y, y], color=GRID, lw=2.5, zorder=1)
            ax.scatter([a], [y], s=150, color=c, edgecolor=SURFACE, lw=1.5, zorder=3)
            ax.scatter([b], [y], s=150, color=SURFACE, edgecolor=c, lw=2, zorder=3)
        ax.set_xlim(*xlim)
        ax.set_title(title)
        ax.set_xlabel("mean z², level removed (lower = better)")
        style_axis(ax)
    axes[0].set_yticks(ypos, names)
    axes[0].set_ylim(-0.6, len(ypos) - 0.4)

    ax = axes[2]
    ax.axvline(1.0, **THEORY)
    for j, (cls, c) in enumerate(CLASS_C.items()):
        for y, (run, _) in zip(ypos, BG_PIPELINE_RUNS):
            ax.scatter(
                [m.loc[f"{run}_ridge", f"slope_{cls}"]],
                [y + (1 - j) * 0.26],
                s=120,
                color=c,
                edgecolor=SURFACE,
                lw=1.2,
                zorder=3,
                label=cls if y == ypos[0] else None,
            )
    ax.set_xscale("log")
    ticks = [0.1, 0.2, 0.5, 1, 2, 4]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0.09, 4.5)
    ax.set_title("c  Slope by context class (default penalty)")
    ax.set_xlabel("← too strong     slope (log)     too weak →")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, 1.06),
        ncol=3,
        frameon=False,
        fontsize=14,
        handletextpad=0.2,
        columnspacing=1.0,
    )
    style_axis(ax)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_background_pipeline.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_background_bias(prefix="rr_logistic_fits.pipeline_bg6_dnm10"):
    """Draw what the current background drops, and how each background changes rr."""
    bins = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.bins.tsv"), sep="\t")
    cls = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.contexts.tsv"), sep="\t")
    bins["cls"] = bins.context.map(cls.set_index("context")["class"])

    def by_bin(g, bg):
        d = g[g.background == bg].assign(w=lambda x: x.possible * x.rr)
        return d.groupby("bin")[["w", "possible", "background_weight"]].sum()

    fig, axes = plt.subplots(1, 2, figsize=(17, 7))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.8, bottom=0.14, wspace=0.42)
    fig.suptitle(
        "What the current background leaves out, and what that does to rr",
        x=0.07,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    fig.text(
        0.07,
        0.875,
        "Pipeline background runs, ridge with a penalty of 10 per DNM; deciles of"
        " possible sites by their rr under the current background.",
        fontsize=13.5,
        color=INK2,
    )
    ax = axes[0]
    for c, col in CLASS_C.items():
        g = bins[bins.cls == c]
        cur, al = by_bin(g, "absent"), by_bin(g, "all_sites")
        x = cur.w / cur.possible
        drop = 100 * (1 - cur.background_weight / al.background_weight)
        ax.plot(x, drop, color=col, lw=2.4, marker="o", ms=7)
        ax.text(
            x.iloc[-1] + 0.01, drop.iloc[-1], c, va="center", fontsize=14, color=INK
        )
    ax.set_xlabel("rr under the current background (decile mean)")
    ax.set_ylabel("% of possible sites left out\n(seen in gnomAD exomes)")
    ax.set_title("a  Share left out, by rr")
    ax.set_ylim(0, 50)
    style_axis(ax, grid_axis="y")

    ax = axes[1]
    others = [b for b in LOGIT_BACKGROUNDS if b[0] != "absent"]
    for i, (c, col) in enumerate(CLASS_C.items()):
        g = bins[bins.cls == c]
        cur = by_bin(g, "absent")
        x = np.log(cur.w / cur.possible)
        for j, (bg, label, bc, ls) in enumerate(others):
            o = by_bin(g, bg)
            slope = np.polyfit(x, np.log(o.w / o.possible), 1)[0]
            ax.scatter(
                [slope],
                [i + (j - 2) * 0.13],
                s=90,
                color=bc,
                marker="o" if ls == "-" else "s" if ls == (0, (4, 2.5)) else "D",
                edgecolor=SURFACE,
                lw=1.2,
                zorder=3,
                label=label if i == 0 else None,
            )
    ax.axvline(1.0, **THEORY)
    ax.set_yticks(range(len(CLASS_C)), list(CLASS_C))
    ax.invert_yaxis()
    ax.set_xlim(0.6, 1.4)
    ax.set_xlabel("steepness relative to the current background")
    ax.set_title("b  How steep each background makes rr")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=14)
    style_axis(ax)
    fig.subplots_adjust(right=0.7)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_background.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_feature_use():
    """Draw, in the abstract, how each version turns the 52 feature columns into rr."""
    W, H = 1500, 690
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    rng = np.random.default_rng(3)
    n_feat, n_scale, cell = 13, 4, 15
    gx = 300
    # Illustrative only: which cells each p-value screen keeps, and ridge weights.
    screen_keep = {(1, 0), (5, 2), (9, 3)}
    loose_keep = screen_keep | {(3, 1), (7, 0), (11, 2)}
    ridge_w = rng.uniform(0.15, 0.6, (n_scale, n_feat))
    rows = [
        (
            "Ruchit's pipeline\n(current method)",
            RUCHIT_C,
            lambda f, k: 1.0 if (f, k) in screen_keep else 0.0,
            [
                "screen each column on its own;",
                "keep it if p < 0.05 / 52",
                "(one scale per feature)",
            ],
            ["fit the kept columns", "with no penalty"],
            ["rr as fitted"],
            False,
        ),
        (
            "Ruchit's latest pick",
            RUCHIT_C,
            lambda f, k: 1.0 if (f, k) in loose_keep else 0.0,
            [
                "the same screen with a",
                "looser cutoff (p < 0.01):",
                "more columns kept",
            ],
            ["fit the kept columns", "with no penalty"],
            ["rr^α, α = 0.625", "(weaker)"],
            False,
        ),
        (
            "Ridge (recommended)",
            REC_C,
            lambda f, k: ridge_w[k, f],
            ["no screen:", "every column goes in"],
            [
                "fit all columns together;",
                "one penalty shrinks every",
                "weight toward 0",
            ],
            ["rr as fitted"],
            False,
        ),
        (
            "Per-variant rates\n+ ridge",
            ROULETTE_C,
            lambda f, k: ridge_w[k, f],
            ["no screen:", "every column goes in"],
            [
                "the same ridge fit, with",
                "log(rate) as a fixed offset:",
                "it learns what rates miss",
            ],
            ["expected × rate × rr"],
            True,
        ),
    ]
    ax.text(
        30,
        34,
        "How each version turns the 52 feature columns into rr",
        fontsize=15,
        weight="bold",
        va="center",
    )
    ax.text(
        30,
        62,
        "Each grid is one context's 13 features (across) × 4 window sizes (down);"
        " shading = how much a column counts (illustrative). In CpG contexts every version"
        " also leaves out the 5 methylation-linked features.",
        fontsize=11,
        color=INK2,
        va="center",
    )
    for r, (name, col, weight, step1, step2, out, offset) in enumerate(rows):
        y0 = 110 + r * 128
        ax.add_patch(
            FancyBboxPatch(
                (30, y0 + 8),
                220,
                70,
                boxstyle="round,pad=0,rounding_size=6",
                facecolor="white",
                edgecolor=col,
                lw=2 if col == REC_C else 1.5,
            )
        )
        ax.text(
            140,
            y0 + 43,
            name,
            ha="center",
            va="center",
            fontsize=12.5,
            weight="bold" if col == REC_C else "normal",
        )
        for f in range(n_feat):
            for k in range(n_scale):
                w = weight(f, k)
                ax.add_patch(
                    Rectangle(
                        (gx + f * (cell + 2), y0 + 12 + k * (cell + 2)),
                        cell,
                        cell,
                        facecolor=col if w > 0 else SURFACE,
                        alpha=max(w, 0.0) if w > 0 else 1.0,
                        edgecolor=AXIS if w == 0 else SURFACE,
                        lw=0.8,
                    )
                )
        gw = n_feat * (cell + 2)
        arrow(ax, gx + gw + 8, y0 + 43, gx + gw + 40, y0 + 43)
        box(ax, gx + gw + 46, y0 + 6, 285, 74, SAME_FILL, SAME_EDGE, step1, size=11.5)
        arrow(ax, gx + gw + 335, y0 + 43, gx + gw + 367, y0 + 43)
        box(ax, gx + gw + 373, y0 + 6, 290, 74, RR_FILL, RR_EDGE, step2, size=11.5)
        if offset:
            ax.add_patch(
                FancyBboxPatch(
                    (gx + gw + 373, y0 + 96),
                    290,
                    28,
                    boxstyle="round,pad=0,rounding_size=5",
                    facecolor=RATE_FILL,
                    edgecolor=RATE_EDGE,
                    lw=1.3,
                )
            )
            ax.text(
                gx + gw + 518,
                y0 + 110,
                "per-variant rate (Roulette or home-built)",
                ha="center",
                va="center",
                fontsize=10.5,
                color=INK,
            )
            arrow(ax, gx + gw + 518, y0 + 95, gx + gw + 518, y0 + 82)
        arrow(ax, gx + gw + 667, y0 + 43, gx + gw + 699, y0 + 43)
        box(
            ax,
            gx + gw + 705,
            y0 + 18,
            190,
            50,
            "white",
            col,
            out,
            size=11.5,
            center=True,
        )
    # Legend.
    ly = H - 30
    items = [
        (INK2, 1.0, "column used at full weight"),
        (INK2, 0.35, "used, shrunk by the penalty"),
        (SURFACE, 1.0, "not used"),
    ]
    lx = gx
    for c, a, label in items:
        ax.add_patch(
            Rectangle(
                (lx, ly - 8),
                cell,
                cell,
                facecolor=c,
                alpha=a,
                edgecolor=AXIS if c == SURFACE else SURFACE,
                lw=0.8,
            )
        )
        ax.text(lx + cell + 8, ly, label, va="center", fontsize=11, color=INK2)
        lx += cell + 8 + len(label) * 6.8 + 30
    fig.savefig(
        os.path.join(OUT_DIR, "rr_feature_use.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_rate_parts():
    """Draw why the home-built rates fall short of Roulette's, part by part."""
    rates = pd.read_csv(os.path.join(OUT_DIR, "rr_rates.metrics.tsv"), sep="\t")
    rates = rates.set_index("version")
    strength = pd.read_csv(os.path.join(OUT_DIR, "rr_rates.strength.tsv"), sep="\t")
    local = pd.read_csv(
        os.path.join(OUT_DIR, "roulette_combo.v2.metrics.tsv"), sep="\t"
    ).set_index("version")
    agree = pd.read_csv(
        os.path.join(OUT_DIR, "compare_rate_parts.agreement.tsv"), sep="\t"
    )
    genes = pd.read_csv(
        os.path.join(OUT_DIR, "compare_rate_parts.genes.tsv.gz"), sep="\t"
    )

    fig, axes = plt.subplots(
        1, 3, figsize=(19, 6.4), gridspec_kw={"width_ratios": [1.25, 1, 0.8]}
    )
    fig.subplots_adjust(left=0.15, right=0.98, top=0.8, bottom=0.15, wspace=0.42)
    fig.suptitle(
        "Why the home-built rates fall short of Roulette's: the regional part",
        x=0.03,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    fig.text(
        0.03,
        0.875,
        "Each part alone, applied to Ruchit's exported counts (no rr; CpG left at 1)."
        " Roulette's parts are split locally from its published rates.",
        fontsize=13.5,
        color=INK2,
    )

    # (a) Noise of each part alone.
    ax = axes[0]
    best = strength[strength.version == "homebuilt_only"].core_noise.min()
    rows = [
        (
            "Roulette, both parts",
            rates.loc["roulette_only", "core_noise"],
            ROULETTE_C,
            None,
        ),
        (
            "Roulette, regional part",
            local.loc["roulette_regional_non_cpg", "core_noise"],
            ROULETTE_C,
            None,
        ),
        (
            "Roulette, sequence part",
            local.loc["roulette_sequence_non_cpg", "core_noise"],
            ROULETTE_C,
            None,
        ),
        (
            "Home-built, both parts",
            rates.loc["homebuilt_only", "core_noise"],
            ROULETTE_C,
            "///",
        ),
        (
            "Home-built, regional part",
            rates.loc["homebuilt_regional", "core_noise"],
            ROULETTE_C,
            "///",
        ),
        (
            "Home-built, sequence part",
            rates.loc["homebuilt_sequence", "core_noise"],
            ROULETTE_C,
            "///",
        ),
    ]
    yy = list(range(len(rows)))[::-1]
    for y, (label, v, c, h) in zip(yy, rows):
        ax.barh(y, v, height=0.6, color=c, hatch=h, hatchcolor="white", linewidth=0)
        note = (
            f"  (| {best:.2f} at α 0.75)" if label == "Home-built, both parts" else ""
        )
        ax.text(v + 0.03, y, f"{v:.2f}{note}", va="center", fontsize=14)
    ax.scatter([best], [yy[3]], marker="|", s=300, color=INK, lw=2.5, zorder=4)
    ax.axvline(rates.loc["none", "core_noise"], color=MUTED, lw=1.2, ls=(0, (3, 3)))
    ax.text(
        rates.loc["none", "core_noise"] - 0.03,
        len(rows) - 0.45,
        "no correction",
        ha="right",
        fontsize=13,
        color=INK2,
    )
    ax.set_yticks(yy, [r[0] for r in rows])
    ax.set_xlim(1.5, 3.1)
    ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.set_title("a  Per-gene noise of each part alone")
    ax.set_xlabel("mean z², level removed (lower = better)")
    style_axis(ax)

    # (b) Gene-level regional factor, home-built against Roulette.
    ax = axes[1]
    g = genes[genes.exp >= 20]
    hot = g.chromosome.eq("chr8") & g.mid.lt(20_000_000)
    ax.scatter(
        g.roulette_regional[~hot],
        g.rr_homebuilt_regional[~hot],
        s=5,
        color=MUTED,
        alpha=0.25,
        lw=0,
        rasterized=True,
    )
    ax.scatter(
        g.roulette_regional[hot],
        g.rr_homebuilt_regional[hot],
        s=36,
        color=RUCHIT_C,
        edgecolor=SURFACE,
        lw=0.6,
        zorder=3,
        label="genes in the first 20 Mb of 8p",
    )
    lim = (0.7, 1.7)
    ax.plot(lim, lim, **THEORY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [0.8, 1.0, 1.25, 1.6]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    a = agree[
        (agree.level == "gene") & (agree["class"] == "all") & (agree.part == "regional")
    ].iloc[0]
    ax.text(
        0.04,
        0.96,
        f"correlation {a['corr']:.2f}\nhome-built spread\n{a.sd_homebuilt / a.sd_roulette:.2f}× Roulette's",
        transform=ax.transAxes,
        va="top",
        fontsize=13.5,
        color=INK2,
    )
    ax.set_title("b  Each gene's regional factor")
    ax.set_xlabel("Roulette (gene mean, non-CpG sites)")
    ax.set_ylabel("home-built")
    ax.legend(loc="lower right", frameon=False, fontsize=13)
    style_axis(ax, grid_axis="both")

    # (c) Spread of the gene-level regional factor by context class.
    ax = axes[2]
    cl = ["A/T", "C non-CpG"]
    for i, c in enumerate(cl):
        r = agree[
            (agree.level == "gene") & (agree["class"] == c) & (agree.part == "regional")
        ].iloc[0]
        for j, (who, v, h) in enumerate(
            (("Roulette", r.sd_roulette, None), ("home-built", r.sd_homebuilt, "///"))
        ):
            ax.bar(
                i + (j - 0.5) * 0.36,
                100 * v,
                width=0.34,
                color=ROULETTE_C,
                hatch=h,
                hatchcolor="white",
                linewidth=0,
                label=who if i == 0 else None,
            )
            ax.text(
                i + (j - 0.5) * 0.36,
                100 * v + 0.2,
                f"{100 * v:.1f}%",
                ha="center",
                fontsize=13,
            )
    ax.set_xticks(range(len(cl)), cl)
    ax.set_ylabel("spread across genes (SD of log factor, %)")
    ax.set_title("c  Spread by context class")
    ax.legend(loc="upper left", frameon=False, fontsize=13)
    ax.set_ylim(0, 12.5)
    style_axis(ax, grid_axis="y")
    fig.savefig(
        os.path.join(OUT_DIR, "rr_rate_parts.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_scorecard():
    """Draw the "what a good outcome looks like" scorecard for every candidate."""
    met = read_cell_scores("metrics").set_index("version")
    rel = pd.read_csv(os.path.join(OUT_DIR, "release.scores.tsv"), sep="\t")
    rel = rel.set_index("version")
    hot = pd.read_csv(os.path.join(OUT_DIR, "downstream.hotspots.tsv"), sep="\t")
    hot = hot[(hot.region == "8p, first 20 Mb") & (hot.set == "control genes")]
    hot = hot.set_index("version").oe
    cols = [
        ("bonferroni_full", "Ruchit's pipeline\n(current method)", RUCHIT_C),
        ("1e-2_a0.625_full", "Ruchit's latest pick\n(1e-2, α 0.625)", RUCHIT_C),
        ("ridge_qcexcl", "Ridge\n(recommended)", REC_C),
        ("roulette_ridge", "Roulette + ridge\n(quick first try)", ROULETTE_C),
        ("homebuilt_ridge", "Home-built + ridge\n(quick first try)", ROULETTE_C),
    ]
    none = met.loc["none"]

    def gain(k, col, best):
        """Share of the no-correction gap closed, relative to the best version."""
        return (none[col] - met.loc[k, col]) / (none[col] - best)

    best_noise = min(met.loc[k, "core_noise"] for k, *_ in cols)
    best_oof = min(met.loc[k, "oof_noise"] for k, *_ in cols)
    best_free = min(met.loc[k, "outcome_free_noise"] for k, *_ in cols)
    best_out = min(met.loc[k, "outliers_no_qc_flag"] for k, *_ in cols)
    rows = []
    rows.append(
        (
            "Keep the level the calibration sets",
            "mean rr over exome synonymous sites = 1",
            [("normalized", 1.0)] * 5,
        )
    )
    rows.append(
        (
            "Reduce per-gene noise",
            f"mean z²; no correction {none.core_noise:.2f}",
            [
                (f"{met.loc[k, 'core_noise']:.2f}", gain(k, "core_noise", best_noise))
                for k, *_ in cols
            ],
        )
    )
    rows.append(
        (
            "Hold up on held-out chromosomes",
            "noise, models fit on the other half",
            [
                (
                    f"{met.loc[k, 'oof_noise']:.2f}",
                    (none.core_noise - met.loc[k, "oof_noise"])
                    / (none.core_noise - best_oof),
                )
                for k, *_ in cols
            ],
        )
    )
    rows.append(
        (
            "Hold up on genes it doesn't select on",
            f"noise, genes without a QC flag; none {none.outcome_free_noise:.2f}",
            [
                (
                    f"{met.loc[k, 'outcome_free_noise']:.2f}",
                    gain(k, "outcome_free_noise", best_free),
                )
                for k, *_ in cols
            ],
        )
    )
    rows.append(
        (
            "Remove synonymous outliers",
            f"outliers, genes without a QC flag; none {int(none.outliers_no_qc_flag)}",
            [
                (
                    f"{int(met.loc[k, 'outliers_no_qc_flag'])}",
                    gain(k, "outliers_no_qc_flag", best_out),
                )
                for k, *_ in cols
            ],
        )
    )
    rows.append(
        (
            "Not overcorrect",
            "calibration slope; 1 = right strength",
            [
                (
                    f"{met.loc[k, 'slope']:.2f}"
                    + (", too strong" if met.loc[k, "slope"] < 0.7 else ""),
                    1 - min(1, abs(1 - met.loc[k, "slope"]) / 0.45),
                )
                for k, *_ in cols
            ],
        )
    )
    rows.append(
        (
            "Correct the known regional excess",
            f"8p genes, synonymous obs/exp; none {hot['none']:.2f}",
            [
                (
                    (f"{hot[k]:.2f}", 1 - abs(hot[k] - 1) / (hot["none"] - 1))
                    if k in hot
                    else ("not run", None)
                )
                for k, *_ in cols
            ],
        )
    )
    lof = {k: 1 - rel.loc[k, "lof_oe"] / rel.loc[k, "lof_oe"] for k in rel.index}

    def lofmis(k):
        if k not in rel.index:
            return ("not run in\nthe pipeline", None)
        d_lof = rel.loc["none", "lof_oe"] / rel.loc[k, "lof_oe"] - 1
        d_mis = rel.loc["none", "mis_oe"] / rel.loc[k, "mis_oe"] - 1
        txt = f"LoF exp {100 * d_lof:+.1f}%,\nmissense {100 * d_mis:+.1f}%"
        return (txt, 1 - min(1, max(abs(d_lof), abs(d_mis)) / 0.07))

    rows.append(
        (
            "Leave LoF and missense levels alone",
            "expected counts against the run without rr\n(LOEUF separates ClinGen HI genes in all)",
            [lofmis(k) for k, *_ in cols],
        )
    )
    rows.append(
        (
            "Be simple to explain and maintain",
            "settings to choose and defend",
            [
                ("one cutoff, but too\nstrong without α", 0.45),
                ("a cutoff and α", 0.55),
                ("none: penalty set\non held-out DNMs", 1.0),
                ("external rates; its\nstrength still to set", 0.45),
                ("our own build; several\nchoices to tune", 0.25),
            ],
        )
    )

    nr, nc = len(rows), len(cols)
    fig = plt.figure(figsize=(17, 1.3 + 0.82 * (nr + 1)))
    ax = fig.add_axes([0, 0, 1, 1])
    W, H = 17, 1.3 + 0.82 * (nr + 1)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    x0, cw, ch, gap = 5.6, 2.2, 0.68, 0.08
    good = np.array([0.105, 0.62, 0.43])
    cmap = LinearSegmentedColormap.from_list(
        "score", ["#fbe1d6", "#fdf0cc", "#d5efe2", "#9fd9bd"]
    )
    for j, (k, name, c) in enumerate(cols):
        x = x0 + j * (cw + gap)
        ax.add_patch(
            FancyBboxPatch(
                (x, 0.25),
                cw,
                0.75,
                boxstyle="round,pad=0,rounding_size=0.08",
                facecolor="white",
                edgecolor=c,
                lw=2 if k == "ridge_qcexcl" else 1.4,
            )
        )
        ax.text(
            x + cw / 2,
            0.625,
            name,
            ha="center",
            va="center",
            fontsize=12.5,
            weight="bold" if k == "ridge_qcexcl" else "normal",
            color=INK,
        )
    for i, (crit, how, cells) in enumerate(rows):
        y = 1.2 + i * 0.82
        ax.text(
            0.25,
            y + ch / 2 - 0.1,
            crit,
            fontsize=13,
            weight="bold",
            va="center",
            color=INK,
        )
        ax.text(
            0.25,
            y + ch / 2 + 0.17,
            how,
            fontsize=10.5,
            va="center",
            color=INK2,
            linespacing=1.1,
        )
        for j, (txt, score) in enumerate(cells):
            x = x0 + j * (cw + gap)
            fill = AXIS if score is None else cmap(max(0.0, min(1.0, score)))
            ax.add_patch(
                FancyBboxPatch(
                    (x, y),
                    cw,
                    ch,
                    boxstyle="round,pad=0,rounding_size=0.06",
                    facecolor=fill if score is not None else SURFACE,
                    edgecolor=SURFACE if score is not None else AXIS,
                    lw=1.2,
                )
            )
            ax.text(
                x + cw / 2,
                y + ch / 2,
                txt,
                ha="center",
                va="center",
                fontsize=12,
                color=INK2 if score is None else INK,
                linespacing=1.15,
            )
    ax.text(
        0.25,
        H - 0.35,
        "Cell colour: green = closest to the goal among these versions, amber ="
        " partway, red = furthest. Noise, outliers and slope: cell-level scoring of the\nfull"
        " pipeline runs on Ruchit's exported counts. 8p and LoF/missense: the constraint"
        " pipeline's own counts and releases.",
        fontsize=10.5,
        color=INK2,
        va="center",
        linespacing=1.4,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "rr_scorecard.png"),
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


FEATURE_NAMES = {
    "GC_content": "GC content",
    "CpG_island": "CpG islands",
    "met_sperm": "sperm methylation",
    "Nucleosome": "nucleosome occupancy",
    "SINE": "SINEs",
    "LINE": "LINEs",
    "LCR": "low-complexity regions",
    "recomb_male": "male recombination",
    "recomb_female": "female recombination",
    "dist2telo": "distance to telomere",
    "dist2cent": "distance to centromere",
    "cDNM_maternal_05M": "maternal clustered DNMs",
    "cDNM_paternal_05M": "paternal clustered DNMs",
}


def draw_penalty(prefix="penalty_path", context="CCT"):
    """Draw what the ridge penalty is, what it does, and how it is chosen."""
    loss = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.loss.tsv"), sep="\t")
    coefs = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.coefs.tsv"), sep="\t")
    rng_ = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.rr_range.tsv"), sep="\t")
    best = float(loss.loc[loss.heldout_log_loss.idxmin(), "penalty"])

    fig = plt.figure(figsize=(16, 8.4))
    top = fig.add_axes([0, 0.6, 1, 0.4])
    top.set_xlim(0, 1600)
    top.set_ylim(350, 0)
    top.axis("off")
    top.text(
        40,
        40,
        "What the ridge penalty is, and how it is chosen",
        fontsize=19,
        weight="bold",
        va="center",
    )
    top.text(
        40,
        80,
        "Each context's model is fit by minimizing two terms; the penalty λ sets"
        " how much the second one counts.",
        fontsize=15,
        color=INK2,
        va="center",
    )
    box(
        top,
        120,
        120,
        470,
        110,
        SAME_FILL,
        SAME_EDGE,
        ["how badly the model predicts", "which training sites are DNMs", "(log loss)"],
        size=14,
        center=True,
    )
    top.text(640, 175, "+", fontsize=40, ha="center", va="center", color=INK)
    top.text(
        720,
        175,
        "λ  ×",
        fontsize=28.5,
        ha="center",
        va="center",
        color=INK,
        weight="bold",
    )
    box(
        top,
        790,
        120,
        470,
        110,
        RR_FILL,
        RR_EDGE,
        [
            "how big the weights are",
            "(sum of squared weights;",
            "the intercept is free)",
        ],
        size=14,
        center=True,
    )
    top.text(
        355,
        262,
        "small when the weights fit the DNMs,\nincluding their noise",
        fontsize=14,
        ha="center",
        va="top",
        color=INK2,
    )
    top.text(
        1025,
        262,
        "small when every weight is near 0,\nso rr is near 1 everywhere",
        fontsize=14,
        ha="center",
        va="top",
        color=INK2,
    )
    top.text(
        1300,
        150,
        "λ = 0: no penalty\n(weights chase noise,\nrr too strong)",
        fontsize=14,
        va="top",
        color=INK,
    )
    top.text(
        1300,
        245,
        "λ very large: rr = 1\n(the signal is thrown\naway too)",
        fontsize=14,
        va="top",
        color=INK,
    )

    axes = [fig.add_axes([0.06 + i * 0.32, 0.1, 0.26, 0.45]) for i in range(3)]
    # (a) Weights as the penalty grows.
    ax = axes[0]
    piv = coefs.pivot(index="penalty", columns="column", values="weight")
    top_cols = piv.iloc[0].abs().sort_values(ascending=False).index[:4]
    for col in piv.columns:
        if col in top_cols:
            continue
        ax.plot(piv.index, piv[col], color=AXIS, lw=1.2)
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
    for col, c in zip(top_cols, palette):
        feat, scale = col.rsplit("_", 1)
        ax.plot(
            piv.index,
            piv[col],
            color=c,
            lw=2.4,
            label=f"{FEATURE_NAMES.get(feat, feat)}, {scale}",
        )
    ax.axhline(0, color=INK2, lw=0.8)
    ax.axvline(best, **THEORY)
    ax.set_xscale("log")
    ax.set_title(f"a  {context}: every weight shrinks toward 0")
    ax.set_xlabel("penalty λ (per DNM, log scale)")
    ax.set_ylabel("weight of each column (log odds per SD)")
    ax.legend(loc="lower right", frameon=False, fontsize=12.5)
    style_axis(ax, grid_axis="y")

    # (b) rr range as the penalty grows.
    ax = axes[1]
    ax.fill_between(
        rng_.penalty,
        rng_.rr_p05,
        rng_.rr_p95,
        color=CLASS_C["C non-CpG"],
        alpha=0.35,
        lw=0,
    )
    ax.plot(rng_.penalty, rng_.rr_p05, color=CLASS_C["C non-CpG"], lw=2)
    ax.plot(rng_.penalty, rng_.rr_p95, color=CLASS_C["C non-CpG"], lw=2)
    ax.axhline(1, color=INK2, lw=0.8)
    ax.axvline(best, **THEORY)
    ax.set_xscale("log")
    ax.set_title(f"b  {context}: so rr's range narrows")
    ax.set_xlabel("penalty λ (per DNM, log scale)")
    ax.set_ylabel("rr of the middle 90% of sites")
    style_axis(ax, grid_axis="y")

    # (c) Held-out DNM log loss.
    ax = axes[2]
    rel = loss.heldout_log_loss - loss.heldout_log_loss.min()
    ax.plot(loss.penalty, rel, color=REC_C, lw=2.6, marker="o", ms=8)
    ax.axvline(best, **THEORY)
    ax.text(
        best * 1.2, rel.max() * 0.3, f"chosen:\nλ = {best:g}", fontsize=14, color=INK
    )
    ax.set_xscale("log")
    ax.set_title("c  Chosen by predicting held-out DNMs")
    ax.set_xlabel("penalty λ (per DNM, log scale)")
    ax.set_ylabel("held-out log loss, all contexts\n(above the best; lower = better)")
    ax.text(
        0.97,
        0.97,
        "fit on odd chromosomes,\nscore on even, and the\nreverse; no gnomAD data",
        transform=ax.transAxes,
        va="top",
        ha="right",
        fontsize=13,
        color=INK2,
    )
    style_axis(ax, grid_axis="y")
    for ax in axes:
        ticks = [0.1, 1, 10, 100, 1000]
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax.minorticks_off()
    fig.text(
        0.06,
        0.015,
        "Local approximation (Satterstrom DNMs against the current background,"
        " built from the exported cells), with the penalty per DNM. The pipeline's default,"
        " 0.3 per site, was chosen the same way.",
        fontsize=13,
        color=INK2,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "rr_penalty.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_rates_diagram():
    """Draw what Roulette is, what we do with its rates, and the home-built version."""
    W, H = 1600, 760
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ax.text(
        30,
        34,
        "Per-variant rates: what Roulette is, and how we tried it (a quick first try)",
        fontsize=16,
        weight="bold",
        va="center",
    )
    cw = 470
    xs = [30, 565, 1100]
    heads = [
        "1  Roulette (Seplyarskiy et al. 2023)",
        "2  What the pipeline does with the rates",
        "3  The home-built version",
    ]
    for x, h in zip(xs, heads):
        ax.text(x, 80, h, fontsize=14, weight="bold", va="center")
    sz = 12

    # 1. Roulette.
    x = xs[0]
    box(
        ax,
        x,
        105,
        cw,
        92,
        RATE_FILL,
        RATE_EDGE,
        [
            "Learned from rare variants",
            "in gnomAD v3 genomes: where",
            "mutations happened in recent history",
        ],
        size=sz,
    )
    arrow(ax, x + cw / 2, 199, x + cw / 2, 223)
    box(
        ax,
        x,
        225,
        cw,
        172,
        RATE_FILL,
        RATE_EDGE,
        [
            "A model of each possible SNV's rate",
            "• the wider sequence around the site",
            "   (beyond our three bases)",
            "• CpG methylation and transcription",
            "• an empirical regional rate: each",
            "   trinucleotide's mutability in 50 kb",
        ],
        size=sz,
    )
    arrow(ax, x + cw / 2, 399, x + cw / 2, 423)
    box(
        ax,
        x,
        425,
        cw,
        72,
        RATE_FILL,
        RATE_EDGE,
        ["A rate for every possible SNV", "(v5.2; low-quality sites are flagged)"],
        size=sz,
    )
    arrow(ax, x + cw + 4, 461, xs[1] - 6, 168, color=RATE_EDGE, lw=1.6)

    # 2. What the pipeline does.
    x = xs[1]
    box(
        ax,
        x,
        105,
        cw,
        72,
        SAME_FILL,
        SAME_EDGE,
        ["Import the rates of coding SNVs", "--import-variant-rates"],
        size=sz,
    )
    arrow(ax, x + cw / 2, 179, x + cw / 2, 203)
    box(
        ax,
        x,
        205,
        cw,
        120,
        SAME_FILL,
        SAME_EDGE,
        [
            "Turn each rate into a multiplier",
            "rate ÷ the average rate of its group",
            "(context × ref × alt × methylation), so it",
            "adds only what our μ leaves out; CpG = 1",
        ],
        size=sz,
    )
    arrow(ax, x + cw * 0.27, 327, x + cw * 0.27, 371)
    arrow(ax, x + cw * 0.75, 327, x + cw * 0.75, 371)
    box(
        ax,
        x,
        373,
        cw * 0.5 - 6,
        150,
        RR_FILL,
        RR_EDGE,
        [
            "Train the rr",
            "log(multiplier) as a",
            "fixed offset: the ridge",
            "learns only what the",
            "rates miss",
        ],
        size=sz,
    )
    box(
        ax,
        x + cw * 0.5 + 6,
        373,
        cw * 0.5 - 6,
        150,
        "white",
        ROULETTE_C,
        [
            "Expected counts",
            "calibrated expected",
            "× multiplier × rr",
            "(rr normalized with the",
            "multiplier included)",
        ],
        size=sz,
        lw=1.8,
    )
    arrow(ax, x + cw * 0.5 - 4, 448, x + cw * 0.5 + 4, 448)
    ax.text(
        x, 545, "--use-variant-rates, --variant-rate-ht-path", fontsize=11, color=INK2
    )

    # 3. Home-built.
    x = xs[2]
    box(
        ax,
        x,
        105,
        cw,
        92,
        SAME_FILL,
        SAME_EDGE,
        [
            "Rare variants in gnomAD v4.1 genomes",
            "outside exons ± 50 bp, low-complexity",
            "regions and segdups",
        ],
        size=sz,
    )
    arrow(ax, x + cw * 0.27, 199, x + cw * 0.27, 223)
    arrow(ax, x + cw * 0.75, 199, x + cw * 0.75, 223)
    box(
        ax,
        x,
        225,
        cw * 0.5 - 6,
        190,
        SAME_FILL,
        SAME_EDGE,
        [
            "Regional part",
            "each trinucleotide's",
            "observed share in",
            "50 kb windows ÷ its",
            "genome-wide share,",
            "shrunk; capped at 3",
        ],
        size=sz,
    )
    box(
        ax,
        x + cw * 0.5 + 6,
        225,
        cw * 0.5 - 6,
        190,
        SAME_FILL,
        SAME_EDGE,
        [
            "Sequence part",
            "each heptamer's",
            "genome-wide share ÷",
            "its trinucleotide's",
        ],
        size=sz,
    )
    arrow(ax, x + cw / 2, 417, x + cw / 2, 441)
    box(
        ax,
        x,
        443,
        cw,
        72,
        SAME_FILL,
        SAME_EDGE,
        ["multiplier = regional × sequence", "then used exactly like Roulette's (2)"],
        size=sz,
    )
    arrow(ax, x - 6, 479, xs[1] + cw + 6, 265, color=MUTED, lw=1.4)

    # What it buys and costs.
    box(
        ax,
        30,
        590,
        1540,
        140,
        "white",
        AXIS,
        [
            "A quick first try, not a worked-out method: one untuned run of each (see How each version performs)",
            "+ about twice the noise reduction of the ridge rr alone, half the unflagged outliers, and the maternal hotspot regions corrected",
            "+ held-out DNMs back it, independently of gnomAD",
            "− with the ridge on top it is too strong; it moves LoF (−5%) and missense (+3%) expected counts relative to synonymous;",
            "   part of its gain is gnomAD fitting gnomAD; the home-built regional part varies too much",
        ],
        size=12.5,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "rr_rates_diagram.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_outliers():
    """Draw synonymous outliers by category, and histone gene z-scores below."""
    out = read_outliers()
    rows = [
        ("No correction", "none", NONE_C, None),
        ("Gnocchi adj_r", "adj_r", ADJ_R_C, None),
        ("Legacy syn_adj_r", "syn_adj_r", LEGACY_C, None),
        ("Ruchit, 17k DNMs", "ruchit_17k", RUCHIT_C, None),
        ("Ruchit, Satterstrom\nDNMs", "ruchit_satterstrom", RUCHIT_C, "///"),
        ("Ridge, Satterstrom\nDNMs", "recommended", REC_C, None),
    ]
    flag_fill = {
        "histone genes": "#cfc8f2",
        "low_exome_coverage": "#e6e5df",
        "low_exome_mapping_quality": "#c3c2b7",
    }
    fig = plt.figure(figsize=(15, 12.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], hspace=0.5, wspace=0.22)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.1)
    ax = fig.add_subplot(gs[0, 0])
    ypos = list(range(len(rows)))[::-1]
    for y, (_, corr, c, h) in zip(ypos, rows):
        for direction, sign in (("deficit", -1), ("excess", 1)):
            start = 0
            for flag in ("unflagged", *flag_fill):
                n = out_count(out, corr, direction, flag)
                if not n:
                    continue
                unflagged = flag == "unflagged"
                ax.barh(
                    y,
                    sign * n,
                    left=sign * start,
                    height=0.56,
                    color=c if unflagged else flag_fill[flag],
                    hatch=h if unflagged else None,
                    hatchcolor="white",
                    edgecolor=SURFACE,
                    linewidth=1.5,
                )
                if flag == "histone genes" and n >= 15:
                    ax.text(
                        sign * (start + n / 2),
                        y,
                        f"{n}",
                        ha="center",
                        va="center",
                        fontsize=13,
                        color=INK,
                    )
                if unflagged and n >= 25:
                    ax.text(
                        sign * (start + n / 2),
                        y,
                        f"{n}",
                        ha="center",
                        va="center",
                        fontsize=14,
                        fontweight="bold",
                        color="white",
                        bbox=dict(boxstyle="round,pad=0.12", facecolor=c, linewidth=0),
                    )
                start += n
            ax.text(
                sign * (start + 8),
                y,
                f"{start}",
                ha="left" if sign > 0 else "right",
                va="center",
                fontsize=14,
                color=INK,
            )
    ax.axvline(0, color=INK2, lw=1)
    ax.text(
        -10,
        len(rows) - 0.38,
        "← obs ≪ exp (deficit)",
        ha="right",
        fontsize=14,
        color=INK2,
    )
    ax.text(
        10,
        len(rows) - 0.38,
        "obs ≫ exp (excess) →",
        ha="left",
        fontsize=14,
        color=INK2,
    )
    ax.text(175, -1.0, "chance alone: ≈ 0", ha="right", fontsize=13, color=INK2)
    ax.set_yticks(ypos, [r[0] for r in rows])
    ax.set_xlim(-440, 180)
    ax.set_ylim(-1.3, len(rows) - 0.05)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    ax.set_title("a  Synonymous outliers (|z| > 8), histone genes split out")
    ax.set_xlabel("transcripts (of 18,406)")
    ax.legend(
        handles=[
            Patch(color=INK2, label="none of these (bar color = version)"),
            Patch(color=flag_fill["histone genes"], label="histone genes"),
            Patch(
                color=flag_fill["low_exome_coverage"], label="low_exome_coverage flag"
            ),
            Patch(
                color=flag_fill["low_exome_mapping_quality"],
                label="low_exome_mapping_quality flag",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.42, -0.13),
        ncol=2,
        frameon=False,
        fontsize=13,
    )
    style_axis(ax)

    # (b) Synonymous z of histone genes vs every other gene, no correction.
    ax = fig.add_subplot(gs[0, 1])
    zs = pd.read_csv(os.path.join(OUT_DIR, "syn_z_by_histone_group.tsv"), sep="\t")
    edges = np.arange(-4, 19, 1.0)
    centers = edges[:-1] + 0.5
    groups = [
        ("histone, chr6p22 cluster", "histone, chr6p22 cluster", "#8f84e0", -0.3),
        ("histone, elsewhere", "histone, elsewhere", "#cfc8f2", 0.0),
        (
            "non-histone, same chr6 region",
            "non-histone, same chr6 region",
            "#5f5e5a",
            0.3,
        ),
    ]
    for group, label, c, off in groups:
        z = zs[zs.group == group].z_none.clip(edges[0], edges[-1] - 1e-9)
        pct = 100 * np.histogram(z, edges)[0] / len(z)
        ax.bar(
            centers + off,
            pct,
            width=0.29,
            color=c,
            linewidth=0,
            label=f"{label} (n = {len(z)})",
        )
    other = zs[zs.group == "other genes"].z_none.clip(edges[0], edges[-1] - 1e-9)
    pct = 100 * np.histogram(other, edges)[0] / len(other)
    ax.step(
        np.append(edges[:-1], edges[-1]),
        np.append(pct, pct[-1]),
        where="post",
        color=NONE_C,
        lw=1.8,
        label=f"all other genes (n = {len(other):,})",
    )
    ax.axvline(8, **THEORY)
    ax.text(8.15, 33, "outlier cutoff\n(z = 8)", fontsize=13, color=INK2, va="top")
    ax.axvline(0, color=AXIS, lw=1)
    ax.set_xlim(-4, 18)
    ax.set_ylim(0, 60)
    ax.set_xticks(range(-4, 19, 2))
    ax.set_title("b  Synonymous z of histone genes (no correction)")
    ax.set_xlabel(
        "synonymous z = (obs − exp) / √exp\npositive = more variants than expected"
    )
    ax.set_ylabel("% of genes in the group")
    leg = ax.legend(
        loc="upper right", fontsize=12.5, bbox_to_anchor=(1.0, 1.0), framealpha=1
    )
    leg.get_frame().set_facecolor(SURFACE)
    leg.get_frame().set_linewidth(0)
    style_axis(ax, grid_axis="y")

    # (c, d) Gene size: a bigger gene crosses |z| > 8 with a smaller relative excess.
    # c uses expected count (what z depends on); d uses coding length.
    zs["oe"] = zs.obs / zs.exp_recommended
    other = zs[zs.group.isin(["other genes", "non-histone, same chr6 region"])]
    hit = other[other.z_recommended.abs() > 8]
    exp_per_bp = float((other.exp_recommended / other.cds_length).median())
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1], sharey=ax_c)
    for ax, xcol, title in (
        (
            ax_c,
            "exp_recommended",
            "c  Gene size by expected count: bigger genes need a smaller excess",
        ),
        (ax_d, "cds_length", "d  Gene size by coding length"),
    ):
        ax.scatter(
            other[xcol],
            other.oe,
            s=4,
            color=NONE_C,
            alpha=0.25,
            linewidth=0,
            label=f"genes without a QC flag (n = {len(other):,})",
        )
        ax.scatter(
            hit[xcol],
            hit.oe,
            s=38,
            color=REC_C,
            edgecolor=SURFACE,
            linewidth=0.8,
            zorder=3,
            label=f"non-histone outliers (n = {len(hit)})",
        )
        for group, label, c in (
            ("histone, chr6p22 cluster", "histone genes, chr6p22 cluster", "#8f84e0"),
            ("histone, elsewhere", "histone genes elsewhere", "#cfc8f2"),
        ):
            h = zs[zs.group == group]
            ax.scatter(
                h[xcol],
                h.oe,
                s=38,
                color=c,
                edgecolor=SURFACE,
                linewidth=0.8,
                zorder=3,
                label=f"{label} (n = {len(h)})",
            )
        ax.set_xscale("log")
        if xcol == "exp_recommended":
            x = np.logspace(np.log10(5), np.log10(6000), 200)
            e = x
            ax.set_xlim(5, 6000)
            ax.set_xticks([10, 100, 1000], ["10", "100", "1,000"])
            ax.set_xlabel("expected synonymous variants in the gene")
            curve_label, label_x = "|z| = 8", 26
        else:
            x = np.logspace(np.log10(30), np.log10(40000), 200)
            e = x * exp_per_bp
            ax.set_xlim(30, 40000)
            ax.set_xticks([100, 1000, 10000], ["100 bp", "1 kb", "10 kb"])
            ax.set_xlabel("coding sequence length")
            curve_label, label_x = "|z| = 8 (at the median\n160 expected per kb)", 170
        ax.plot(x, 1 + 8 / np.sqrt(e), **THEORY)
        ax.plot(x, 1 - 8 / np.sqrt(e), **THEORY)
        ax.text(label_x, 2.95, curve_label, fontsize=12.5, color=INK2, va="top")
        for gene, dx, dy in (("CSMD1", -40, 10), ("H1-4", 8, 6), ("MUC4", 8, -12)):
            r = zs[zs.gene == gene]
            if len(r):
                r = r.iloc[0]
                ax.annotate(
                    gene,
                    (r[xcol], r.oe),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=12.5,
                    color=INK2,
                )
        ax.set_title(title)
        style_axis(ax, grid_axis="y")
    ax_c.set_ylim(0.3, 3.0)
    ax_c.set_ylabel("synonymous obs/exp (ridge rr)")
    plt.setp(ax_d.get_yticklabels(), visible=False)
    handles, labels = ax_c.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=13,
        bbox_to_anchor=(0.54, 0.0),
        markerscale=1.4,
    )

    fig.savefig(
        os.path.join(OUT_DIR, "rr_outliers.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_alternatives():
    """Draw out-of-fold noise and calibration slope for the rr alternatives."""
    res = pd.concat(
        [
            pd.read_csv(os.path.join(OUT_DIR, f), sep="\t")
            for f in (
                "rr_shrinkage_experiments.tsv",
                "rr_shrinkage_experiments.pass2.tsv",
            )
        ]
    ).drop_duplicates(["variant", "dnm_set"], keep="last")
    res["dnm_set"] = res.dnm_set.fillna("")

    def get(variant, dnm_set=""):
        return res[(res.variant == variant) & (res.dnm_set == dnm_set)].iloc[0]

    none = get("no correction")
    rows = [
        ("Now: Bonferroni, no penalty", "Bonferroni, unpenalized (his method)", "17k"),
        ("", "Bonferroni, unpenalized (his method)", "Satterstrom"),
        ("Ridge, tuned on held-out DNMs", "ridge, all features, C=0.0005", "17k"),
        ("", "ridge, all features, C=0.0005", "Satterstrom"),
        ("Sample split", "sample split: screen half A, fit half B", "17k"),
        ("", "sample split: screen half A, fit half B", "Satterstrom"),
        (
            "Background without gnomAD sites*",
            "background: absent from gnomAD (current design)",
            "17k",
        ),
        ("Background with all sites*", "background: all possible sites", "17k"),
        ("Fit on gnomAD synonymous", "no DNMs: Poisson on gnomAD syn, alpha=0.1", ""),
    ]
    style = {
        "17k": (RUCHIT_C, "o"),
        "Satterstrom": (RUCHIT_C, "s"),
        "": ("#6250d6", "D"),
    }
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 5.6), gridspec_kw={"width_ratios": [1.6, 1]}, sharey=True
    )
    fig.subplots_adjust(left=0.3, right=0.97, top=0.88, bottom=0.2, wspace=0.08)
    ypos = list(range(len(rows)))[::-1]
    labels = []
    for y, (label, variant, dnm_set) in zip(ypos, rows):
        r = get(variant, dnm_set)
        is_ridge_satt = (
            variant == "ridge, all features, C=0.0005" and dnm_set == "Satterstrom"
        )
        c, m = (REC_C, "D") if is_ridge_satt else style[dnm_set]
        set_label = {"17k": "17k", "Satterstrom": "Satterstrom", "": "no DNMs"}[dnm_set]
        labels.append(f"{label}  ·  {set_label}" if label else set_label)
        ax = axes[0]
        ax.plot([r.core_noise, r.crossfit_noise], [y, y], color=c, lw=1.5, zorder=2)
        ax.scatter(
            [r.core_noise],
            [y],
            s=112,
            marker=m,
            color=c,
            edgecolor=SURFACE,
            linewidth=1.5,
            zorder=3,
        )
        ax.scatter(
            [r.crossfit_noise],
            [y],
            s=96,
            marker=m,
            color=SURFACE,
            edgecolor=c,
            linewidth=1.8,
            zorder=3,
        )
        axes[1].scatter(
            [r.slope],
            [y],
            s=112,
            marker=m,
            color=c,
            edgecolor=SURFACE,
            linewidth=1.5,
            zorder=3,
        )
        axes[1].text(
            r.slope + (0.05 if r.slope > 1 else -0.05),
            y,
            f"{r.slope:.2f}",
            ha="left" if r.slope > 1 else "right",
            va="center",
            fontsize=13,
            color=INK2,
        )
    ax = axes[0]
    ax.axvline(none.core_noise, **THEORY)
    ax.text(
        none.core_noise - 0.01,
        len(rows) - 0.35,
        "no correction",
        ha="right",
        fontsize=13,
        color=INK2,
    )
    ax.set_yticks(ypos, labels)
    ax.set_ylim(-0.7, len(rows) - 0.1)
    ax.set_xlim(2.3, 2.95)
    ax.set_title("Noise on held-out chromosomes")
    ax.set_xlabel(
        "mean z², overall level removed (lower = better)\nfilled = rr as fitted; open = with the best α added"
    )
    style_axis(ax)
    ax = axes[1]
    ax.axvline(1.0, **THEORY)
    ax.text(1.02, len(rows) - 0.35, "1 = no α needed", fontsize=13, color=INK2)
    ax.set_xlim(0.3, 1.6)
    ax.set_title("Calibration slope")
    ax.set_xlabel("1 = right strength")
    style_axis(ax)
    fig.text(
        0.01,
        0.02,
        "* background built from the exported MANE Select cells (approximation); compare the two rows with each other only",
        fontsize=12.5,
        color=INK2,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "rr_alternatives.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_leftover_noise():
    """Draw where the noise left after the ridge rr sits, from residual_noise.py."""
    size = pd.read_csv(os.path.join(OUT_DIR, "residual_noise.by_size.tsv"), sep="\t")
    cov = pd.read_csv(os.path.join(OUT_DIR, "residual_noise.covariates.tsv"), sep="\t")
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1, 1.25]}
    )
    fig.subplots_adjust(left=0.07, right=0.97, top=0.86, bottom=0.17, wspace=0.62)

    ax = axes[0]
    x = size.median_exp
    ax.plot(x, size.chance, label="Chance alone", **THEORY)
    ax.plot(x, size.none, color=NONE_C, lw=2, marker="o", ms=7, label="No correction")
    ax.plot(x, size.ridge_rr, color=REC_C, lw=2.6, marker="D", ms=7, label="Ridge rr")
    ax.set_xscale("log")
    ax.set_xticks(
        [40, 80, 150, 300, 500, 1000], ["40", "80", "150", "300", "500", "1,000"]
    )
    ax.set_ylim(0, 6.2)
    ax.set_title("Leftover noise grows with gene size")
    ax.set_xlabel("expected synonymous variants per gene (bin median)")
    ax.set_ylabel("mean z²")
    ax.legend(loc="upper left", frameon=False, fontsize=13)
    style_axis(ax, grid_axis="y")

    ax = axes[1]
    order = [
        ("no correction", "No correction", NONE_C),
        ("none (ridge rr only)", "Ridge rr", REC_C),
        ("gene size", "+ gene size", "#6250d6"),
        ("expression (GTEx)", "+ expression (GTEx)", "#6250d6"),
        ("gene quality metrics", "+ gene quality metrics", "#6250d6"),
        (
            "13 regional features, nonlinear",
            "+ regional features, nonlinear",
            "#6250d6",
        ),
        ("all of the above", "+ all of the above", "#6250d6"),
    ]
    ypos = list(range(len(order)))[::-1]
    for y, (key, label, c) in zip(ypos, order):
        v = float(cov[cov.covariates == key].noise.iloc[0])
        ax.barh(y, v, height=0.56, color=c, linewidth=0)
        ax.text(v + 0.03, y, f"{v:.2f}", va="center", fontsize=14)
    ax.axvline(CHANCE_FLOOR, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.04,
        len(order) - 0.45,
        "chance ≤ 0.70",
        fontsize=13,
        color=INK2,
    )
    ax.set_yticks(ypos, [o[1] for o in order])
    ax.set_xlim(0, 3.0)
    ax.set_ylim(-0.6, len(order) - 0.1)
    ax.set_title("What gene-level covariates remove")
    ax.set_xlabel(
        "mean z², overall level removed (lower = better)\ncovariates fit on the other chromosome half"
    )
    style_axis(ax)

    fig.savefig(
        os.path.join(OUT_DIR, "rr_leftover_noise.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_outlier_drivers():
    """Draw the DNM, regional and downsampling checks of the outliers, from outlier_drivers.py."""
    groups = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.groups.tsv"), sep="\t"
    ).set_index("set")
    dec = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.dnm_by_oe_decile.tsv"), sep="\t"
    )
    regions = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.regions.tsv"), sep="\t"
    )
    ds = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.downsampling.tsv"), sep="\t"
    )
    genes = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.all_genes.tsv"), sep="\t"
    )
    hot_c, other_c, repeat_c, histone_c = "#d9622b", "#1a9e6e", "#a8a7a0", "#8f84e0"
    group_style = {
        "maternal hotspot region": (hot_c, "in a maternal\nhotspot region"),
        "other": (other_c, "other"),
        "repeat, segdup or low coverage": (repeat_c, "repeat, segdup or\nlow coverage"),
    }

    fig = plt.figure(figsize=(14, 15))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.78, 1], hspace=0.5, wspace=0.25)
    fig.subplots_adjust(left=0.1, right=0.97, top=0.96, bottom=0.06)

    # (a) DNMs vs gnomAD for the outlier groups.
    ax = fig.add_subplot(gs[0, 0])
    sets = [
        (
            "maternal hotspot region",
            hot_c,
            "excess outliers in a\nmaternal hotspot region",
            None,
        ),
        ("other", other_c, "other excess\noutliers", (-12, 0)),
        (
            "repeat, segdup or low coverage",
            repeat_c,
            "repeat, segdup or\nlow-coverage outliers",
            (12, 12),
        ),
        (
            "histone genes, chr6p22 cluster",
            histone_c,
            "histone genes,\nchr6p22 cluster",
            (-12, 0),
        ),
        (
            "size-matched controls",
            INK2,
            "size-matched genes\nwithout a flag",
            (10, -14),
        ),
    ]
    x = np.array([0.9, 2.0])
    ax.plot(x, x, **THEORY)
    ax.text(1.93, 1.83, "DNMs agree", fontsize=13, color=INK2, ha="right", rotation=36)
    ax.axhline(1, color=AXIS, lw=1)
    ax.text(
        1.98,
        0.95,
        "DNMs show no excess",
        fontsize=13,
        color=INK2,
        ha="right",
        va="top",
    )
    for key, c, label, offset_pt in sets:
        dx, dy = offset_pt or (None, None)
        r = groups.loc[key]
        ax.errorbar(
            r.gnomad_oe,
            r.dnm_oe,
            yerr=[[r.dnm_oe - r.dnm_lo], [r.dnm_hi - r.dnm_oe]],
            fmt="o",
            ms=12,
            color=c,
            mec=SURFACE,
            mew=1.5,
            elinewidth=2,
            capsize=0,
            zorder=3,
        )
        n = f"n = {int(r.n):,}"
        if (dx, dy) == (None, None):
            ax.annotate(
                f"{label} ({n})",
                (r.gnomad_oe, r.dnm_oe),
                xytext=(1.02, 2.3),
                fontsize=13,
                color=INK,
                va="center",
                arrowprops=dict(arrowstyle="-", color=c, lw=1, shrinkB=6),
            )
            continue
        ax.annotate(
            f"{label} ({n})",
            (r.gnomad_oe, r.dnm_oe),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=13,
            color=INK,
            va="center",
            ha="right" if dx < 0 else "left",
        )
    ax.set_xlim(0.9, 2.0)
    ax.set_ylim(0.7, 3.7)
    ax.set_title("a  DNMs show the same excess, except in repeats")
    ax.set_xlabel("gnomAD synonymous obs/exp (ridge rr)")
    ax.set_ylabel("synonymous DNM obs/exp (ridge rr), 95% CI")
    style_axis(ax, grid_axis="y")

    # (b) DNMs by gnomAD obs/exp decile, all genes.
    ax = fig.add_subplot(gs[0, 1])
    x = np.array([0.82, 1.24])
    ax.axhline(1, color=AXIS, lw=1)
    ax.plot(x, 1 + 0.75 * (x - 1), **THEORY)
    ax.text(
        1.15,
        1.13,
        "if all the spread above\nchance were mutation rate",
        fontsize=13,
        color=INK2,
        ha="right",
        va="bottom",
    )
    ax.text(
        1.235,
        0.975,
        "if none of it were",
        fontsize=13,
        color=INK2,
        ha="right",
        va="top",
    )
    ax.plot(x, 1 + 0.51 * (x - 1), color=REC_C, lw=1.6)
    ax.text(
        0.9,
        0.885,
        "solid line: fit, slope 0.51 ± 0.07",
        fontsize=13,
        color=INK,
        va="top",
    )
    ax.errorbar(
        dec.gnomad_oe,
        dec.dnm_oe,
        yerr=[dec.dnm_oe - dec.dnm_lo, dec.dnm_hi - dec.dnm_oe],
        fmt="o",
        ms=9.5,
        color=REC_C,
        mec=SURFACE,
        mew=1.2,
        elinewidth=1.6,
        capsize=0,
        zorder=3,
    )
    ax.set_xlim(0.82, 1.24)
    ax.set_ylim(0.82, 1.2)
    ax.set_title("b  DNMs follow synonymous obs/exp")
    ax.set_xlabel(
        "gnomAD synonymous obs/exp (ridge rr), decile mean\n15,027 genes without a flag, ≥ 100 expected"
    )
    ax.set_ylabel("synonymous DNM obs/exp, 95% CI")
    style_axis(ax, grid_axis="y")

    # (c) Where the outliers are: genome-wide z with the maternal hotspot regions shaded.
    ax = fig.add_subplot(gs[1, :])
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    genes = genes[genes.chromosome.isin(chroms)]
    length = genes.groupby("chromosome").mid.max().reindex(chroms).fillna(0) + 2e6
    offset = pd.Series(
        np.concatenate([[0], np.cumsum(length.values[:-1])]), index=chroms
    )
    genes["x"] = genes.mid + genes.chromosome.map(offset)
    hotspots = [
        ("chr8", 0, 20e6),
        ("chr2", 31.8e6, 61.0e6),
        ("chr7", 0, 60.1e6),
        ("chr9", 0, 43.0e6),
        ("chr16", 0, 90.4e6),
    ]
    for c, s, e in hotspots:
        ax.axvspan(
            offset[c] + s,
            offset[c] + min(e, length[c]),
            color="#fdf0cc",
            lw=0,
            zorder=0,
        )
    plain = genes[~genes.outlier & ~genes.histone_6p22]
    ax.scatter(plain.x, plain.z.clip(-12, 24), s=3, color=NONE_C, alpha=0.3, lw=0)
    h = genes[genes.histone_6p22]
    ax.scatter(h.x, h.z.clip(-12, 24), s=26, color=histone_c, lw=0, zorder=2)
    for key, (c, _) in group_style.items():
        g = genes[genes.outlier & (genes.group == key)]
        ax.scatter(
            g.x, g.z.clip(-12, 24), s=48, color=c, edgecolor=SURFACE, lw=0.7, zorder=3
        )
    ax.axhline(8, **THEORY)
    ax.axhline(-8, **THEORY)
    ax.text(offset["chr4"] + 4e7, 8.6, "|z| = 8", fontsize=13, color=INK2)
    for gene, dx, dy in (
        ("CSMD1", 6, 2),
        ("PIEZO1", 6, 2),
        ("MUC5AC", 6, 2),
        ("FLG", 6, 2),
        ("SRRM2", -6, 2),
    ):
        r = genes[genes.gene == gene]
        if len(r):
            r = r.iloc[0]
            ax.annotate(
                gene,
                (r.x, min(r.z, 24)),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=12.5,
                color=INK2,
                ha="right" if dx < 0 else "left",
            )
    centers = offset + length / 2
    ax.set_xticks(centers.values, [c.replace("chr", "") for c in chroms], fontsize=12.5)
    ax.set_xlim(0, offset["chrY"] + length["chrY"])
    ax.set_ylim(-12, 25)
    ax.set_ylabel("synonymous z (ridge rr)")
    ax.set_title(
        "c  23 of the 46 excess outliers sit in the maternal hotspot regions (shaded); 3.9 expected"
    )
    ax.legend(
        handles=[
            Patch(
                color="#fdf0cc",
                label="maternal C>G hotspot region (Jónsson et al. 2017)",
            ),
            plt.Line2D(
                [],
                [],
                marker="o",
                ls="",
                color=hot_c,
                label="outlier in a hotspot region",
            ),
            plt.Line2D([], [], marker="o", ls="", color=other_c, label="other outlier"),
            plt.Line2D(
                [],
                [],
                marker="o",
                ls="",
                color=repeat_c,
                label="repeat, segdup or low-coverage outlier",
            ),
            plt.Line2D(
                [],
                [],
                marker="o",
                ls="",
                color=histone_c,
                label="histone gene, chr6p22 cluster",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=3,
        frameon=False,
        fontsize=13,
    )
    style_axis(ax, grid_axis="y")

    # (d) Ordinary genes in each hotspot region: the rr corrects only part of it.
    ax = fig.add_subplot(gs[2, 0])
    order = [
        "8p, first 20 Mb",
        "2p16-22",
        "7p",
        "9p",
        "16p",
        "16q",
        "rest of the genome",
    ]
    reg = regions.set_index("set").loc[order]
    names = ["8p, 0–20 Mb", "2p16–22", "7p", "9p", "16p", "16q", "elsewhere"]
    labels = [f"{n}\n{int(k):,} genes" for n, k in zip(names, reg.n)]
    ypos = np.arange(len(order))[::-1]
    for y, (_, r) in zip(ypos, reg.iterrows()):
        ax.plot([r.dnm_lo, r.dnm_hi], [y - 0.18] * 2, color=INK2, lw=1.6)
        ax.scatter(
            r.dnm_oe, y - 0.18, s=88, color=SURFACE, edgecolor=INK2, lw=1.6, zorder=3
        )
        ax.scatter(
            r.gnomad_oe_no_rr,
            y + 0.14,
            s=88,
            color=NONE_C,
            edgecolor=SURFACE,
            lw=1,
            zorder=3,
        )
        ax.scatter(
            r.gnomad_oe,
            y + 0.14,
            s=96,
            marker="D",
            color=REC_C,
            edgecolor=SURFACE,
            lw=1,
            zorder=3,
        )
    ax.axvline(1, color=AXIS, lw=1)
    ax.set_yticks(ypos, labels)
    ax.set_xlim(0.85, 1.62)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title("d  Genes in the hotspot regions that aren't outliers")
    ax.set_xlabel("synonymous obs/exp (outliers left out)")
    ax.legend(
        handles=[
            plt.Line2D([], [], marker="o", ls="", color=NONE_C, label="gnomAD, no rr"),
            plt.Line2D(
                [], [], marker="D", ls="", color=REC_C, label="gnomAD, ridge rr"
            ),
            plt.Line2D(
                [],
                [],
                marker="o",
                ls="-",
                color=INK2,
                mfc=SURFACE,
                label="DNMs, ridge rr (95% CI)",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.45, -0.13),
        ncol=3,
        frameon=False,
        fontsize=13,
    )
    style_axis(ax)

    # (e) The excess by downsampling and genetic ancestry group.
    ax = fig.add_subplot(gs[2, 1])
    d = groups.loc["all excess outliers"]
    ax.axhspan(d.dnm_lo, d.dnm_hi, color="#f0efec", lw=0)
    ax.text(
        1100,
        d.dnm_lo + 0.012,
        "DNM obs/exp of the same genes (95% CI)",
        fontsize=13,
        color=INK2,
        va="bottom",
    )
    ax.axhline(1, color=AXIS, lw=1)
    for anc, g in ds.groupby("ancestry"):
        g = g[g.n_samples >= 1000]
        if anc == "global":
            ax.plot(
                g.n_samples,
                g.rel_oe,
                color=REC_C,
                lw=2.6,
                marker="o",
                ms=7,
                zorder=3,
                label="all samples",
            )
        else:
            ax.plot(
                g.n_samples,
                g.rel_oe,
                color=NONE_C,
                lw=1.3,
                marker="o",
                ms=4,
                label=(
                    "one ancestry group each: African, Admixed American,\nEast Asian, European (non-Finnish), South Asian"
                    if anc == "afr"
                    else None
                ),
            )
    ax.legend(loc="upper right", frameon=False, fontsize=13)
    ax.set_xscale("log")
    ax.set_xticks([1000, 10000, 100000, 500000], ["1k", "10k", "100k", "500k"])
    ax.set_xlim(800, 8e5)
    ax.set_ylim(0.9, 1.8)
    ax.set_title("e  Same excess in every ancestry and size")
    ax.set_xlabel("downsampled exomes (log scale)")
    ax.set_ylabel(
        "excess outliers' synonymous obs/exp,\nrelative to genes without a flag"
    )
    style_axis(ax, grid_axis="y")

    fig.savefig(
        os.path.join(OUT_DIR, "rr_outlier_drivers.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_features():
    """Draw the ridge coefficients of each regional feature, from rr_background_experiments.py."""
    coef = pd.read_csv(
        os.path.join(OUT_DIR, "rr_background.coefficients.tsv"), sep="\t"
    )
    coef = coef[coef.background == "absent"]
    features = [
        ("GC_content", "GC content"),
        ("CpG_island", "CpG islands"),
        ("met_sperm", "sperm methylation"),
        ("Nucleosome", "nucleosome occupancy"),
        ("SINE", "SINEs"),
        ("LINE", "LINEs"),
        ("LCR", "low-complexity regions"),
        ("recomb_male", "male recombination"),
        ("recomb_female", "female recombination"),
        ("dist2telo", "distance to telomere"),
        ("dist2cent", "distance to centromere"),
        ("cDNM_maternal_05M", "maternal clustered DNMs"),
        ("cDNM_paternal_05M", "paternal clustered DNMs"),
    ]
    scales = ["1k", "10k", "100k", "1M"]
    classes = ["A/T", "C non-CpG", "CpG"]
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#e1e0d9", "#d9622b"])
    # Percent change in the DNM rate per SD of the feature, averaged over contexts.
    mean = 100 * (np.exp(coef.groupby(["class", "feature", "scale"]).coef.mean()) - 1)
    total = mean.groupby(["class", "feature"]).sum()
    vmax = float(np.nanpercentile(np.abs(total.to_numpy()), 95))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 6.2), sharey=True)
    fig.subplots_adjust(left=0.17, right=0.9, top=0.86, bottom=0.12, wspace=0.08)
    for ax, c in zip(axes, classes):
        m = np.full((len(features), len(scales) + 1), np.nan)
        for i, (f, _) in enumerate(features):
            for j, s in enumerate(scales):
                if (c, f, s) in mean.index:
                    m[i, j] = mean[(c, f, s)]
            if (c, f) in total.index and not np.isnan(m[i, :4]).all():
                m[i, 4] = total[(c, f)]
        im = ax.imshow(m, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        ax.axvline(3.5, color=SURFACE, lw=4)
        for i in range(len(features)):
            for j in range(len(scales) + 1):
                if np.isnan(m[i, j]):
                    ax.text(
                        j,
                        i,
                        "not used",
                        ha="center",
                        va="center",
                        fontsize=10.5,
                        color=MUTED,
                    )
                elif j == 4 or abs(m[i, j]) >= 0.5 * vmax:
                    ax.text(
                        j,
                        i,
                        f"{m[i, j]:+.1f}",
                        ha="center",
                        va="center",
                        fontsize=12,
                        color=INK,
                        fontweight="bold" if j == 4 else "normal",
                    )
        ax.set_xticks(
            range(len(scales) + 1), ["1 kb", "10 kb", "100 kb", "1 Mb", "all"]
        )
        ax.set_yticks(range(len(features)), [f[1] for f in features])
        ax.set_title(c, loc="center")
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(False)
    cax = fig.add_axes([0.915, 0.25, 0.012, 0.5])
    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_visible(False)
    cb.set_label(
        "% change in DNM rate per SD of the feature\n(red = more mutable)", fontsize=13
    )
    fig.suptitle(
        "What each feature does in the ridge rr: % change in rate per SD, by window scale ('all' = sum)",
        x=0.17,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.17,
        0.04,
        "Satterstrom DNMs against the absent-from-gnomAD background (exported MANE Select cells). "
        "Correlated features share signal, so read directions, not single values.",
        fontsize=13,
        color=INK2,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "rr_features.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_dnm_sets():
    """Draw held-out noise of the two DNM sets under each fit, from rr_shrinkage_experiments.py."""
    res = pd.concat(
        [
            pd.read_csv(os.path.join(OUT_DIR, f), sep="\t")
            for f in (
                "rr_shrinkage_experiments.tsv",
                "rr_shrinkage_experiments.pass2.tsv",
            )
        ]
    ).drop_duplicates(["variant", "dnm_set"], keep="last")
    rows = [
        (
            "Bonferroni screen, no penalty\n(current method)",
            "Bonferroni, unpenalized (his method)",
            "core_noise",
        ),
        (
            "Bonferroni screen, no penalty,\nwith the best α",
            "Bonferroni, unpenalized (his method)",
            "crossfit_noise",
        ),
        (
            "Screen on one half,\nfit on the other",
            "sample split: screen half A, fit half B",
            "core_noise",
        ),
        (
            "Ridge on all features,\nno screen",
            "ridge, all features, C=0.0005",
            "core_noise",
        ),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    fig.subplots_adjust(left=0.3, right=0.96, top=0.86, bottom=0.17)
    ypos = np.arange(len(rows))[::-1]
    for y, (label, variant, col) in zip(ypos, rows):
        v = {
            s: float(res[(res.variant == variant) & (res.dnm_set == s)][col].iloc[0])
            for s in ("17k", "Satterstrom")
        }
        ax.plot([v["17k"], v["Satterstrom"]], [y, y], color=RUCHIT_C, lw=1.5, zorder=2)
        ax.scatter(
            v["17k"],
            y,
            s=112,
            marker="o",
            color=RUCHIT_C,
            edgecolor=SURFACE,
            lw=1.4,
            zorder=3,
            label="17k DNMs" if y == ypos[0] else None,
        )
        ax.scatter(
            v["Satterstrom"],
            y,
            s=112,
            marker="s",
            color=SURFACE,
            edgecolor=RUCHIT_C,
            lw=1.8,
            zorder=3,
            label="Satterstrom DNMs" if y == ypos[0] else None,
        )
        gap = v["Satterstrom"] - v["17k"]
        ax.text(
            max(v.values()) + 0.012,
            y,
            f"gap {gap:+.3f}",
            va="center",
            fontsize=14,
            color=INK2,
        )
    ax.set_yticks(ypos, [r[0] for r in rows])
    ax.set_xlim(2.4, 2.86)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.45, 1.16),
        ncol=2,
        frameon=False,
        fontsize=14,
    )
    ax.set_xlabel(
        "per-gene synonymous noise on held-out chromosomes, mean z² (lower = better)"
    )
    style_axis(ax)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_dnm_sets.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_histones():
    """Draw the histone cluster checks, from histone_checks.py and outlier_drivers.py."""
    ds = pd.read_csv(os.path.join(OUT_DIR, "histone_checks.downsampling.tsv"), sep="\t")
    od = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.downsampling.tsv"), sep="\t"
    )
    groups = pd.read_csv(
        os.path.join(OUT_DIR, "outlier_drivers.groups.tsv"), sep="\t"
    ).set_index("set")
    fac = pd.read_csv(os.path.join(OUT_DIR, "histone_checks.dnm_factor.tsv"), sep="\t")
    histone_c, other_c = "#8f84e0", "#1a9e6e"
    fig, axes = plt.subplots(
        1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1.15, 1]}
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.2, wspace=0.42)

    ax = axes[0]
    for key, c in (
        ("histone genes, chr6p22 cluster", histone_c),
        ("all excess outliers", other_c),
    ):
        r = groups.loc[key]
        ax.axhspan(r.dnm_lo, r.dnm_hi, color=c, alpha=0.13, lw=0)
    ax.text(
        900,
        groups.loc["histone genes, chr6p22 cluster"].dnm_hi - 0.04,
        "DNM obs/exp (95% CI), histone cluster",
        fontsize=13,
        color=INK2,
        ha="left",
        va="top",
    )
    ax.text(
        1.9e5,
        groups.loc["all excess outliers"].dnm_hi + 0.03,
        "DNM obs/exp, the 46 non-histone outliers",
        fontsize=13,
        color=INK2,
        ha="right",
        va="bottom",
    )
    for anc, g in ds.groupby("ancestry"):
        g = g[g.n_samples >= 1000]
        if anc == "global":
            ax.plot(
                g.n_samples,
                g.rel_oe,
                color=histone_c,
                lw=2.8,
                marker="o",
                ms=7,
                zorder=3,
                label="histone cluster, all samples",
            )
        else:
            ax.plot(
                g.n_samples,
                g.rel_oe,
                color=histone_c,
                lw=1,
                alpha=0.6,
                marker="o",
                ms=3.5,
                label=(
                    "histone cluster, one ancestry group each" if anc == "afr" else None
                ),
            )
    g = od[(od.ancestry == "global") & (od.n_samples >= 1000)]
    ax.plot(
        g.n_samples,
        g.rel_oe,
        color=other_c,
        lw=2.2,
        marker="o",
        ms=5.5,
        label="46 non-histone outliers",
    )
    ax.axhline(1, color=AXIS, lw=1)
    ax.set_xscale("log")
    ax.set_xticks([1000, 10000, 100000, 500000], ["1k", "10k", "100k", "500k"])
    ax.set_xlim(800, 7e5)
    ax.set_ylim(0.6, 3.7)
    ax.legend(
        loc="upper left", bbox_to_anchor=(0.0, -0.2), ncol=2, frameon=False, fontsize=13
    )
    ax.set_title(
        "a  The histone excess is about 3× while few sites are seen,\n    and shrinks as sites saturate"
    )
    ax.set_xlabel("downsampled exomes (log scale)")
    ax.set_ylabel("synonymous obs/exp (no rr),\nrelative to genes without a flag")
    style_axis(ax, grid_axis="y")

    ax = axes[1]
    labels = [
        "no correction",
        "× 2.84\n(DNM excess)",
        "× 2.84,\nbefore\nsaturation",
        "× 5.2 non-CpG,\n× 1.15 CpG,\nbefore\nsaturation",
    ]
    x = np.arange(len(fac))
    for i, (col, name, c) in enumerate(
        (
            ("oe", "all sites", REC_C),
            ("oe_noncpg", "non-CpG", "#5f5e5a"),
            ("oe_cpg", "CpG", "#c3c2b7"),
        )
    ):
        ax.bar(x + (i - 1) * 0.26, fac[col], width=0.24, color=c, lw=0, label=name)
    ax.axhline(1, **THEORY)
    ax.text(len(fac) - 0.55, 1.03, "obs/exp = 1", fontsize=13, color=INK2, ha="right")
    ax.set_xticks(x, labels, fontsize=13)
    ax.set_ylim(0, 2.3)
    ax.legend(loc="upper right", frameon=False, fontsize=13)
    ax.set_title("b  A histone factor learned from DNMs\n    would overcorrect gnomAD")
    ax.set_ylabel("histone cluster synonymous obs/exp")
    style_axis(ax, grid_axis="y")

    fig.savefig(
        os.path.join(OUT_DIR, "rr_histones.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_concept():
    """Draw the conceptual diagram of the problem and the rr fix."""
    W, H = 1500, 720
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ramp = LinearSegmentedColormap.from_list("rate", ["#f6f3ea", "#e8c39e", "#c8692c"])

    # Left: the problem.
    ax.text(30, 42, "1  The problem", fontsize=17, fontweight="bold", color=INK)
    ax.text(
        30,
        76,
        "Expected counts come from each site's local sequence only.\nThe real rate also varies by region.",
        fontsize=12.5,
        color=INK2,
        va="top",
    )
    x0, x1, yb = 60, 690, 160
    grad = np.linspace(0, 1, 300)[None, :]
    ax.imshow(grad, extent=(x0, x1, yb + 38, yb), cmap=ramp, aspect="auto", zorder=1)
    ax.add_patch(
        Rectangle((x0, yb), x1 - x0, 38, fill=False, edgecolor=AXIS, lw=1, zorder=2)
    )
    ax.text(
        x0, yb + 50, "lower regional mutation rate", fontsize=11.5, color=INK2, va="top"
    )
    ax.text(
        x1,
        yb + 50,
        "higher regional mutation rate",
        fontsize=11.5,
        color=INK2,
        va="top",
        ha="right",
    )
    ax.text(
        (x0 + x1) / 2,
        yb - 12,
        "a chromosome (GC content, recombination, distance to the end, ...)",
        fontsize=11.5,
        color=INK2,
        ha="center",
    )
    genes = [("gene A", 175, 0.9), ("gene B", 575, 1.1)]
    base_y, unit = 540, 190
    for name, gx, oe in genes:
        ax.add_patch(
            Rectangle(
                (gx - 30, yb + 4),
                60,
                30,
                facecolor=SURFACE,
                edgecolor=INK,
                lw=1.5,
                zorder=3,
            )
        )
        ax.text(
            gx,
            yb + 19,
            name.split()[1],
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            zorder=4,
        )
        ax.text(
            gx,
            248,
            f"{name}: same sequence,\nso the same expected count",
            ha="center",
            va="top",
            fontsize=11.5,
            color=INK2,
        )
        for dx, h, fill, edge, lab in (
            (-40, 1.0, SURFACE, INK2, "expected"),
            (40, oe, NONE_C, NONE_C, "observed"),
        ):
            ax.add_patch(
                Rectangle(
                    (gx + dx - 30, base_y - unit * h * 0.85),
                    60,
                    unit * h * 0.85,
                    facecolor=fill,
                    edgecolor=edge,
                    lw=1.5,
                    ls=(0, (3, 2)) if lab == "expected" else "-",
                )
            )
            ax.text(gx + dx, base_y + 20, lab, ha="center", fontsize=11.5, color=INK2)
        ax.text(
            gx,
            base_y + 56,
            f"synonymous obs/exp = {oe:.1f}",
            ha="center",
            fontsize=13.5,
            fontweight="bold",
        )
    ax.text(
        175, base_y + 88, "looks more constrained", ha="center", fontsize=12, color=INK2
    )
    ax.text(
        575, base_y + 88, "looks less constrained", ha="center", fontsize=12, color=INK2
    )
    ax.text(
        375,
        base_y + 128,
        "Neither is under selection: synonymous obs/exp should be 1 everywhere.",
        ha="center",
        fontsize=12.5,
        color=INK,
    )
    ax.plot([735, 735], [30, H - 30], color=GRID, lw=1.2)

    # Right: the fix.
    ax.text(
        770,
        42,
        "2  The fix: a regional factor (rr)",
        fontsize=17,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        770,
        76,
        "learned from where de novo mutations land",
        fontsize=12.5,
        color=INK2,
        va="top",
    )
    bw, bh, fs = 335, 100, 11.5
    box(
        ax,
        775,
        120,
        bw,
        bh,
        SYN_FILL,
        SYN_EDGE,
        [
            "De novo mutations (DNMs)",
            "where mutation happened in one",
            "generation, in people not in gnomAD",
        ],
        size=fs,
    )
    box(
        ax,
        1140,
        120,
        bw,
        bh,
        SAME_FILL,
        SAME_EDGE,
        [
            "Background sites",
            "possible synonymous sites",
            "(now: no variant in gnomAD)",
        ],
        size=fs,
    )
    box(
        ax,
        880,
        265,
        490,
        bh,
        RR_FILL,
        RR_EDGE,
        [
            "One model per trinucleotide context",
            "which regional features predict where DNMs land",
            "13 features × 4 window sizes",
        ],
        size=fs,
        center=True,
    )
    box(
        ax,
        880,
        405,
        490,
        76,
        RR_FILL,
        RR_EDGE,
        ["rr for every 1-kb window and context", "about 0.85 to 1.25"],
        size=fs,
        center=True,
    )
    box(
        ax,
        880,
        520,
        490,
        76,
        REC_FILL,
        REC_EDGE,
        ["expected × rr", "for every variant class (synonymous, missense, LoF)"],
        size=fs,
        center=True,
    )
    arrow(ax, 940, 120 + bh, 1050, 265)
    arrow(ax, 1310, 120 + bh, 1200, 265)
    arrow(ax, 1125, 265 + bh, 1125, 405)
    arrow(ax, 1125, 481, 1125, 520)
    arrow(ax, 1125, 596, 1125, 626)
    ax.text(
        1125,
        652,
        "check on synonymous: genes A and B both at obs/exp ≈ 1",
        ha="center",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        1125,
        686,
        "open: which DNMs, which background, how strong, what's left over",
        ha="center",
        fontsize=11.5,
        color=INK2,
    )

    fig.savefig(
        os.path.join(OUT_DIR, "rr_concept.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_model():
    """Draw how the rr model is set up: training rows, per-context fit, and rr."""
    W, H = 1560, 660
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    # 1. Training rows.
    ax.text(
        30,
        40,
        "1  Training rows (illustrative, one context: ACA)",
        fontsize=15,
        fontweight="bold",
        color=INK,
    )
    cols = [
        ("site", 30),
        ("y", 175),
        ("GC, 1 kb", 225),
        ("male recomb., 1 Mb", 330),
        ("…", 480),
    ]
    for name, x in cols:
        ax.text(x + 8, 92, name, fontsize=11, fontweight="bold", color=INK2)
    rows = [
        ("chr1:1,204,331", 1, "61%", "1.9", SYN_FILL),
        ("chr7:5,530,612", 1, "57%", "2.4", SYN_FILL),
        ("chr3:48,101,920", 0, "44%", "0.6", SAME_FILL),
        ("chr12:6,441,075", 0, "49%", "1.1", SAME_FILL),
        ("chr5:140,212,804", 0, "39%", "0.4", SAME_FILL),
    ]
    for i, (site, y, gc, rec, fill) in enumerate(rows):
        top = 108 + i * 44
        ax.add_patch(
            FancyBboxPatch(
                (30, top),
                470,
                36,
                boxstyle="round,pad=0,rounding_size=4",
                facecolor=fill,
                edgecolor=SURFACE,
                lw=1,
            )
        )
        for val, (_, x) in zip((site, str(y), gc, rec, "…"), cols):
            ax.text(
                x + 8,
                top + 18,
                val,
                fontsize=11,
                va="center",
                color=INK,
                fontweight="bold" if val in ("0", "1") else "normal",
            )
    ax.add_patch(
        Rectangle((32, 352), 18, 14, facecolor=SYN_FILL, edgecolor=SYN_EDGE, lw=1)
    )
    ax.text(
        58,
        359,
        "y = 1: a synonymous site where a de novo mutation was seen",
        fontsize=11,
        va="center",
        color=INK,
    )
    ax.add_patch(
        Rectangle((32, 380), 18, 14, facecolor=SAME_FILL, edgecolor=SAME_EDGE, lw=1)
    )
    ax.text(
        58,
        387,
        "y = 0: a sampled possible synonymous site with no DNM",
        fontsize=11,
        va="center",
        color=INK,
    )
    ax.text(
        58,
        407,
        "and no variant in gnomAD exomes (5% sample, one row per locus)",
        fontsize=11,
        va="center",
        color=INK,
    )
    ax.text(
        32,
        445,
        "Each row carries the features of its 1-kb window:\n13 features × 4 window sizes = 52 columns.\nCpG contexts drop 5 methylation-linked features\n(32 columns), in every version.",
        fontsize=11,
        color=INK2,
        va="top",
    )
    arrow(ax, 515, 250, 565, 250, color=MUTED, lw=1.6)

    # 2. One logistic regression per context.
    ax.text(
        580,
        40,
        "2  One logistic regression per context",
        fontsize=15,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        580,
        70,
        "32 contexts: a base and its two neighbours, strand-collapsed",
        fontsize=11,
        color=INK2,
        va="top",
    )
    box(
        ax,
        580,
        110,
        420,
        92,
        RR_FILL,
        RR_EDGE,
        [
            "log odds that y = 1",
            "= β₀ + β₁z₁ + β₂z₂ + … + β₅₂z₅₂",
            "z = the window's features, z-scored in the context",
        ],
        size=12,
        center=True,
    )
    box(
        ax,
        580,
        236,
        200,
        150,
        NEW_FILL,
        NEW_EDGE,
        [
            "Current method",
            "screen each feature",
            "by p-value (Bonferroni)",
            "keep the passing ones",
            "fit with no penalty",
        ],
        size=11,
    )
    box(
        ax,
        800,
        236,
        200,
        150,
        REC_FILL,
        REC_EDGE,
        [
            "Ridge",
            "no screen: all 52 columns",
            "(32 at CpG contexts)",
            "penalty shrinks every β",
            "toward 0 (none to exactly 0)",
        ],
        size=11,
    )
    box(
        ax,
        580,
        402,
        420,
        76,
        RATE_FILL,
        RATE_EDGE,
        [
            "Per-variant rates + ridge (a quick first try)",
            "ridge, with log(per-variant rate) added to the log",
            "odds as a fixed offset: rr learns what the rates miss",
        ],
        size=11,
    )
    ax.text(
        790,
        496,
        "β₀ is the log odds at average features. How many y = 0\nrows are sampled only moves β₀, which cancels in rr.",
        fontsize=11,
        color=INK2,
        ha="center",
        va="top",
    )
    arrow(ax, 1015, 250, 1065, 250, color=MUTED, lw=1.6)

    # 3. From the model to rr (curve inset).
    ax.text(
        1080, 40, "3  From the model to rr", fontsize=15, fontweight="bold", color=INK
    )
    inset = fig.add_axes([1095 / W, 1 - 400 / H, 430 / W, 300 / H])
    eta = np.linspace(-4.2, -1.4, 200)
    p = 1 / (1 + np.exp(-eta))
    inset.plot(eta, p, color=REC_C, lw=2.2)
    b0, d = -3.0, 0.6
    for x, lab, c in (
        (b0, "average window", NONE_C),
        (b0 + d, "more mutable window", REC_C),
    ):
        y = 1 / (1 + np.exp(-x))
        inset.plot([x, x], [0, y], color=c, lw=1, ls=(0, (2, 2)))
        inset.scatter([x], [y], s=60, color=c, zorder=3, edgecolor=SURFACE, lw=1)
        inset.annotate(
            lab,
            (x, y),
            xytext=(-10, 10) if x == b0 else (12, -14),
            textcoords="offset points",
            fontsize=10.5,
            ha="right" if x == b0 else "left",
            color=INK,
        )
    rr = (1 / (1 + np.exp(-(b0 + d)))) / (1 / (1 + np.exp(-b0)))
    inset.text(
        0.03,
        0.95,
        f"rr = p(window) / p(average window)\n     = {rr:.2f} here (exaggerated)",
        transform=inset.transAxes,
        fontsize=10.5,
        va="top",
        color=INK,
    )
    inset.set_xticks([b0, b0 + d], ["β₀", "β₀ + β·z"])
    inset.set_yticks([])
    inset.set_ylabel("probability that y = 1", fontsize=10.5, color=INK2)
    inset.set_ylim(0, p.max() * 1.08)
    for side in ("top", "right"):
        inset.spines[side].set_visible(False)
    ax.text(
        1090,
        440,
        "rr = σ(β₀ + β·z) / σ(β₀) ≈ exp(β·z)",
        fontsize=12.5,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax.text(
        1090,
        470,
        "how much more (or less) mutable this window is than an\naverage window, for this context; "
        "real values are about 0.85–1.25",
        fontsize=11,
        color=INK2,
        va="top",
    )
    ax.text(
        1090,
        526,
        "Then: normalize so rr averages 1 over exome synonymous\nsites in each context, and multiply each variant's expected\n"
        "count by its window's rr (and, in the rate versions,\nby the variant's rate).",
        fontsize=11,
        color=INK2,
        va="top",
    )

    fig.savefig(
        os.path.join(OUT_DIR, "rr_model.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_problem_summary():
    """Draw the one-row summary of the problem for the top of the report."""
    tx = pd.read_csv(
        os.path.join(OUT_DIR, "data/rr_eval_export/transcripts.tsv.bgz"),
        sep="\t",
        compression="gzip",
    )
    zs = pd.read_csv(os.path.join(OUT_DIR, "syn_z_by_histone_group.tsv"), sep="\t")
    fig, axes = plt.subplots(
        1, 3, figsize=(18, 5.6), gridspec_kw={"width_ratios": [1, 1, 1.05]}
    )
    fig.subplots_adjust(left=0.05, right=0.985, top=0.85, bottom=0.17, wspace=0.5)

    # (a) Synonymous z of core genes vs chance.
    ax = axes[0]
    unflagged = tx.unflagged.astype(str).str.lower().eq("true")
    core = tx[unflagged & (tx.exp_unadj >= 20)]
    e = core.exp_unadj * core.obs.sum() / core.exp_unadj.sum()
    z = ((core.obs - e) / np.sqrt(e)).clip(-8, 12)
    edges = np.arange(-8, 12.5, 0.5)
    counts = np.histogram(z, edges)[0]
    ax.bar(edges[:-1] + 0.25, counts, width=0.46, color=NONE_C, lw=0)
    x = np.linspace(-8, 12, 400)
    chance = (
        len(z)
        * 0.5
        * np.exp(-(x**2) / (2 * CHANCE_FLOOR))
        / np.sqrt(2 * np.pi * CHANCE_FLOOR)
    )
    ax.plot(x, chance, **THEORY)
    top = max(chance.max(), counts.max()) * 1.08
    ax.text(
        2.6,
        top * 0.8,
        f"dotted: chance alone\n(mean z² ≤ {CHANCE_FLOOR:.1f})",
        fontsize=13,
        color=INK2,
    )
    ax.text(
        2.6,
        top * 0.55,
        f"bars: observed\n(mean z² = {np.mean(z ** 2):.1f})",
        fontsize=13,
        color=INK,
    )
    ax.set_ylim(0, top)
    ax.set_xlim(-8, 12)
    ax.set_title("a  Synonymous counts scatter\n    far more than chance")
    ax.set_xlabel("synonymous z = (obs − exp) / √exp, no correction")
    ax.set_ylabel(f"genes (of {len(z):,} without a flag)")
    style_axis(ax, grid_axis="y")

    # (b) obs/exp by the gene's DNM-predicted regional rate.
    ax = axes[1]
    d = zs.merge(tx[["transcript", "unflagged"]], on="transcript")
    d = d[d.unflagged.astype(str).str.lower().eq("true") & (d.exp_none >= 20)]
    d = d.assign(rr=d.exp_recommended / d.exp_none)
    d["bin"] = pd.qcut(d.rr.rank(method="first"), 10, labels=False)
    g = d.groupby("bin").agg(
        rr=("rr", "mean"),
        obs=("obs", "sum"),
        e0=("exp_none", "sum"),
        e1=("exp_recommended", "sum"),
    )
    k0, k1 = d.obs.sum() / d.exp_none.sum(), d.obs.sum() / d.exp_recommended.sum()
    lim = (0.94, 1.12)
    ax.plot(lim, lim, **THEORY, label="if rr were exactly right")
    ax.axhline(1, color=AXIS, lw=1)
    ax.plot(
        g.rr,
        g.obs / g.e0 / k0,
        color=NONE_C,
        lw=2.2,
        marker="o",
        ms=7,
        label="no correction",
    )
    ax.plot(
        g.rr,
        g.obs / g.e1 / k1,
        color=REC_C,
        lw=2.6,
        marker="D",
        ms=7,
        label="with the ridge rr",
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.legend(loc="upper left", frameon=False, fontsize=14)
    ax.set_title(
        "b  Genes in more mutable regions carry\n    more variants than the model expects"
    )
    ax.set_xlabel("gene's regional rate rr, learned from DNMs (decile mean)")
    ax.set_ylabel("synonymous obs/exp")
    style_axis(ax, grid_axis="y")

    # (c) How much of the scatter each correction removes.
    ax = axes[2]
    rec = pd.read_csv(os.path.join(OUT_DIR, "recommended_rr_metrics.tsv"), sep="\t")
    rows = [
        ("No correction", 2.908, NONE_C),
        ("Gnocchi adj_r", 5.450, ADJ_R_C),
        ("Legacy syn_adj_r", 2.509, LEGACY_C),
        (
            "Ridge rr",
            float(rec[rec.alpha == 1.0].core_noise.iloc[0]),
            REC_C,
        ),
    ]
    ypos = np.arange(len(rows))[::-1]
    for y, (label, v, c) in zip(ypos, rows):
        ax.barh(y, min(v, 3.35), height=0.58, color=c, lw=0)
        ax.text(
            min(v, 3.35) + 0.04,
            y,
            f"{v:.2f}" + (" (worse)" if v > 3.35 else ""),
            va="center",
            fontsize=14,
        )
    ax.axvline(CHANCE_FLOOR, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.04,
        len(rows) - 0.55,
        f"chance alone ≤ {CHANCE_FLOOR:.2f}",
        fontsize=13,
        color=INK2,
    )
    ax.annotate(
        "",
        xy=(CHANCE_FLOOR, -0.55),
        xytext=(rows[-1][1], -0.55),
        arrowprops=dict(arrowstyle="<->", color=INK2, lw=1),
    )
    ax.text(
        rows[-1][1] + 0.08,
        -0.55,
        "left after ridge: mostly\nreal gene-level rate\ndifferences",
        ha="left",
        va="center",
        fontsize=13,
        color=INK2,
    )
    ax.set_yticks(ypos, [r[0] for r in rows])
    ax.set_xlim(0, 4.6)
    ax.set_ylim(-1.2, len(rows) - 0.3)
    ax.set_title("c  What each correction removes")
    ax.set_xlabel("per-gene noise, mean z² (lower = better)")
    style_axis(ax)

    fig.savefig(
        os.path.join(OUT_DIR, "rr_problem_summary.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_background():
    """Draw the rr training background comparison, from rr_background_experiments.py."""
    met = pd.read_csv(os.path.join(OUT_DIR, "rr_background.metrics.tsv"), sep="\t")
    cls = pd.read_csv(os.path.join(OUT_DIR, "rr_background.class_slopes.tsv"), sep="\t")
    seen = pd.read_csv(
        os.path.join(OUT_DIR, "rr_background.frac_seen_by_rr.tsv"), sep="\t"
    )
    pairs = pd.read_csv(
        os.path.join(OUT_DIR, "rr_background.rr_absent_vs_all.tsv"), sep="\t"
    )
    class_c = CLASS_C
    backgrounds = [
        ("absent", "Absent from gnomAD (current rule)"),
        (
            "absent_context_matched",
            "Absent, sampled by context\n(in proportion to all sites)",
        ),
        ("absent_equal_per_context", "Absent, sampled evenly by context"),
        (
            "absent_reweighted",
            "Absent, reweighted by the model's\nexpected chance of being seen",
        ),
        ("all_sites", "All possible synonymous variants"),
    ]

    def get(bg, app="multiplicative"):
        return met[(met.background == bg) & (met.application == app)].iloc[0]

    fig = plt.figure(figsize=(15, 12.5))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.18)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.95, bottom=0.08)

    # (a) Share of the background the current rule drops, by regional rate.
    ax = fig.add_subplot(gs[0, 0])
    for c, g in seen.groupby("class"):
        g = g.sort_values("rr_all_sites")
        x = np.linspace(g.rr_all_sites.min(), g.rr_all_sites.max(), 50)
        q = g.class_frac_seen.iloc[0]
        ax.plot(x, 100 * (1 - (1 - q) ** x), **THEORY)
        ax.plot(
            g.rr_all_sites, 100 * g.frac_seen, color=class_c[c], lw=2, marker="o", ms=7
        )
        ax.text(
            g.rr_all_sites.iloc[-1] + 0.006,
            100 * g.frac_seen.iloc[-1],
            c,
            color=INK,
            fontsize=14,
            va="center",
        )
    ax.text(
        0.02,
        0.97,
        "dotted: if sites were seen at a rate\nproportional to rr alone",
        transform=ax.transAxes,
        va="top",
        fontsize=13,
        color=INK2,
    )
    ax.set_xlim(right=seen.rr_all_sites.max() + 0.05)
    ax.set_title(
        "a  How much background the current rule drops\n    changes with the regional rate"
    )
    ax.set_xlabel("regional rate rr (all-sites model, decile mean)")
    ax.set_ylabel(
        "% of possible synonymous variants seen in gnomAD\n(= dropped from the current background)"
    )
    style_axis(ax, grid_axis="y")

    # (b) rr learned against the current background vs against all sites.
    ax = fig.add_subplot(gs[0, 1])
    lo, hi = 0.86, 1.22
    x = np.array([lo, hi])
    ax.plot(x, x, **THEORY)
    ax.text(1.2, 1.12, "dotted: same rr", ha="right", va="top", fontsize=13, color=INK2)
    for c, g in pairs.groupby("class"):
        b = np.polyfit(g.log_rr_all_sites, g.log_rr_absent, 1)[0]
        ax.plot(
            np.exp(g.log_rr_all_sites),
            np.exp(g.log_rr_absent),
            color=class_c[c],
            lw=2,
            marker="o",
            ms=5.5,
            label=f"{c}: {b:.2f}× as steep",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ticks = [0.9, 1.0, 1.1, 1.2]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.legend(loc="upper left", frameon=False, fontsize=14)
    ax.set_title(
        "b  So the rr it learns is a little steeper\n    than against all sites"
    )
    ax.set_xlabel("rr against all possible sites (20 bins per class)")
    ax.set_ylabel("rr against the current background")
    style_axis(ax, grid_axis="y")

    # (c) Per-gene noise, held-out chromosomes.
    ax = fig.add_subplot(gs[1, 0])
    none = met[met.option == "No correction"].iloc[0]
    ypos = np.arange(len(backgrounds))[::-1]
    for y, (bg, _) in zip(ypos, backgrounds):
        r = get(bg)
        ax.scatter(
            r.oof_noise, y, s=112, color=REC_C, edgecolor=SURFACE, lw=1.2, zorder=3
        )
        txt = f"{r.oof_noise:.2f}"
        pc = met[(met.background == bg) & (met.application == "per-class strength")]
        if len(pc):
            v = pc.oof_noise.iloc[0]
            ax.plot([v, r.oof_noise], [y, y], color=REC_C, lw=1.2)
            ax.scatter(
                v,
                y,
                s=112,
                marker="D",
                color=SURFACE,
                edgecolor=REC_C,
                lw=1.6,
                zorder=3,
            )
            txt += f" → {v:.2f}"
        ax.text(
            max(r.oof_noise, none.core_noise * 0 + r.oof_noise) + 0.012,
            y + 0.18,
            txt,
            fontsize=13,
            va="bottom",
        )
    ax.axvline(none.core_noise, **THEORY)
    ax.text(
        none.core_noise - 0.008,
        len(backgrounds) - 0.45,
        "no correction",
        ha="right",
        fontsize=13,
        color=INK2,
    )
    ax.set_yticks(ypos, [b[1] for b in backgrounds])
    ax.set_xlim(2.45, 2.95)
    ax.set_ylim(-0.6, len(backgrounds) - 0.2)
    ax.set_title("c  Per-gene noise on held-out chromosomes")
    ax.set_xlabel(
        "mean z², overall level removed (lower = better)\n"
        "dot = rr as fitted; open diamond = with one strength per context class"
    )
    style_axis(ax)

    # (d) Calibration slope by context class.
    ax = fig.add_subplot(gs[1, 1])
    for y, (bg, _) in zip(ypos, backgrounds):
        for i, c in enumerate(class_c):
            s = cls[
                (cls.background == bg)
                & (cls.application == "multiplicative")
                & (cls["class"] == c)
            ].slope.iloc[0]
            ax.scatter(
                np.clip(s, -0.5, 3.0),
                y + (i - 1) * 0.17,
                s=72,
                color=class_c[c],
                edgecolor=SURFACE,
                lw=0.8,
                zorder=3,
                label=c if y == ypos[0] else None,
            )
            if s > 3.0:
                ax.text(
                    3.0,
                    y + (i - 1) * 0.17,
                    f"  {s:.0f}",
                    va="center",
                    fontsize=12.5,
                    color=INK2,
                )
    ax.axvline(1, **THEORY)
    ax.text(
        1.03, len(backgrounds) - 0.45, "1 = right strength", fontsize=13, color=INK2
    )
    ax.set_yticks(ypos, [])
    ax.set_xlim(0, 1.8)
    ax.set_ylim(-0.6, len(backgrounds) - 0.2)
    ax.set_title(
        "d  Calibration slope by context class\n    (below 1 = too strong, above 1 = too weak)"
    )
    ax.set_xlabel("calibration slope on gnomAD synonymous counts")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=14,
    )
    style_axis(ax)

    fig.savefig(
        os.path.join(OUT_DIR, "rr_background.local.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


LOGIT_EXAMPLES = ["CAG", "CCT", "CCG"]
"""The context with the most DNMs in each class."""

LOGIT_BACKGROUNDS = [
    ("absent", "Absent from gnomAD, per locus (current)", REC_C, "-"),
    ("absent_variant", "Absent from gnomAD, per variant", REC_C, (0, (4, 2.5))),
    (
        "absent_context_matched",
        "Absent, sampled by context in proportion",
        "#2a78d6",
        (0, (5, 2)),
    ),
    (
        "absent_equal_per_context",
        "Absent, sampled evenly by context",
        "#1baf7a",
        (0, (1, 1.5)),
    ),
    ("all_sites", "All sites, per locus", "#eb6834", "-"),
    ("all_variant", "All sites, per variant", "#eb6834", (0, (4, 2.5))),
]
"""The pipeline background runs (cluster/export_background_fits.py): key, label,
colour, line style."""


def read_logistic_fits(prefix="rr_logistic_fits"):
    """
    Read the per-context logistic fit tables.

    "rr_logistic_fits" is the local version (rr_logistic_fits.py) and
    "rr_logistic_fits.pipeline" the full pipeline run (cluster/export_logistic_fits.py).
    """
    return [
        pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.{t}.tsv"), sep="\t")
        for t in ("bins", "points", "contexts")
    ]


def logistic_panel(ax, context, bins, points, ctx, rng, small=False):
    """
    Draw one context's logistic regression: sites at y = 0 / 1 and the fitted curve.

    The curve is drawn with equal total weight for DNMs and background, which only
    moves the intercept, so that it is centered.
    """
    row = ctx.set_index("context").loc[context]
    c = CLASS_C[row["class"]]
    lo, hi = row.score_p05, row.score_p95
    xmax = 0.5
    ax.axvspan(lo, hi, color=GRID, alpha=0.7, lw=0, zorder=0)
    s = np.linspace(-xmax, xmax, 400)
    p = 1 / (1 + np.exp(-(row.intercept_equal_weight + s)))
    inside = (s >= lo) & (s <= hi)
    ax.plot(s, p, color=INK2, lw=1.1 if small else 1.4, ls=(0, (2, 2)))
    ax.plot(s[inside], p[inside], color=INK, lw=2 if small else 2.8)
    if not small:
        # The whole S, to show that sites only use its nearly straight middle.
        ins = ax.inset_axes([0.04, 0.5, 0.27, 0.24])
        sw = np.linspace(-6, 6, 300)
        ins.plot(sw, 1 / (1 + np.exp(-(row.intercept_equal_weight + sw))), color=INK2)
        ins.add_patch(
            Rectangle(
                (-xmax, 0.3), 2 * xmax, 0.4, facecolor=GRID, edgecolor=INK2, lw=0.8
            )
        )
        ins.set_xlim(-6, 6)
        ins.set_ylim(-0.05, 1.05)
        ins.set_xticks([-5, 0, 5])
        ins.set_yticks([0, 1])
        ins.tick_params(labelsize=11, length=0, labelcolor=MUTED)
        ins.set_facecolor(SURFACE)
        for side in ins.spines.values():
            side.set_color(AXIS)
        ins.set_title(
            "whole curve; box = this panel", fontsize=12, color=INK2, weight="normal"
        )
    pts = points[points.context == context]
    for y, (y0, y1) in ((1, (1.03, 1.13)), (0, (-0.13, -0.03))):
        v = pts[pts.y == y].score.clip(-xmax, xmax)
        ax.scatter(
            v,
            rng.uniform(y0, y1, len(v)),
            s=4 if small else 9,
            color=c,
            alpha=0.35,
            lw=0,
        )
    b = bins[(bins.context == context) & (bins.background == "absent")]
    ax.scatter(
        b.score,
        b.observed_p_equal_weight,
        s=19 if small else 51,
        color=SURFACE,
        edgecolor=INK,
        lw=0.9 if small else 1.2,
        zorder=3,
    )
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.2 if small else -0.3, 1.2)
    ax.set_yticks([0, 0.5, 1])
    style_axis(ax, grid_axis="y")
    if small:
        ax.set_title(
            f"{context}  ({row.dnms:.0f})", fontsize=14, color=INK, weight="bold"
        )
        ax.tick_params(labelsize=12)
        return
    ax.text(
        -xmax + 0.05,
        1.16,
        "DNM sites (y = 1)",
        fontsize=13,
        color=INK2,
        va="bottom",
    )
    ax.text(
        -xmax + 0.05,
        -0.16,
        "background sites (y = 0)",
        fontsize=13,
        color=INK2,
        va="top",
    )
    ax.text(
        min(hi + 0.04, xmax - 0.36),
        0.2,
        f"shaded: middle 90%\nof sites, rr {np.exp(lo):.2f}–{np.exp(hi):.2f}",
        fontsize=13,
        color=INK2,
    )


def rr_curve_panel(ax, context, bins, ctx, small=False):
    """Draw one context's DNM rate per score decile against the model's rr."""
    r = ctx.set_index("context").loc[context]
    b = bins[(bins.context == context) & (bins.background == "absent")]
    n0 = b.background_sites
    k = b.dnms.sum() / n0.sum()
    rel = b.dnms / n0 / k
    se = rel * np.sqrt(1 / b.dnms.clip(lower=1) + 1 / n0)
    s = np.linspace(b.score.min(), b.score.max(), 100)
    ax.plot(s, np.exp(s), color=INK, lw=1.6 if small else 2.4, label="model rr (≈ eˢ)")
    ax.errorbar(
        b.score,
        rel,
        yerr=1.96 * se,
        fmt="o",
        color=CLASS_C[r["class"]],
        mfc=SURFACE,
        mew=1 if small else 1.4,
        ms=4.5 if small else 9.5,
        elinewidth=0.8 if small else 1.1,
        capsize=0,
        label="DNM rate, relative",
    )
    ax.axhline(1, color=AXIS, lw=1, zorder=0)
    ax.set_yscale("log")
    ticks = [0.8, 1.0, 1.25, 1.6] if small else [0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6]
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(0.6, 1.75)
    style_axis(ax, grid_axis="y")
    if small:
        ax.tick_params(labelsize=12)


def score_dist_panel(ax, context, hist, ctx, small=False):
    """Draw one context's score distributions for DNM and background sites."""
    r = ctx.set_index("context").loc[context]
    for y, color, label in (
        (0, MUTED, "background sites"),
        (1, CLASS_C[r["class"]], "DNM sites"),
    ):
        h = hist[(hist.context == context) & (hist.y == y)]
        ax.plot(
            (h.lo + h.hi) / 2,
            h.share / (h.hi - h.lo),
            color=color,
            lw=(1.6 if y else 1.2) if small else (2.2 if y else 1.8),
            label=label,
        )
        ax.axvline(
            r.mean_dnm if y else r.mean_background,
            color=color,
            lw=1 if small else 1.4,
            ls=(0, (4, 2)),
        )
    shift = (r.mean_dnm - r.mean_background) / r.sd_background
    if small:
        ax.text(
            0.03,
            0.97,
            f"Δ {shift:.2f} SD",
            transform=ax.transAxes,
            va="top",
            fontsize=12,
            color=INK2,
        )
        ax.set_ylim(bottom=0)
        ax.set_yticks([])
        style_axis(ax, grid_axis="y")
        ax.tick_params(labelsize=12)
        return
    ax.text(
        0.0,
        1.02,
        f"means (dashed) differ by {shift:.2f} SD of the background",
        transform=ax.transAxes,
        va="bottom",
        fontsize=13,
        color=INK2,
    )
    ax.set_ylim(bottom=0)
    style_axis(ax, grid_axis="y")


def draw_logistic_fits(prefix="rr_logistic_fits", sites="possible synonymous sites"):
    """
    Draw the per-context logistic regressions.

    :param prefix: Table prefix (see ``read_logistic_fits``), also the figure names.
    :param sites: Which sites the shaded middle 90% and panel d are over.
    """
    bins, points, ctx = read_logistic_fits(prefix)
    hist_path = os.path.join(OUT_DIR, f"{prefix}.hist.tsv")
    # With score histograms (cluster/export_logistic_fits.py), two more rows show
    # the DNM rate per score decile against rr, and the two score distributions.
    stacked = os.path.exists(hist_path)
    rng = np.random.default_rng(0)
    if stacked:
        hist = pd.read_csv(hist_path, sep="\t")
        fig = plt.figure(figsize=(15, 18.5))
        gs = fig.add_gridspec(
            4, 3, height_ratios=[1, 0.85, 0.85, 0.9], hspace=0.6, wspace=0.16
        )
        fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.07)
    else:
        fig = plt.figure(figsize=(15, 9.6))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.9], hspace=0.62, wspace=0.16)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.12)
    fig.suptitle(
        "One logistic regression per context: how far the regional score moves a "
        "site's chance of being a DNM",
        x=0.07,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    for i, context in enumerate(LOGIT_EXAMPLES):
        ax = fig.add_subplot(gs[0, i])
        logistic_panel(ax, context, bins, points, ctx, rng)
        row = ctx.set_index("context").loc[context]
        ax.set_title(
            f"{'abc'[i]}  {context} ({row['class']}): {row.dnms:.0f} DNMs", pad=14
        )
        ax.set_xlabel("regional score s = β·z (log-odds; rr ≈ eˢ)")
        if i == 0:
            ax.set_ylabel("chance a site is a DNM\n(equal DNMs and background)")
    note_y = ax.get_position().y0 - 0.045
    fig.text(
        0.07,
        note_y,
        "solid: the fitted curve where sites fall · dashed: the same curve beyond them"
        " · open circles: share of DNMs in each score decile, on held-out chromosomes",
        fontsize=13,
        color=INK2,
    )
    if stacked:
        for i, context in enumerate(LOGIT_EXAMPLES):
            ax = fig.add_subplot(gs[1, i])
            rr_curve_panel(ax, context, bins, ctx)
            ax.set_title(f"{'def'[i]}  {context}: DNMs against rr, per score decile")
            ax.set_xlabel("regional score s (decile mean, held-out)")
            if i == 0:
                ax.set_ylabel("relative DNM rate")
            if i == 2:
                ax.legend(loc="upper left", frameon=False, fontsize=12)
            ax = fig.add_subplot(gs[2, i])
            score_dist_panel(ax, context, hist, ctx)
            ax.set_title(f"{'ghi'[i]}  {context}: score distributions", pad=22)
            ax.set_xlabel("regional score s (held-out)")
            if i == 0:
                ax.set_ylabel("density")
                ax.legend(loc="upper right", frameon=False, fontsize=13)

    ax = fig.add_subplot(gs[3 if stacked else 1, :])
    ctx = ctx.assign(
        order=ctx["class"].map({c: i for i, c in enumerate(CLASS_C)})
    ).sort_values(["order", "context"])
    x = np.arange(len(ctx)) + 0.8 * ctx.order.to_numpy()
    for xi, (_, r) in zip(x, ctx.iterrows()):
        c = CLASS_C[r["class"]]
        ax.plot(
            [xi, xi],
            [np.exp(r.score_p05), np.exp(r.score_p95)],
            color=c,
            lw=4,
            solid_capstyle="round",
        )
    for c in CLASS_C:
        xs = x[ctx["class"].to_numpy() == c]
        ax.text(xs.min() - 0.3, 1.49, c, fontsize=15, color=INK, weight="bold")
    ax.axhline(1, color=AXIS, lw=1, zorder=0)
    ax.set_xticks(
        x, [f"{r.context} ({r.dnms:.0f})" for _, r in ctx.iterrows()], rotation=90
    )
    ax.tick_params(axis="x", labelsize=13, labelcolor=INK2)
    ax.set_ylim(0.75, 1.55)
    ax.set_xlim(x.min() - 0.8, x.max() + 0.8)
    ax.set_ylabel(f"rr of the middle 90%\nof {sites.split()[0]} sites")
    ax.set_title(
        f"{'j' if stacked else 'd'}  Every context's model: the range of rr it gives"
        f" the middle 90% of {sites} (DNMs per context in brackets)"
    )
    style_axis(ax, grid_axis="y")
    fig.savefig(
        os.path.join(OUT_DIR, f"{prefix}.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)

    if stacked:
        draw_logistic_fits_all_stacked(prefix, bins, points, ctx, hist, rng)
        return
    fig, axes = plt.subplots(4, 8, figsize=(17, 9.5), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.05, right=0.99, top=0.9, bottom=0.08, hspace=0.45, wspace=0.12
    )
    fig.suptitle(
        "Every context's logistic regression (DNMs in brackets; colour = context"
        " class; drawn with equal numbers of DNMs and background)",
        x=0.05,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    for ax, (_, r) in zip(axes.ravel(), ctx.iterrows()):
        logistic_panel(ax, r.context, bins, points, ctx, rng, small=True)
    fig.supxlabel("regional score s = β·z (log-odds; rr ≈ eˢ)", fontsize=15.5)
    fig.supylabel("chance a site is a DNM", fontsize=15.5)
    fig.savefig(
        os.path.join(OUT_DIR, f"{prefix}.all_contexts.png"),
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_logistic_fits_all_stacked(prefix, bins, points, ctx, hist, rng):
    """
    Draw every context's logistic regression with the main figure's three stacked rows.

    For each context, top to bottom: the fitted curve (as panels a–c of the main
    figure), DNMs against rr per score decile (d–f), and the two score distributions
    (g–i). Contexts are in class order, eight per block.

    :param prefix: Table prefix, also the figure name.
    :param bins: Per-context score-decile table.
    :param points: Sampled sites.
    :param ctx: Per-context table, sorted in class order.
    :param hist: Score histograms.
    :param rng: Random generator for the jitter.
    """
    n_col, n_block = 8, int(np.ceil(len(ctx) / 8))
    height = 5.0 * n_block + 1.5
    fig = plt.figure(figsize=(18, height))
    gs = fig.add_gridspec(
        4 * n_block,
        n_col,
        height_ratios=([1, 0.8, 0.7, 0.3] * n_block)[:-1] + [0.01],
        hspace=0.2,
        wspace=0.14,
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=1 - 1.75 / height, bottom=0.03)
    fig.text(
        0.06,
        1 - 0.3 / height,
        "Every context's logistic regression, in the main figure's three rows"
        " (DNMs in brackets; colour = context class)",
        fontsize=17.5,
        weight="bold",
        va="top",
    )
    fig.text(
        0.06,
        1 - 0.68 / height,
        "top: the fitted curve, drawn with equal numbers of DNMs and background; open"
        " circles = share of DNMs per score decile\nmiddle: DNMs per background site"
        " (relative, ± 95% CI) against the model's rr (black)\nbottom: the score"
        " distributions of DNM (colour) and background (gray) sites; Δ = difference of"
        " the means, in SD of the background scores. All on held-out chromosomes.",
        fontsize=14,
        color=INK2,
        va="top",
        linespacing=1.5,
    )
    shared = None
    for k, (_, r) in enumerate(ctx.iterrows()):
        blk, col = divmod(k, n_col)
        axes = []
        for row in range(3):
            ax = fig.add_subplot(gs[4 * blk + row, col], sharex=shared)
            shared = shared or ax
            axes.append(ax)
        logistic_panel(axes[0], r.context, bins, points, ctx, rng, small=True)
        rr_curve_panel(axes[1], r.context, bins, ctx, small=True)
        score_dist_panel(axes[2], r.context, hist, ctx, small=True)
        for ax in axes[:2]:
            ax.tick_params(labelbottom=False)
        if col == 0:
            axes[0].set_ylabel("chance of DNM", fontsize=13)
            axes[1].set_ylabel("DNM rate,\nrelative", fontsize=13)
            axes[2].set_ylabel("density", fontsize=13)
        else:
            for ax in axes[:2]:
                ax.tick_params(labelleft=False)
        axes[2].set_xlabel("regional score s", fontsize=12.5, labelpad=1)
    shared.set_xlim(-0.35, 0.5)
    fig.savefig(
        os.path.join(OUT_DIR, f"{prefix}.all_contexts.png"),
        dpi=140,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_logistic_backgrounds(prefix="rr_logistic_fits.pipeline_bg6_dnm10"):
    """
    Draw how the logistic fits change with the training background.

    "rr_logistic_fits.pipeline_bg6_dnm10" is the pipeline runs with a penalty of 10
    per DNM, the same shrinkage for every background size ("..._bg6" is the default
    per-site penalty; cluster/export_background_fits.py).
    """
    bins, ctx = (
        pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.{t}.tsv"), sep="\t")
        for t in ("bins", "contexts")
    )
    cls = ctx.set_index("context")["class"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 11.5))
    fig.subplots_adjust(
        left=0.08, right=0.98, top=0.86, bottom=0.17, hspace=0.38, wspace=0.2
    )
    fig.suptitle(
        "The same logistic regressions under each training background",
        x=0.08,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    fig.text(
        0.08,
        0.94,
        "Ridge with a penalty of 10 per DNM, so backgrounds of different sizes get the"
        " same shrinkage. a–c: lines = fitted curves, dots = observed share of DNMs\nper"
        " score decile. d–f: lines = rr under each background; circles = DNMs per"
        " possible site in each decile, held out, ± 95% CI; dotted = rr equal to the"
        " DNM rate.",
        fontsize=13.5,
        color=INK2,
        va="top",
        linespacing=1.5,
    )
    for i, context in enumerate(LOGIT_EXAMPLES):
        ax = axes[0, i]
        g = bins[bins.context == context]
        for bg, label, c, ls in LOGIT_BACKGROUNDS:
            b = g[g.background == bg]
            ax.plot(b.score, b.fitted_p, color=c, ls=ls, lw=2, label=label)
            ax.scatter(
                b.score,
                b.observed_p,
                s=26,
                color=c,
                edgecolor=SURFACE,
                lw=0.6,
                zorder=3,
            )
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.yaxis.set_minor_formatter(plt.FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.tick_params(axis="y", which="both", labelsize=13, labelcolor=INK2)
        ax.set_title(f"{'abc'[i]}  {context} ({cls[context]}): the fitted curves")
        ax.set_xlabel("regional score under the current rule's model")
        if i == 0:
            ax.set_ylabel("share of training rows that are DNMs\n(log scale)")
        style_axis(ax, grid_axis="y")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower left",
        bbox_to_anchor=(0.07, 0.0),
        ncol=3,
        frameon=False,
        fontsize=14,
        handlelength=3,
    )

    for i, c in enumerate(CLASS_C):
        ax = axes[1, i]
        g = bins[bins.context.map(cls) == c]
        ref = g[g.background == "absent"]
        tot = ref.groupby("context")[["dnms", "possible"]].sum()
        ref = ref.assign(
            exp_dnms=ref.possible * ref.context.map(tot.dnms / tot.possible),
            w_rr=ref.possible * ref.rr,
        )
        by = ref.groupby("bin")[["dnms", "exp_dnms", "possible", "w_rr"]].sum()
        x = by.w_rr / by.possible
        rel = by.dnms / by.exp_dnms
        se = np.sqrt(by.dnms) / by.exp_dnms
        lo, hi = 0.78, 1.3
        ax.plot([lo, hi], [lo, hi], **THEORY)
        for bg, label, col, ls in LOGIT_BACKGROUNDS:
            b = g[g.background == bg].assign(w_rr=lambda d: d.possible * d.rr)
            b = b.groupby("bin")[["w_rr", "possible"]].sum()
            ax.plot(x, b.w_rr / b.possible, color=col, ls=ls, lw=2)
        ax.errorbar(
            x,
            rel,
            yerr=1.96 * se,
            fmt="o",
            ms=8,
            color=INK,
            mfc=SURFACE,
            mew=1.3,
            elinewidth=1,
            capsize=0,
            zorder=3,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ticks = [0.8, 0.9, 1.0, 1.1, 1.2]
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
        ax.minorticks_off()
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(f"{'def'[i]}  {c} ({len(tot)} contexts): rr and DNMs")
        ax.set_xlabel("rr under the current rule (decile mean, held-out)")
        if i == 0:
            ax.set_ylabel(
                "rr under each background (lines)\nDNMs per possible site, relative"
                " (circles)"
            )
        style_axis(ax, grid_axis="y")
    fig.savefig(
        os.path.join(OUT_DIR, "rr_logistic_fits.backgrounds.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


def draw_histone_conversion():
    """Draw the histone paralog-copying test, from histone_gene_conversion.py."""
    ex = pd.read_csv(
        os.path.join(OUT_DIR, "histone_gene_conversion.excess.tsv"), sep="\t"
    )
    tests = pd.read_csv(
        os.path.join(OUT_DIR, "histone_gene_conversion.tests.tsv"), sep="\t"
    )
    classes = [
        ("copies no paralog", "alt copies no paralog"),
        ("copies a paralog", "alt copies a paralog"),
        ("copies the closest paralog", "alt copies the closest paralog"),
        ("all", "all synonymous alleles"),
    ]
    fig, axes = plt.subplots(
        1, 3, figsize=(16.5, 4.8), gridspec_kw={"width_ratios": [1, 1, 1.05]}
    )
    fig.subplots_adjust(left=0.14, right=0.98, top=0.8, bottom=0.17, wspace=0.5)
    fig.suptitle(
        "Do the histone genes' extra variants copy a paralog's base?",
        x=0.14,
        ha="left",
        fontsize=17.5,
        weight="bold",
    )
    y = np.arange(len(classes))[::-1]
    for ax, (source, title, xmax) in zip(
        axes[:2],
        [
            (
                "gnomAD v4.1 exomes",
                "a  gnomAD: observed / expected\n    (expected from nearby non-histone genes)",
                2.8,
            ),
            (
                "DNMs",
                "b  DNMs: observed / expected\n    (expected from Roulette's rate per variant)",
                13.5,
            ),
        ],
    ):
        for yi, (cls, _) in zip(y, classes):
            r = ex[(ex.source == source) & (ex.alleles == cls)].iloc[0]
            c = INK2 if cls == "all" else REC_C
            ax.plot([r.ratio_lo, r.ratio_hi], [yi, yi], color=c, lw=1.6)
            ax.scatter(r.ratio, yi, s=96, color=c, edgecolor=SURFACE, lw=1.2, zorder=3)
            ax.text(
                r.ratio_hi + xmax * 0.02,
                yi + 0.16,
                f"{r.ratio:.2f} ({r.observed:.0f} / {r.expected:.0f})",
                fontsize=13,
                color=INK2,
                va="bottom",
            )
        ax.axvline(1, **THEORY)
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.6, len(classes) - 0.3)
        ax.set_title(title)
        ax.set_xlabel("observed / expected (95% interval)")
        style_axis(ax)
    axes[0].set_yticks(y, [c[1] for c in classes])
    axes[1].set_yticks(y, [])

    ax = axes[2]
    fams = ["all histone genes", "H1", "H2A", "H2B", "H3", "H4"]
    yf = np.arange(len(fams))[::-1]
    for yi, fam in zip(yf, fams):
        for dy, flag, filled in (
            (0.13, "alt_in_paralog", True),
            (-0.13, "alt_in_paralog, with Roulette site rate", False),
        ):
            r = tests[(tests.set == fam) & (tests.flag == flag)].iloc[0]
            ax.plot(
                [r.or_flagged_lo, r.or_flagged_hi], [yi + dy] * 2, color=REC_C, lw=1.4
            )
            ax.scatter(
                r.or_flagged,
                yi + dy,
                s=74,
                color=REC_C if filled else SURFACE,
                edgecolor=REC_C,
                lw=1.4,
                zorder=3,
            )
    ax.axvline(1, **THEORY)
    ax.set_yticks(
        yf, [f"{f} ({tests[tests.set == f].transcripts.iloc[0]})" for f in fams]
    )
    ax.set_xlim(0.7, 2.1)
    ax.set_ylim(-0.6, len(fams) - 0.4)
    ax.set_title("c  Within the genes: alleles that copy\n    a paralog against others")
    ax.set_xlabel(
        "odds ratio (95% interval)\nfilled: by trinucleotide · open: + Roulette site rate"
    )
    style_axis(ax)
    fig.savefig(
        os.path.join(OUT_DIR, "rr_histone_conversion.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


if __name__ == "__main__":
    draw_diagram()
    draw_scorecard()
    draw_new_metrics()
    draw_rate_parts()
    draw_full_run_release()
    draw_performance()
    draw_deciles()
    draw_deciles(
        "syn_oe_by_decile.full_run.explained_removed.tsv",
        "rr_by_decile.explained_removed.png",
        "Explained synonymous outliers removed: 91 histone, repeat or segdup, and maternal-hotspot genes",
    )
    draw_outliers()
    draw_alternatives()
    draw_leftover_noise()
    draw_outlier_drivers()
    draw_background()
    draw_background_pipeline()
    draw_background_bias()
    draw_concept()
    draw_model()
    draw_dnm_sets()
    draw_histones()
    draw_problem_summary()
    draw_features()
    draw_feature_use()
    draw_penalty()
    draw_rates_diagram()
    draw_logistic_fits()
    draw_logistic_fits("rr_logistic_fits.pipeline", "background sites")
    draw_logistic_backgrounds()
    draw_histone_conversion()
