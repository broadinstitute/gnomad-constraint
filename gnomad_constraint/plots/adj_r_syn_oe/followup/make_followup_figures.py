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

SAME_FILL, SAME_EDGE = "#f0efec", "#c3c2b7"
SYN_FILL, SYN_EDGE = "#cde2fb", "#2a78d6"
NEW_FILL, NEW_EDGE = "#fbe1d6", "#eb6834"
RR_FILL, RR_EDGE = "#fdf0cc", "#eda100"
REC_FILL, REC_EDGE = "#ffffff", REC_C

CLASS_C = {"A/T": "#1a9e6e", "C non-CpG": "#8f84e0", "CpG": "#d9622b"}

THEORY = dict(color=INK2, lw=1.3, ls=(0, (1, 2)))
CHANCE_FLOOR = 0.70
"""Upper bound on mean z² from sampling alone for the core set: each site is
observed or not, so a transcript's variance is at most sum(n p (1 - p)) per
(window, context) cell, not its expected count."""

plt.rcParams.update(
    {
        "font.family": ["Arial", "DejaVu Sans"],
        "font.size": 10.5,
        "text.color": INK,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": INK2,
        "axes.titlesize": 11.5,
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
    """Draw the shared pipeline, the three existing rr builds and the recommended one."""
    W, H = 1300, 955
    fig = plt.figure(figsize=(W / 100, H / 100))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    # A. Shared pipeline.
    ax.text(
        10,
        18,
        "A. Constraint pipeline: identical in every version",
        fontsize=11.5,
        fontweight="bold",
        va="center",
    )
    ax.text(
        10,
        38,
        "No correction = steps 1 → 4 → 6. Every rr version adds step 5 and changes"
        " nothing else.",
        fontsize=9.5,
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
            "per variant, before",
            "summing (skipped in",
            "'no correction')",
        ],
        ["6  Sum per transcript", "obs/exp, LOEUF,", "z-scores, flags"],
    ]
    bw, gap, by, bh = 190, 24, 56, 82
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
            size=8.8,
            dash=is_rr,
        )
        if i < len(top) - 1:
            arrow(ax, x + bw + 2, by + bh / 2, x + bw + gap - 2, by + bh / 2)

    # B. rr recipes, three lanes.
    lane_x0, lane_w, lane_gap = 165, 270, 12
    lane_x = [lane_x0 + i * (lane_w + lane_gap) for i in range(4)]
    lane_cx = [x + lane_w / 2 for x in lane_x]
    rr_cx = 10 + 4 * (bw + gap) + bw / 2
    yb = 172
    ax.text(
        10, yb, "B. How the rr table", fontsize=11.5, fontweight="bold", va="center"
    )
    ax.text(
        10,
        yb + 18,
        "for step 5 is built",
        fontsize=11.5,
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
    ]
    hy, hh = 188, 44
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
        ax.scatter([x + 16], [hy + 15], s=60, color=c, zorder=3)
        ax.text(x + 28, hy + 15, name, fontsize=10, fontweight="bold", va="center")
        ax.text(x + 12, hy + 32, sub[:44], fontsize=8.5, color=INK2, va="center")

    S, B, N, R = "same", "syn", "new", "rec"
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
                (
                    N,
                    [
                        "Satterstrom 2026 DNMs",
                        "or the union with the 17k set,",
                        "pending the 282 dropped DNMs",
                    ],
                ),
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
            ],
        ),
        (
            "Code",
            [
                (S, ["atgu/gnomad_nc_constraint", "plus released model pickles"]),
                (B, ["Sandbox notebooks", "not in a repo"]),
                (N, ["regional_rate_pipeline.py", "in gnomad-constraint"]),
                (N, ["Ruchit's pipeline", "ridge in fit_rr_context_model"]),
            ],
        ),
    ]
    styles = {
        S: (SAME_FILL, SAME_EDGE),
        B: (SYN_FILL, SYN_EDGE),
        N: (NEW_FILL, NEW_EDGE),
        R: (REC_FILL, REC_EDGE),
    }
    ry, rh, rgap = hy + hh + 14, 62, 12
    for r, (label, cells) in enumerate(rows):
        y = ry + r * (rh + rgap)
        ax.text(
            150,
            y + rh / 2,
            label,
            ha="right",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color=INK2,
            linespacing=1.3,
        )
        for x, (kind, lines) in zip(lane_x, cells):
            fill, edge = styles[kind]
            lw = 2.2 if kind == R else 1.3
            if len(lines) == 1:
                box(ax, x, y, lane_w, rh, fill, edge, lines, size=9, center=True, lw=lw)
            else:
                box(ax, x, y, lane_w, rh, fill, edge, lines, size=9, lw=lw)
            if r < len(rows) - 1:
                arrow(ax, x + lane_w / 2, y + rh + 1, x + lane_w / 2, y + rh + rgap - 1)
        if r == 0:
            y_first = y
    y_end = ry + len(rows) * (rh + rgap)
    for x in lane_x:
        arrow(ax, x + lane_w / 2, hy + hh + 1, x + lane_w / 2, y_first - 1)

    # Legend.
    ly = y_end + 14
    items = [
        (SAME_FILL, SAME_EDGE, "Same as Gnocchi's original"),
        (SYN_FILL, SYN_EDGE, "Changed for syn_adj_r (Ruchit keeps it)"),
        (NEW_FILL, NEW_EDGE, "New in Ruchit's pipeline"),
        (REC_FILL, REC_EDGE, "Changed for the ridge fit"),
        (RR_FILL, RR_EDGE, "Where rr enters"),
    ]
    lx = 165
    for fill, edge, text in items:
        ax.add_patch(
            FancyBboxPatch(
                (lx, ly - 8),
                22,
                16,
                boxstyle="round,pad=0,rounding_size=3",
                facecolor=fill,
                edgecolor=edge,
                linewidth=1.3,
            )
        )
        ax.text(lx + 30, ly, text, fontsize=9.5, va="center", color=INK2)
        lx += 30 + len(text) * 6.6 + 26

    fig.savefig(os.path.join(OUT_DIR, "rr_pipelines_diagram.png"), dpi=200)
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


def draw_performance():
    """Draw the rr performance panels, mostly from Ruchit's v2 report."""
    rec = pd.read_csv(os.path.join(OUT_DIR, "recommended_rr_metrics.tsv"), sep="\t")
    rec1 = rec[rec.alpha == 1.0].iloc[0]
    fig, axes = plt.subplots(3, 2, figsize=(11, 12.4))
    fig.subplots_adjust(
        left=0.2, right=0.97, top=0.96, bottom=0.06, hspace=0.48, wspace=0.62
    )

    names = [
        "No correction",
        "Gnocchi adj_r",
        "Legacy syn_adj_r",
        "Ruchit, 17k DNMs",
        "Ruchit, Satterstrom\nDNMs",
        "Ridge, Satterstrom\nDNMs",
    ]
    colors = [NONE_C, ADJ_R_C, LEGACY_C, RUCHIT_C, RUCHIT_C, REC_C]
    hatches = [None, None, None, None, "///", None]
    ypos = list(range(len(names)))[::-1]

    # (a) Per-gene noise.
    ax = axes[0, 0]
    vals = [2.908, 5.450, 2.509, 2.557, 2.696, rec1.core_noise]
    for y, v, c, h in zip(ypos, vals, colors, hatches):
        ax.barh(y, v, height=0.56, color=c, hatch=h, hatchcolor="white", linewidth=0)
        ax.text(
            v - 0.1,
            y,
            f"{v:.2f}",
            ha="right",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color="white",
            bbox=dict(boxstyle="round,pad=0.15", facecolor=c, linewidth=0),
        )
    ax.axvline(2.908, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
    ax.axvline(CHANCE_FLOOR, zorder=0, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.08,
        len(names) - 0.45,
        f"chance ≤ {CHANCE_FLOOR:.2f}",
        fontsize=8.8,
        color=INK2,
        va="center",
    )
    ax.set_yticks(ypos, names)
    ax.set_xlim(0, 6.0)
    ax.set_ylim(-0.6, len(names) - 0.15)
    ax.set_title("a  Per-gene synonymous noise")
    ax.set_xlabel("mean z², overall level removed (lower = better)")
    style_axis(ax)

    # (b) Synonymous outliers, release-flagged vs unflagged genes.
    ax = axes[0, 1]
    out = read_outliers()
    corrs = [
        "none",
        "adj_r",
        "syn_adj_r",
        "ruchit_17k",
        "ruchit_satterstrom",
        "recommended",
    ]
    for y, corr, c, h in zip(ypos, corrs, colors, hatches):
        flagged = sum(
            out_count(out, corr, d, f)
            for d in ("deficit", "excess")
            for f in ("low_exome_mapping_quality", "low_exome_coverage", "other flag")
        )
        unflagged = (
            sum(out_count(out, corr, d) for d in ("deficit", "excess")) - flagged
        )
        ax.barh(
            y, flagged, height=0.56, color="#c3c2b7", edgecolor=SURFACE, linewidth=2
        )
        ax.text(
            flagged / 2,
            y,
            f"{flagged}",
            ha="center",
            va="center",
            fontsize=9.5,
            color=INK,
        )
        ax.barh(
            y,
            unflagged,
            left=flagged,
            height=0.56,
            color=c,
            hatch=h,
            hatchcolor="white",
            edgecolor=SURFACE,
            linewidth=2,
        )
        ax.text(
            flagged + unflagged + 8,
            y,
            f"{unflagged}",
            va="center",
            fontsize=9.5,
            color=INK,
        )
    ax.set_yticks(ypos, names)
    ax.set_xlim(0, 520)
    ax.set_ylim(-0.6, len(names) - 0.05)
    ax.text(
        74,
        len(names) - 0.52,
        "flagged",
        ha="center",
        va="bottom",
        fontsize=8.8,
        color=INK2,
    )
    ax.text(
        205,
        len(names) - 0.52,
        "unflagged",
        ha="center",
        va="bottom",
        fontsize=8.8,
        color=INK2,
    )
    ax.text(330, 1, "chance alone: ≈ 0", va="center", fontsize=8.8, color=INK2)
    ax.set_title("b  Synonymous outliers (|z| > 8)")
    ax.set_xlabel(
        "transcripts of 18,406; gray = genes with a release flag,\ncolor = unflagged genes (detail in the outlier figure)"
    )
    style_axis(ax)

    # (c) Calibration slope.
    ax = axes[1, 0]
    snames = names[1:]
    svals = [0.28, 0.75, 0.68, 0.60, rec1.slope]
    sy = list(range(len(snames)))[::-1]
    for y, v, c, h in zip(sy, svals, colors[1:], hatches[1:]):
        ax.barh(y, v, height=0.56, color=c, hatch=h, hatchcolor="white", linewidth=0)
        ax.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=9.5)
    ax.axvline(1.0, **THEORY)
    ax.text(
        0.98,
        len(snames) - 0.45,
        "1 = gnomAD shows all of\nrr's predicted variation",
        ha="right",
        fontsize=8.8,
        color=INK2,
        va="bottom",
    )
    ax.set_yticks(sy, snames)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(-0.6, len(snames) + 0.2)
    ax.set_title("c  Calibration slope (below 1 = rr too strong)")
    ax.set_xlabel("log-log slope of uncorrected obs/exp on gene rr")
    style_axis(ax)

    # (d) Strength.
    ax = axes[1, 1]
    alpha = [0, 0.25, 0.5, 0.625, 0.75, 0.875, 1, 1.25]
    curves = [
        (
            "Legacy syn_adj_r",
            LEGACY_C,
            "-",
            [2.908, 2.675, 2.522, 2.479, 2.461, 2.471, 2.509, 2.683],
        ),
        (
            "Ruchit, 17k DNMs",
            RUCHIT_C,
            "-",
            [2.908, 2.671, 2.521, 2.485, 2.476, 2.500, 2.557, 2.784],
        ),
        (
            "Ruchit, Satterstrom DNMs",
            RUCHIT_C,
            (0, (4, 2.5)),
            [2.908, 2.656, 2.526, 2.512, 2.534, 2.594, 2.696, 3.032],
        ),
    ]
    ax.axvline(1.0, color=AXIS, lw=1)
    ax.text(1.02, 3.0, "in use\n(α = 1)", fontsize=8.8, color=INK2, va="bottom")
    # Theory: same starting noise and same best noise as legacy, but bottoming out at
    # alpha = 1, i.e. what a correctly scaled rr would look like.
    grid = [i / 100 for i in range(0, 126)]
    best, start = 2.461, 2.908
    ax.plot(
        grid,
        [best + (start - best) * (a - 1) ** 2 for a in grid],
        label="Perfectly scaled rr (theory)",
        **THEORY,
    )
    for (label, c, ls, v), slope in zip(curves, (0.75, 0.68, 0.60)):
        ax.plot(alpha, v, color=c, lw=2, ls=ls, label=label)
        i = v.index(min(v))
        ax.scatter(
            [alpha[i]],
            [v[i]],
            s=40,
            color=c,
            edgecolor=SURFACE,
            linewidth=1.5,
            zorder=3,
        )
        ax.scatter(
            [slope],
            [2.415],
            marker="^",
            s=55,
            color=SURFACE if ls != "-" else c,
            edgecolor=c,
            linewidth=1.5,
            zorder=3,
            clip_on=False,
        )
    ax.plot(
        rec.alpha, rec.core_noise, color=REC_C, lw=2.8, label="Ridge, Satterstrom DNMs"
    )
    i = rec.core_noise.idxmin()
    ax.scatter(
        [rec.alpha[i]],
        [rec.core_noise[i]],
        s=46,
        marker="D",
        color=REC_C,
        edgecolor=SURFACE,
        linewidth=1.5,
        zorder=4,
    )
    ax.scatter(
        [rec1.slope],
        [2.415],
        marker="^",
        s=55,
        color=REC_C,
        edgecolor=REC_C,
        zorder=3,
        clip_on=False,
    )
    ax.text(0.02, 2.43, "▲ = slope from panel c", fontsize=8.8, color=INK2, va="bottom")
    ax.set_xlim(-0.03, 1.3)
    ax.set_ylim(2.4, 3.08)
    ax.set_title("d  Every rr is too strong at full strength")
    ax.set_xlabel("strength α (rr^α; 0 = no correction); dots = best α")
    ax.set_ylabel("per-gene noise (mean z²)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.13, 1.0), frameon=False, fontsize=8.8)
    style_axis(ax, grid_axis="y")

    # (e) Best strength vs calibration slope (Ruchit's v2 sections 4 and 9).
    ax = axes[2, 0]
    ax.plot([0, 1.05], [0, 1.05], **THEORY)
    ax.text(
        0.40,
        0.45,
        "best α = slope",
        rotation=45,
        ha="center",
        va="center",
        fontsize=8.8,
        color=INK2,
    )
    pts = [
        ("Gnocchi adj_r", 0.28, 0.25, ADJ_R_C, "o", True, (8, -4)),
        ("Legacy syn_adj_r", 0.75, 0.75, LEGACY_C, "o", True, (10, -16)),
        ("Ruchit 17k", 0.68, 0.75, RUCHIT_C, "o", True, (-72, 6)),
        ("Ruchit 17k, 1e-4", 0.84, 0.875, RUCHIT_C, "s", True, (8, -26)),
        ("Ruchit Satterstrom", 0.60, 0.625, RUCHIT_C, "o", False, (8, -10)),
        ("Ruchit Satterstrom, 1e-4", 0.82, 0.875, RUCHIT_C, "s", False, (-150, 4)),
        (
            "Ridge, Satterstrom",
            rec1.slope,
            float(rec.alpha[rec.core_noise.idxmin()]),
            REC_C,
            "D",
            True,
            (12, -1),
        ),
    ]
    for label, x, y, c, m, filled, off in pts:
        ax.scatter(
            [x],
            [y],
            s=70,
            marker=m,
            color=c if filled else SURFACE,
            edgecolor=c,
            linewidth=1.8,
            zorder=3,
        )
        ax.annotate(
            label,
            (x, y),
            xytext=off,
            textcoords="offset points",
            fontsize=8.8,
            color=INK2,
        )
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("e  Panels c and d agree: best α ≈ slope")
    ax.set_xlabel("calibration slope (panel c)")
    ax.set_ylabel("best strength α (panel d)")
    style_axis(ax, grid_axis="both")

    # (f) chrX/Y.
    ax = axes[2, 1]
    xnames = [
        "No correction",
        "Legacy syn_adj_r\n(no rr on X/Y)",
        "Ruchit 17k, unmasked\n(3 of 13 contexts)",
        "Ruchit 17k, masked",
        "Ridge, Satterstrom\n(no rr on X/Y)",
    ]
    xvals = [4.364, 4.364, 4.332, 4.435, rec1.xy_noise]
    xcols = [NONE_C, LEGACY_C, RUCHIT_C, RUCHIT_C, REC_C]
    xy = list(range(len(xnames)))[::-1]
    ax.axvline(4.364, color=MUTED, lw=1, ls=(0, (3, 3)))
    # Reference: the same relative gain rr gives on all transcripts (5.081 -> 4.447).
    full = 4.364 * 4.447 / 5.081
    ax.axvline(full, **THEORY)
    ax.text(
        full + 0.01,
        len(xnames) - 0.7,
        "if rr helped X/Y as much\nas autosomes (−12%)",
        fontsize=8.8,
        color=INK2,
        va="center",
    )
    for y, v, c in zip(xy, xvals, xcols):
        ax.scatter([v], [y], s=70, color=c, edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.text(v + 0.006, y + 0.22, f"{v:.3f}", fontsize=9.5)
    ax.set_yticks(xy, xnames)
    ax.set_xlim(3.78, 4.48)
    ax.set_ylim(-0.6, len(xnames) - 0.4)
    ax.set_title("f  chrX/Y noise (769 transcripts)")
    ax.set_xlabel("mean z², overall level removed (lower = better)")
    style_axis(ax)

    fig.savefig(os.path.join(OUT_DIR, "rr_performance.png"), dpi=170)
    plt.close(fig)


def draw_deciles(
    summary="syn_oe_by_decile.rr_eval_export.tsv", out="rr_by_decile.png", note=None
):
    """
    Draw syn obs/exp by LOEUF and s_het decile from Ruchit's evaluation export.

    :param summary: Decile summary TSV (remap_rr_to_shet.py --summarize or
        decile_explained_removed.py).
    :param out: Output PNG name.
    :param note: Optional line drawn above the panels.
    """
    df = pd.read_csv(os.path.join(OUT_DIR, summary), sep="\t")
    n = int(df[(df.binning == "shet") & (df.correction == "none")].n.sum())

    def oe(binning, correction):
        d = df[(df.binning == binning) & (df.correction == correction)]
        return d.sort_values("decile").oe.tolist()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.8))
    fig.subplots_adjust(
        left=0.08, right=0.97, top=0.93, bottom=0.08, hspace=0.42, wspace=0.28
    )
    dec = list(range(10))
    series = [
        ("No correction", NONE_C, "-", "none"),
        ("Gnocchi adj_r", ADJ_R_C, "-", "adj_r"),
        ("Legacy syn_adj_r (not normalized)", LEGACY_C, "-", "syn_adj_r"),
        ("Ruchit, 17k DNMs", RUCHIT_C, "-", "ruchit_17k"),
        ("Ruchit, Satterstrom DNMs", RUCHIT_C, (0, (4, 2.5)), "ruchit_satterstrom"),
        ("Ridge, Satterstrom DNMs", REC_C, "-", "recommended"),
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
        ax.text(9.35, 0.997, "expected", fontsize=8.8, color=INK2, va="top")
        for label, c, ls, corr in series:
            v = oe(binning, corr)
            ax.plot(
                dec,
                v,
                color=c,
                lw=2.8 if corr == "recommended" else 2,
                ls=ls,
                label=label,
            )
            ax.scatter(
                [9], [v[-1]], s=36, color=c, edgecolor=SURFACE, linewidth=1.5, zorder=3
            )
        ax.set_xticks(dec)
        ax.set_xlim(-0.3, 10.6)
        ax.set_ylim(0.84, 1.13)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        style_axis(ax, grid_axis="y")
    axes[0, 0].legend(loc="lower left", frameon=False, fontsize=8.5, ncol=2)
    axes[0, 0].annotate(
        "normalized vs not:\nabout 1.3% apart in every decile",
        xy=(3, 0.997),
        xytext=(2.2, 0.950),
        fontsize=8.8,
        color=INK2,
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8),
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
                ("ruchit_17k", RUCHIT_C, "Ruchit 17k"),
                ("recommended", REC_C, "Ridge, Satterstrom"),
            ],
            "j  With rr: LOEUF vs s_het binning",
        ),
    ):
        ax.axhline(1.0, **THEORY)
        ax.text(0, 0.996, "expected = 1", fontsize=8.8, color=INK2, va="top")
        for corr, c, name in corrs:
            lw = 2.8 if corr == "recommended" else 2
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
            ax.text(9.2, loeuf[-1], "by LOEUF", va="center", fontsize=9.5, color=INK2)
            ax.text(
                9.2, shet[-1] - 0.004, "by s_het", va="center", fontsize=9.5, color=INK2
            )
        else:
            ax.legend(
                loc="upper left",
                frameon=False,
                fontsize=8.5,
                bbox_to_anchor=(0.0, 0.97),
            )
        ax.set_xticks(dec)
        ax.set_xlim(-0.3, 11.4)
        ax.set_ylim(0.97, 1.13)
        ax.set_title(title)
        ax.set_xlabel(f"decile; 0 = most constrained (same {n:,} transcripts)")
        style_axis(ax, grid_axis="y")

    if note:
        fig.text(0.08, 0.985, note, fontsize=10, color=INK2, va="top")
    fig.savefig(os.path.join(OUT_DIR, out), dpi=170)
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
                        fontsize=9,
                        color=INK,
                    )
                if unflagged and n >= 25:
                    ax.text(
                        sign * (start + n / 2),
                        y,
                        f"{n}",
                        ha="center",
                        va="center",
                        fontsize=9.5,
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
                fontsize=9.5,
                color=INK,
            )
    ax.axvline(0, color=INK2, lw=1)
    ax.text(
        -10,
        len(rows) - 0.38,
        "← obs ≪ exp (deficit)",
        ha="right",
        fontsize=9.5,
        color=INK2,
    )
    ax.text(
        10,
        len(rows) - 0.38,
        "obs ≫ exp (excess) →",
        ha="left",
        fontsize=9.5,
        color=INK2,
    )
    ax.text(175, -0.45, "chance alone: ≈ 0", ha="right", fontsize=8.8, color=INK2)
    ax.set_yticks(ypos, [r[0] for r in rows])
    ax.set_xlim(-440, 180)
    ax.set_ylim(-0.6, len(rows) - 0.05)
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
        fontsize=9,
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
    ax.text(8.15, 50, "outlier cutoff\n(z = 8)", fontsize=8.8, color=INK2, va="top")
    ax.axvline(0, color=AXIS, lw=1)
    ax.set_xlim(-4, 18)
    ax.set_ylim(0, 52)
    ax.set_xticks(range(-4, 19, 2))
    ax.set_title("b  Synonymous z of histone genes (no correction)")
    ax.set_xlabel(
        "synonymous z = (obs − exp) / √exp\npositive = more variants than expected"
    )
    ax.set_ylabel("% of genes in the group")
    leg = ax.legend(
        loc="upper left", fontsize=9, bbox_to_anchor=(0.0, 1.02), framealpha=1
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
            s=24,
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
                s=24,
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
        ax.text(label_x, 2.95, curve_label, fontsize=8.5, color=INK2, va="top")
        for gene, dx, dy in (("CSMD1", -40, 10), ("H1-4", 8, 6), ("MUC4", 8, -12)):
            r = zs[zs.gene == gene]
            if len(r):
                r = r.iloc[0]
                ax.annotate(
                    gene,
                    (r[xcol], r.oe),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=8.5,
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
        fontsize=9,
        bbox_to_anchor=(0.54, 0.0),
        markerscale=1.4,
    )

    fig.savefig(os.path.join(OUT_DIR, "rr_outliers.png"), dpi=170)
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
            s=70,
            marker=m,
            color=c,
            edgecolor=SURFACE,
            linewidth=1.5,
            zorder=3,
        )
        ax.scatter(
            [r.crossfit_noise],
            [y],
            s=60,
            marker=m,
            color=SURFACE,
            edgecolor=c,
            linewidth=1.8,
            zorder=3,
        )
        axes[1].scatter(
            [r.slope],
            [y],
            s=70,
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
            fontsize=9,
            color=INK2,
        )
    ax = axes[0]
    ax.axvline(none.core_noise, **THEORY)
    ax.text(
        none.core_noise - 0.01,
        len(rows) - 0.35,
        "no correction",
        ha="right",
        fontsize=8.8,
        color=INK2,
    )
    ax.set_yticks(ypos, labels)
    ax.set_ylim(-0.7, len(rows) - 0.1)
    ax.set_xlim(2.3, 2.95)
    ax.set_title("Per-gene synonymous noise, held-out chromosomes")
    ax.set_xlabel(
        "mean z², overall level removed (lower = better)\nfilled = rr as fitted; open = with the best α added"
    )
    style_axis(ax)
    ax = axes[1]
    ax.axvline(1.0, **THEORY)
    ax.text(1.02, len(rows) - 0.35, "1 = no α needed", fontsize=8.8, color=INK2)
    ax.set_xlim(0.3, 1.6)
    ax.set_title("Calibration slope")
    ax.set_xlabel("below 1 = too strong; above 1 = too weak")
    style_axis(ax)
    fig.text(
        0.01,
        0.02,
        "* background built from the exported MANE Select cells (approximation); compare the two rows with each other only",
        fontsize=8.5,
        color=INK2,
    )
    fig.savefig(os.path.join(OUT_DIR, "rr_alternatives.png"), dpi=170)
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
    ax.plot(x, size.none, color=NONE_C, lw=2, marker="o", ms=5, label="No correction")
    ax.plot(x, size.ridge_rr, color=REC_C, lw=2.6, marker="D", ms=5, label="Ridge rr")
    ax.set_xscale("log")
    ax.set_xticks(
        [40, 80, 150, 300, 500, 1000], ["40", "80", "150", "300", "500", "1,000"]
    )
    ax.set_ylim(0, 6.2)
    ax.set_title("Leftover noise grows with gene size")
    ax.set_xlabel("expected synonymous variants per gene (bin median)")
    ax.set_ylabel("mean z²")
    ax.legend(loc="upper left", frameon=False, fontsize=8.8)
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
        ax.text(v + 0.03, y, f"{v:.2f}", va="center", fontsize=9.5)
    ax.axvline(CHANCE_FLOOR, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.04,
        len(order) - 0.45,
        "chance ≤ 0.70",
        fontsize=8.8,
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

    fig.savefig(os.path.join(OUT_DIR, "rr_leftover_noise.png"), dpi=170)
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
    fig.subplots_adjust(left=0.07, right=0.97, top=0.96, bottom=0.06)

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
            (12, -6),
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
    ax.text(1.93, 1.83, "DNMs agree", fontsize=8.8, color=INK2, ha="right", rotation=36)
    ax.axhline(1, color=AXIS, lw=1)
    ax.text(
        1.98,
        0.97,
        "DNMs show no excess",
        fontsize=8.8,
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
            ms=9,
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
                xytext=(1.47, 1.72),
                fontsize=8.8,
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
            fontsize=8.8,
            color=INK,
            va="center",
            ha="right" if dx < 0 else "left",
        )
    ax.set_xlim(0.9, 2.0)
    ax.set_ylim(0.7, 3.7)
    ax.set_title("a  De novo mutations show the same excess, except in repeats")
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
        fontsize=8.8,
        color=INK2,
        ha="right",
        va="bottom",
    )
    ax.text(
        1.235,
        0.975,
        "if none of it were",
        fontsize=8.8,
        color=INK2,
        ha="right",
        va="top",
    )
    ax.plot(x, 1 + 0.51 * (x - 1), color=REC_C, lw=1.6)
    ax.text(
        0.9,
        0.885,
        "solid line: fit, slope 0.51 ± 0.07",
        fontsize=9,
        color=INK,
        va="top",
    )
    ax.errorbar(
        dec.gnomad_oe,
        dec.dnm_oe,
        yerr=[dec.dnm_oe - dec.dnm_lo, dec.dnm_hi - dec.dnm_oe],
        fmt="o",
        ms=7,
        color=REC_C,
        mec=SURFACE,
        mew=1.2,
        elinewidth=1.6,
        capsize=0,
        zorder=3,
    )
    ax.set_xlim(0.82, 1.24)
    ax.set_ylim(0.82, 1.2)
    ax.set_title("b  Across all genes, DNMs follow synonymous obs/exp")
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
    ax.scatter(h.x, h.z.clip(-12, 24), s=16, color=histone_c, lw=0, zorder=2)
    for key, (c, _) in group_style.items():
        g = genes[genes.outlier & (genes.group == key)]
        ax.scatter(
            g.x, g.z.clip(-12, 24), s=30, color=c, edgecolor=SURFACE, lw=0.7, zorder=3
        )
    ax.axhline(8, **THEORY)
    ax.axhline(-8, **THEORY)
    ax.text(offset["chr4"] + 4e7, 8.6, "|z| = 8", fontsize=8.8, color=INK2)
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
                fontsize=8.5,
                color=INK2,
                ha="right" if dx < 0 else "left",
            )
    centers = offset + length / 2
    ax.set_xticks(centers.values, [c.replace("chr", "") for c in chroms], fontsize=8.5)
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
        ncol=5,
        frameon=False,
        fontsize=8.8,
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
            r.dnm_oe, y - 0.18, s=55, color=SURFACE, edgecolor=INK2, lw=1.6, zorder=3
        )
        ax.scatter(
            r.gnomad_oe_no_rr,
            y + 0.14,
            s=55,
            color=NONE_C,
            edgecolor=SURFACE,
            lw=1,
            zorder=3,
        )
        ax.scatter(
            r.gnomad_oe,
            y + 0.14,
            s=60,
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
        fontsize=8.8,
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
        fontsize=8.8,
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
                ms=5,
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
                ms=3,
                label=(
                    "one ancestry group each: African, Admixed American,\nEast Asian, European (non-Finnish), South Asian"
                    if anc == "afr"
                    else None
                ),
            )
    ax.legend(loc="upper right", frameon=False, fontsize=8.8)
    ax.set_xscale("log")
    ax.set_xticks([1000, 10000, 100000, 500000], ["1k", "10k", "100k", "500k"])
    ax.set_xlim(800, 8e5)
    ax.set_ylim(0.9, 1.8)
    ax.set_title("e  Same excess in every ancestry group and sample size")
    ax.set_xlabel("downsampled exomes (log scale)")
    ax.set_ylabel(
        "excess outliers' synonymous obs/exp,\nrelative to genes without a flag"
    )
    style_axis(ax, grid_axis="y")

    fig.savefig(os.path.join(OUT_DIR, "rr_outlier_drivers.png"), dpi=170)
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
                        fontsize=7.2,
                        color=MUTED,
                    )
                elif j == 4 or abs(m[i, j]) >= 0.5 * vmax:
                    ax.text(
                        j,
                        i,
                        f"{m[i, j]:+.1f}",
                        ha="center",
                        va="center",
                        fontsize=7.8,
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
        "% change in DNM rate per SD of the feature\n(red = more mutable)", fontsize=9
    )
    fig.suptitle(
        "What each feature does in the ridge rr: % change in rate per SD, by window scale ('all' = sum)",
        x=0.17,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.text(
        0.17,
        0.04,
        "Satterstrom DNMs against the absent-from-gnomAD background (exported MANE Select cells). "
        "Correlated features share signal, so read directions, not single values.",
        fontsize=8.8,
        color=INK2,
    )
    fig.savefig(os.path.join(OUT_DIR, "rr_features.png"), dpi=170)
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
            s=70,
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
            s=70,
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
            fontsize=9.5,
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
        fontsize=9.5,
    )
    ax.set_xlabel(
        "per-gene synonymous noise on held-out chromosomes, mean z² (lower = better)"
    )
    style_axis(ax)
    fig.savefig(os.path.join(OUT_DIR, "rr_dnm_sets.png"), dpi=170)
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
        fontsize=8.8,
        color=INK2,
        ha="left",
        va="top",
    )
    ax.text(
        1.9e5,
        groups.loc["all excess outliers"].dnm_lo + 0.04,
        "DNM obs/exp, the 46 non-histone outliers",
        fontsize=8.8,
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
                ms=5,
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
                ms=2.5,
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
        ms=4,
        label="46 non-histone outliers",
    )
    ax.axhline(1, color=AXIS, lw=1)
    ax.set_xscale("log")
    ax.set_xticks([1000, 10000, 100000, 500000], ["1k", "10k", "100k", "500k"])
    ax.set_xlim(800, 7e5)
    ax.set_ylim(0.6, 3.7)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    ax.set_title(
        "a  The histone excess is about 3× while few sites are seen,\n    and shrinks as sites saturate"
    )
    ax.set_xlabel("downsampled exomes (log scale)")
    ax.set_ylabel("synonymous obs/exp (no rr),\nrelative to genes without a flag")
    style_axis(ax, grid_axis="y")

    ax = axes[1]
    labels = [
        "no correction",
        "× 2.84 (DNM excess),\nmultiplied",
        "× 2.84, before\nsaturation",
        "× 5.2 non-CpG,\n× 1.15 CpG, before\nsaturation",
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
    ax.text(len(fac) - 0.55, 1.03, "obs/exp = 1", fontsize=8.8, color=INK2, ha="right")
    ax.set_xticks(x, labels, fontsize=8.8)
    ax.set_ylim(0, 2.3)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title("b  A histone factor learned from DNMs\n    would overcorrect gnomAD")
    ax.set_ylabel("histone cluster synonymous obs/exp")
    style_axis(ax, grid_axis="y")

    fig.savefig(os.path.join(OUT_DIR, "rr_histones.png"), dpi=170)
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
        ["rr for every 1-kb window and context", "about 0.9 to 1.1"],
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

    fig.savefig(os.path.join(OUT_DIR, "rr_concept.png"), dpi=170)
    plt.close(fig)


def draw_model():
    """Draw how the rr model is set up: training rows, per-context fit, and rr."""
    W, H = 1560, 600
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
        "Each row carries the features of its 1-kb window:\n13 features × 4 window sizes = 52 columns.",
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
            "no screen",
            "all 52 columns",
            "penalty shrinks every β",
            "toward 0 (none to exactly 0)",
        ],
        size=11,
    )
    ax.text(
        790,
        410,
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
    b0, d = -3.0, 0.35
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
        "real values are about 0.9–1.1",
        fontsize=11,
        color=INK2,
        va="top",
    )
    ax.text(
        1090,
        526,
        "Then: normalize so rr averages 1 over exome synonymous\nsites in each context, and multiply each variant's expected\n"
        "count by its window's rr.",
        fontsize=11,
        color=INK2,
        va="top",
    )

    fig.savefig(os.path.join(OUT_DIR, "rr_model.png"), dpi=170)
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
        1, 3, figsize=(16, 4.9), gridspec_kw={"width_ratios": [1, 1, 1.1]}
    )
    fig.subplots_adjust(left=0.05, right=0.98, top=0.86, bottom=0.17, wspace=0.42)

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
        fontsize=9,
        color=INK2,
    )
    ax.text(
        2.6,
        top * 0.55,
        f"bars: observed\n(mean z² = {np.mean(z ** 2):.1f})",
        fontsize=9,
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
        ms=5,
        label="no correction",
    )
    ax.plot(
        g.rr,
        g.obs / g.e1 / k1,
        color=REC_C,
        lw=2.6,
        marker="D",
        ms=5,
        label="with the ridge rr",
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.legend(loc="upper left", frameon=False, fontsize=9.5)
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
        ("Gnocchi adj_r (non-coding)", 5.450, ADJ_R_C),
        ("Legacy syn_adj_r", 2.509, LEGACY_C),
        (
            "Ridge rr, Satterstrom DNMs",
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
            fontsize=9.5,
        )
    ax.axvline(CHANCE_FLOOR, **THEORY)
    ax.text(
        CHANCE_FLOOR + 0.04,
        len(rows) - 0.55,
        f"chance alone ≤ {CHANCE_FLOOR:.2f}",
        fontsize=8.8,
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
        "left after ridge: mostly real\ngene-level mutation-rate\ndifferences (DNMs agree)",
        ha="left",
        va="center",
        fontsize=8.6,
        color=INK2,
    )
    ax.set_yticks(ypos, [r[0] for r in rows])
    ax.set_xlim(0, 4.2)
    ax.set_ylim(-1.0, len(rows) - 0.3)
    ax.set_title("c  What each regional correction removes")
    ax.set_xlabel("per-gene synonymous noise, mean z² (lower = better)")
    style_axis(ax)

    fig.savefig(os.path.join(OUT_DIR, "rr_problem_summary.png"), dpi=170)
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
            g.rr_all_sites, 100 * g.frac_seen, color=class_c[c], lw=2, marker="o", ms=5
        )
        ax.text(
            g.rr_all_sites.iloc[-1] + 0.006,
            100 * g.frac_seen.iloc[-1],
            c,
            color=INK,
            fontsize=9.5,
            va="center",
        )
    ax.text(
        0.02,
        0.97,
        "dotted: if sites were seen at a rate\nproportional to rr alone",
        transform=ax.transAxes,
        va="top",
        fontsize=8.8,
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
    ax.text(
        1.2, 1.12, "dotted: same rr", ha="right", va="top", fontsize=8.8, color=INK2
    )
    for c, g in pairs.groupby("class"):
        b = np.polyfit(g.log_rr_all_sites, g.log_rr_absent, 1)[0]
        ax.plot(
            np.exp(g.log_rr_all_sites),
            np.exp(g.log_rr_absent),
            color=class_c[c],
            lw=2,
            marker="o",
            ms=4,
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
    ax.legend(loc="upper left", frameon=False, fontsize=9.5)
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
            r.oof_noise, y, s=70, color=REC_C, edgecolor=SURFACE, lw=1.2, zorder=3
        )
        txt = f"{r.oof_noise:.2f}"
        pc = met[(met.background == bg) & (met.application == "per-class strength")]
        if len(pc):
            v = pc.oof_noise.iloc[0]
            ax.plot([v, r.oof_noise], [y, y], color=REC_C, lw=1.2)
            ax.scatter(
                v, y, s=70, marker="D", color=SURFACE, edgecolor=REC_C, lw=1.6, zorder=3
            )
            txt += f" → {v:.2f}"
        ax.text(
            max(r.oof_noise, none.core_noise * 0 + r.oof_noise) + 0.012,
            y + 0.18,
            txt,
            fontsize=9,
            va="bottom",
        )
    ax.axvline(none.core_noise, **THEORY)
    ax.text(
        none.core_noise - 0.008,
        len(backgrounds) - 0.45,
        "no correction",
        ha="right",
        fontsize=8.8,
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
                s=45,
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
                    fontsize=8.5,
                    color=INK2,
                )
    ax.axvline(1, **THEORY)
    ax.text(
        1.03, len(backgrounds) - 0.45, "1 = right strength", fontsize=8.8, color=INK2
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
        fontsize=9.5,
    )
    style_axis(ax)

    fig.savefig(os.path.join(OUT_DIR, "rr_background.png"), dpi=170)
    plt.close(fig)


LOGIT_EXAMPLES = ["CAG", "CCT", "CCG"]
"""The context with the most DNMs in each class."""

LOGIT_BACKGROUNDS = [
    ("absent", "Absent from gnomAD (current rule)", REC_C, "-"),
    (
        "absent_context_matched",
        "Absent, sampled by context in proportion to all sites",
        "#2a78d6",
        (0, (5, 2)),
    ),
    (
        "absent_equal_per_context",
        "Absent, sampled evenly by context",
        "#1baf7a",
        (0, (1, 1.5)),
    ),
    ("all_sites", "All possible synonymous variants", "#eb6834", "-"),
]


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
        ins.tick_params(labelsize=7.5, length=0, labelcolor=MUTED)
        ins.set_facecolor(SURFACE)
        for side in ins.spines.values():
            side.set_color(AXIS)
        ins.set_title(
            "whole curve; box = this panel", fontsize=8, color=INK2, weight="normal"
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
        s=12 if small else 32,
        color=SURFACE,
        edgecolor=INK,
        lw=0.9 if small else 1.2,
        zorder=3,
    )
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.2, 1.2)
    ax.set_yticks([0, 0.5, 1])
    style_axis(ax, grid_axis="y")
    if small:
        ax.set_title(
            f"{context}  ({row.dnms:.0f})", fontsize=9.5, color=INK, weight="bold"
        )
        ax.tick_params(labelsize=8)
        return
    ax.text(
        -xmax + 0.05,
        1.16,
        "DNM sites (y = 1)",
        fontsize=8.8,
        color=INK2,
        va="bottom",
    )
    ax.text(
        -xmax + 0.05,
        -0.16,
        "background sites (y = 0)",
        fontsize=8.8,
        color=INK2,
        va="top",
    )
    ax.text(
        min(hi + 0.04, xmax - 0.27),
        0.2,
        f"shaded: middle 90%\nof sites, rr {np.exp(lo):.2f}–{np.exp(hi):.2f}",
        fontsize=8.8,
        color=INK2,
    )


def rr_curve_panel(ax, context, bins, ctx):
    """Draw one context's DNM rate per score decile against the model's rr."""
    r = ctx.set_index("context").loc[context]
    b = bins[(bins.context == context) & (bins.background == "absent")]
    n0 = b.background_sites
    k = b.dnms.sum() / n0.sum()
    rel = b.dnms / n0 / k
    se = rel * np.sqrt(1 / b.dnms.clip(lower=1) + 1 / n0)
    s = np.linspace(b.score.min(), b.score.max(), 100)
    ax.plot(s, np.exp(s), color=INK, lw=2.4, label="model rr (≈ eˢ)")
    ax.errorbar(
        b.score,
        rel,
        yerr=1.96 * se,
        fmt="o",
        color=CLASS_C[r["class"]],
        mfc=SURFACE,
        mew=1.4,
        ms=7,
        elinewidth=1.1,
        capsize=0,
        label="DNMs per background site, relative",
    )
    ax.axhline(1, color=AXIS, lw=1, zorder=0)
    ax.set_yscale("log")
    ticks = [0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6]
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(0.6, 1.75)
    style_axis(ax, grid_axis="y")


def score_dist_panel(ax, context, hist, ctx):
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
            lw=2.2 if y else 1.8,
            label=label,
        )
        ax.axvline(
            r.mean_dnm if y else r.mean_background,
            color=color,
            lw=1.4,
            ls=(0, (4, 2)),
        )
    shift = (r.mean_dnm - r.mean_background) / r.sd_background
    ax.text(
        0.0,
        1.02,
        f"means (dashed) differ by {r.mean_dnm - r.mean_background:+.3f}"
        f" = {shift:.2f} SD of the background scores",
        transform=ax.transAxes,
        va="bottom",
        fontsize=9,
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
        fontsize=13,
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
            ax.set_ylabel(
                "chance a site is a DNM\n(drawn with equal numbers of DNMs and"
                " background)"
            )
    note_y = ax.get_position().y0 - 0.045
    fig.text(
        0.07,
        note_y,
        "solid: the fitted curve where sites fall · dashed: the same curve beyond them"
        " · open circles: share of DNMs in each score decile, on held-out chromosomes",
        fontsize=9,
        color=INK2,
    )
    if stacked:
        for i, context in enumerate(LOGIT_EXAMPLES):
            ax = fig.add_subplot(gs[1, i])
            rr_curve_panel(ax, context, bins, ctx)
            ax.set_title(f"{'def'[i]}  {context}: DNMs against rr, per score decile")
            ax.set_xlabel("regional score s (decile mean, held-out)")
            if i == 0:
                ax.set_ylabel(
                    "relative DNM rate\n(DNMs per background site ÷ context" " average)"
                )
                ax.legend(loc="upper left", frameon=False, fontsize=9)
            ax = fig.add_subplot(gs[2, i])
            score_dist_panel(ax, context, hist, ctx)
            ax.set_title(f"{'ghi'[i]}  {context}: score distributions", pad=22)
            ax.set_xlabel("regional score s (held-out)")
            if i == 0:
                ax.set_ylabel("density")
                ax.legend(loc="center left", frameon=False, fontsize=9)

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
        ax.text(xs.min() - 0.3, 1.49, c, fontsize=10, color=INK, weight="bold")
    ax.axhline(1, color=AXIS, lw=1, zorder=0)
    ax.set_xticks(
        x, [f"{r.context} ({r.dnms:.0f})" for _, r in ctx.iterrows()], rotation=90
    )
    ax.tick_params(axis="x", labelsize=8.8, labelcolor=INK2)
    ax.set_ylim(0.75, 1.55)
    ax.set_xlim(x.min() - 0.8, x.max() + 0.8)
    ax.set_ylabel(f"rr of the middle 90%\nof {sites.split()[0]} sites")
    ax.set_title(
        f"{'j' if stacked else 'd'}  Every context's model: the range of rr it gives"
        f" the middle 90% of {sites} (DNMs per context in brackets)"
    )
    style_axis(ax, grid_axis="y")
    fig.savefig(os.path.join(OUT_DIR, f"{prefix}.png"), dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(4, 8, figsize=(17, 9.5), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.05, right=0.99, top=0.9, bottom=0.08, hspace=0.45, wspace=0.12
    )
    fig.suptitle(
        "Every context's logistic regression (DNMs in brackets; colour = context"
        " class; drawn with equal numbers of DNMs and background)",
        x=0.05,
        ha="left",
        fontsize=12.5,
        weight="bold",
    )
    for ax, (_, r) in zip(axes.ravel(), ctx.iterrows()):
        logistic_panel(ax, r.context, bins, points, ctx, rng, small=True)
    fig.supxlabel("regional score s = β·z (log-odds; rr ≈ eˢ)", fontsize=10.5)
    fig.supylabel("chance a site is a DNM", fontsize=10.5)
    fig.savefig(os.path.join(OUT_DIR, f"{prefix}.all_contexts.png"), dpi=150)
    plt.close(fig)


def draw_logistic_backgrounds():
    """Draw how the logistic fits change with the training background."""
    bins, _, ctx = read_logistic_fits()
    cls = ctx.set_index("context")["class"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.subplots_adjust(
        left=0.08, right=0.98, top=0.87, bottom=0.08, hspace=0.42, wspace=0.2
    )
    fig.suptitle(
        "The same logistic regressions under each training background",
        x=0.08,
        ha="left",
        fontsize=13,
        weight="bold",
    )
    fig.text(
        0.08,
        0.925,
        "a–c: lines = fitted curves, dots = observed share of DNMs per score decile."
        "   d–f: lines = rr under each background; circles = DNMs per possible site"
        " in each decile, held-out, ± 95% CI; dotted = rr equal to the DNM rate.",
        fontsize=9.2,
        color=INK2,
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
                s=16,
                color=c,
                edgecolor=SURFACE,
                lw=0.6,
                zorder=3,
            )
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.yaxis.set_minor_formatter(plt.FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.tick_params(axis="y", which="both", labelsize=8.8, labelcolor=INK2)
        ax.set_title(f"{'abc'[i]}  {context} ({cls[context]}): the fitted curves")
        ax.set_xlabel("regional score under the current rule's model")
        if i == 0:
            ax.set_ylabel("share of training rows that are DNMs\n(log scale)")
        style_axis(ax, grid_axis="y")
    axes[0, 0].legend(
        loc="upper left",
        bbox_to_anchor=(0, -0.2),
        ncol=4,
        frameon=False,
        fontsize=9.5,
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
            ms=6,
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
        ax.set_title(
            f"{'def'[i]}  {c}, all {len(tot)} contexts: rr against where DNMs fall"
        )
        ax.set_xlabel("rr under the current rule (decile mean, held-out)")
        if i == 0:
            ax.set_ylabel(
                "rr under each background (lines)\nDNMs per possible site, relative"
                " (circles)"
            )
        style_axis(ax, grid_axis="y")
    axes[1, 0].text(
        0.03,
        0.95,
        "the four backgrounds' rr lines\nlie almost on top of each other",
        transform=axes[1, 0].transAxes,
        va="top",
        fontsize=8.8,
        color=INK2,
    )
    fig.savefig(os.path.join(OUT_DIR, "rr_logistic_fits.backgrounds.png"), dpi=170)
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
        fontsize=13,
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
            ax.scatter(r.ratio, yi, s=60, color=c, edgecolor=SURFACE, lw=1.2, zorder=3)
            ax.text(
                r.ratio_hi + xmax * 0.02,
                yi + 0.16,
                f"{r.ratio:.2f} ({r.observed:.0f} / {r.expected:.0f})",
                fontsize=8.8,
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
                s=46,
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
    ax.set_title(
        "c  Within the genes: paralog-copying alleles\n    against others of the same type"
    )
    ax.set_xlabel(
        "odds ratio (95% interval)\nfilled: by trinucleotide · open: + Roulette site rate"
    )
    style_axis(ax)
    fig.savefig(os.path.join(OUT_DIR, "rr_histone_conversion.png"), dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    draw_diagram()
    draw_performance()
    draw_deciles()
    draw_deciles(
        "syn_oe_by_decile.explained_removed.tsv",
        "rr_by_decile.explained_removed.png",
        "Explained synonymous outliers removed: 91 histone, repeat or segdup, and maternal-hotspot genes",
    )
    draw_outliers()
    draw_alternatives()
    draw_leftover_noise()
    draw_outlier_drivers()
    draw_background()
    draw_concept()
    draw_model()
    draw_dnm_sets()
    draw_histones()
    draw_problem_summary()
    draw_features()
    draw_logistic_fits()
    draw_logistic_fits("rr_logistic_fits.pipeline", "background sites")
    draw_logistic_backgrounds()
    draw_histone_conversion()
