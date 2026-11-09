"""Figures for the notes ``notes_one_model_for_expected.md``, linked from the rr report.

``draw_stack``: a graphical abstract of how expected counts are built now, step by
step, with the data each step is fit on and the patches between them. Uses the
plateau training groups (``plateau_models.groups.tsv``), the stored plateau lines,
the v4.1.1 coverage model, the full run's ridge rr cells along chr8, and Roulette's
per-gene regional factor along chr8 (``compare_rate_parts.genes.tsv.gz``).

``draw_fit_results``: how the single cloglog models from ``one_model_fit.py`` compare
with the current stack: per-gene noise, and synonymous obs/exp by coverage and by
context class.

``draw_release_checks``: the same fits applied to missense and pLoF variants: pLoF and
missense obs/exp by LOEUF decile, levels relative to synonymous, LOEUF changes, and NDD
and HI precision-recall (``one_model_release_eval.py``).

``draw_issue_figures``: one figure per problem of the current expected counts, with
every tried fix (``one_model_issues.py``).

``draw_coverage_check``: the stored coverage model against the share of synonymous
variants seen at each exome AN % (``coverage_check.py``).

``draw_fit_checks``: the report's other synonymous checks on the same fits (outliers by
category, calibration slope, obs/exp by LOEUF and s_het decile, maternal hotspot
regions; ``one_model_eval.py``).

``draw_coverage_oe_current``: synonymous obs/exp by exome AN % under the current models
(the released expected counts), as panel b of ``one_model_fit.png`` with only the
current stack, plus each AN % value.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from make_followup_figures import (
    AXIS,
    CLASS_C,
    COVERAGE_MODEL,
    INK,
    INK2,
    LEGACY_C,
    NONE_C,
    OUT_DIR,
    REC_C,
    ROULETTE_C,
    RR_EDGE,
    RR_FILL,
    SUB_C,
    SURFACE,
    arrow,
    draw_disease_gene_pr,
    style_axis,
)
from matplotlib.patches import FancyBboxPatch
from one_model_release_eval import KEYS

PATCH_C = "#c4302b"


def draw_stack() -> None:
    """Draw the current stack as a graphical abstract."""
    g = pd.read_csv(os.path.join(OUT_DIR, "plateau_models.groups.tsv"), sep="\t")
    fits = pd.read_csv(os.path.join(OUT_DIR, "plateau_models.fits.tsv"), sep="\t")
    rr = pd.read_csv(
        os.path.join(OUT_DIR, "data/full_run/syn_cells_rr_full_run.tsv"),
        sep="\t",
        usecols=["win", "context", "rr_ridge_qcexcl"],
    ).drop_duplicates(["win", "context"])
    genes = pd.read_csv(
        os.path.join(OUT_DIR, "compare_rate_parts.genes.tsv.gz"), sep="\t"
    )

    W, H = 1800, 990
    fig = plt.figure(figsize=(W / 100, H / 100))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_xlim(0, W)
    bg.set_ylim(H, 0)
    bg.axis("off")

    def inset(x, y, w, h):
        a = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
        a.set_facecolor("none")
        a.tick_params(labelsize=10.5, length=0)
        return a

    def card(x, y, w, h, title, sub, data, fill="white", edge=AXIS):
        bg.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0,rounding_size=10",
                facecolor=fill,
                edgecolor=edge,
                lw=1.5,
            )
        )
        bg.text(x + 12, y + 22, title, fontsize=14, weight="bold", va="center")
        bg.text(x + 12, y + 46, sub, fontsize=11.5, color=INK2, va="center")
        bg.text(
            x + 12,
            y + h - 20,
            "fit on: " + data,
            fontsize=11,
            color=INK,
            va="center",
            style="italic",
        )

    def patch(x, y, w, text):
        bg.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                74,
                boxstyle="round,pad=0,rounding_size=8",
                facecolor="#fdecea",
                edgecolor=PATCH_C,
                lw=1.2,
            )
        )
        bg.text(
            x + 10,
            y + 37,
            text,
            fontsize=11,
            color=PATCH_C,
            va="center",
            linespacing=1.3,
        )

    bg.text(
        30,
        34,
        "The stack today: five models, each fit on its own data, each correcting what"
        " the others miss",
        fontsize=17,
        weight="bold",
        va="center",
    )
    xs = [30, 385, 740, 1095, 1450]
    cw, ch, y0 = 320, 360, 80

    # 1. Mutation rates.
    card(
        xs[0],
        y0,
        cw,
        ch,
        "1  Mutation rate μ",
        "one per context × ref × alt × methylation",
        "gnomAD genomes, rare intergenic variants",
    )
    a = inset(xs[0] + 50, y0 + 70, cw - 70, ch - 175)
    subs = list(SUB_C)
    for i, s in enumerate(subs):
        d = g[g["sub"] == s].drop_duplicates(
            ["context", "ref", "alt", "methylation_level"]
        )
        a.scatter(
            np.full(len(d), i) + np.random.default_rng(i).uniform(-0.2, 0.2, len(d)),
            d.mu * 1e8,
            s=14,
            color=SUB_C[s],
            alpha=0.8,
            lw=0,
        )
    a.set_yscale("log")
    a.set_xticks(
        range(len(subs)),
        [s.replace("CpG transition", "CpG ts") for s in subs],
        rotation=45,
        ha="right",
        fontsize=9.5,
    )
    a.set_ylabel("μ (× 10⁻⁸)", fontsize=10.5)
    style_axis(a, grid_axis="y")

    # 2. Plateau model.
    card(
        xs[1],
        y0,
        cw,
        ch,
        "2  Plateau model f(μ)",
        "μ → share of possible variants seen",
        "exome synonymous sites, high coverage",
    )
    # One straight line per model group, each on its own linear axes (the two groups
    # sit at very different rates), dots colored by substitution as in card 1.
    gp = g[g.possible >= 100]
    for cpg, x0, w, name, title in (
        (False, xs[1] + 48, 120, "everything else", "everything else"),
        (True, xs[1] + 200, 96, "CpG transition", "CpG transitions"),
    ):
        a = inset(x0, y0 + 112, w, ch - 205)
        d = gp[gp.cpg == cpg]
        for s, c in SUB_C.items():
            e = d[d["sub"] == s]
            a.scatter(
                e.mu * 1e8,
                e.prop,
                s=8 + 70 * np.sqrt(e.possible / gp.possible.max()),
                color=c,
                alpha=0.6,
                lw=0,
                zorder=3,
            )
        st = fits[(fits.fit == "stored") & (fits.line == name)].iloc[0]
        x = np.linspace(0, d.mu.max() * 1.05, 20)
        a.plot(x * 1e8, st.intercept + st.slope * x, color=INK, lw=2, zorder=4)
        a.set_xlim(0, d.mu.max() * 1.05e8)
        a.set_title(title, fontsize=10.5, color=INK2, pad=4)
        a.tick_params(labelsize=9.5)
        style_axis(a, grid_axis="both")
        if not cpg:
            a.set_ylim(0, 0.5)
            a.set_ylabel("share seen", fontsize=10.5)
        else:
            a.set_ylim(0.5, 1.05)
    bg.text(
        xs[1] + cw / 2 + 10,
        y0 + ch - 50,
        "μ (× 10⁻⁸)",
        fontsize=10.5,
        ha="center",
        va="center",
    )
    bg.text(
        xs[1] + 12,
        y0 + 70,
        "a straight line per group, share = a + b·μ",
        fontsize=10.5,
        color=INK,
        va="center",
    )

    # 3. Coverage model.
    card(
        xs[2],
        y0,
        cw,
        ch,
        "3  Coverage model h(AN %)",
        "scales expected down at low coverage",
        "exome synonymous sites, low coverage",
    )
    a = inset(xs[2] + 50, y0 + 70, cw - 70, ch - 160)
    cov = np.linspace(19, 100, 200)
    hcov = np.where(cov >= 90, 1.0, COVERAGE_MODEL[0] + COVERAGE_MODEL[1] * cov)
    a.plot(cov, hcov, color=INK, lw=2.2)
    a.axvline(90, color=AXIS, lw=1, ls=(0, (3, 3)))
    a.axhline(1, color=AXIS, lw=1)
    a.text(91, 0.15, "high coverage:\nh = 1", fontsize=10, color=INK2)
    a.text(62, 1.2, "overshoots 1\nbelow 90", fontsize=10, color=PATCH_C, ha="right")
    a.set_xlabel("exome AN %", fontsize=10.5)
    a.set_ylabel("h", fontsize=10.5)
    a.set_ylim(0, 1.3)
    style_axis(a, grid_axis="both")

    # 4. rr.
    card(
        xs[3],
        y0,
        cw,
        ch,
        "4  Regional rate rr",
        "a logistic model per context on window features",
        "DNMs vs background, 52 features",
        fill=RR_FILL,
        edge=RR_EDGE,
    )
    a = inset(xs[3] + 50, y0 + 70, cw - 70, ch - 160)
    d = rr[(rr.context == "CCT") & rr.win.str.startswith("chr8-")].copy()
    d["pos"] = d.win.str.split("-").str[1].astype(int) / 1e6
    d = d[d.pos < 60].sort_values("pos")
    a.axvspan(0, 20, color=RR_FILL, lw=0)
    a.scatter(d.pos, d.rr_ridge_qcexcl, s=5, color=RR_EDGE, alpha=0.5, lw=0)
    a.plot(
        d.pos,
        d.rr_ridge_qcexcl.rolling(41, center=True, min_periods=5).median(),
        color=INK,
        lw=1.8,
    )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.set_xlabel("chr8 position (Mb), context CCT", fontsize=10.5)
    a.set_ylabel("rr per 1-kb window", fontsize=10.5)
    a.set_ylim(0.8, 1.25)
    style_axis(a, grid_axis="y")

    # 5. Per-variant rates.
    card(
        xs[4],
        y0,
        cw,
        ch,
        "5  Per-variant rate m (trial)",
        "Roulette: sequence and regional parts",
        "gnomAD v3 genomes, rare variants",
        edge=ROULETTE_C,
    )
    a = inset(xs[4] + 50, y0 + 70, cw - 70, ch - 160)
    c8 = genes[(genes.chromosome == "chr8") & (genes.mid < 60e6)].sort_values("mid")
    a.axvspan(0, 20, color=RR_FILL, lw=0)
    a.scatter(
        c8.mid / 1e6, c8.roulette_multiplier, s=7, color=ROULETTE_C, alpha=0.5, lw=0
    )
    a.plot(
        c8.mid / 1e6,
        c8.roulette_multiplier.rolling(21, center=True, min_periods=5).median(),
        color=INK,
        lw=1.8,
    )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.set_xlabel("chr8 position (Mb)", fontsize=10.5)
    a.set_ylabel("m per gene (outside CpG)", fontsize=10.5)
    a.set_ylim(0.8, 1.9)
    style_axis(a, grid_axis="y")

    for i in range(4):
        arrow(
            bg,
            xs[i] + cw + 4,
            y0 + ch / 2,
            xs[i + 1] - 4,
            y0 + ch / 2,
            color=INK2,
            lw=1.8,
        )

    # The patches between the steps.
    bg.text(
        30,
        478,
        "Patches the steps need because they're fit apart",
        fontsize=13.5,
        weight="bold",
        va="center",
        color=PATCH_C,
    )
    patch(
        xs[1],
        500,
        cw,
        "unweighted lines: total obs/exp 1.016;\n" "each context keeps a ±5–20% offset",
    )
    patch(
        xs[2],
        500,
        cw,
        "a separate line in AN %, fit on low-\ncoverage sites only; it reaches 1.16\njust below 90, then drops to 1",
    )
    patch(
        xs[3],
        500,
        cw,
        "too strong on gnomAD: needs α or a\n"
        "penalty; methylation features dropped\nat CpG (overlap with μ)",
    )
    patch(
        xs[4],
        500,
        cw,
        "CpG set to 1 (overlap with μ);\n" "too strong with the ridge on top",
    )

    # How it comes together.
    bg.add_patch(
        FancyBboxPatch(
            (30, 610),
            W - 60,
            340,
            boxstyle="round,pad=0,rounding_size=10",
            facecolor="white",
            edgecolor=REC_C,
            lw=1.6,
        )
    )
    bg.text(
        52,
        642,
        "Put together, per possible variant v, then summed",
        fontsize=15,
        weight="bold",
        va="center",
    )
    bg.text(
        W / 2,
        720,
        r"$E_v \;=\; f_{g(v)}(\mu_v)\;\times\; h(c_v)\;\times\; r_{w(v),\,k(v)}\;\times\; m_v$",
        fontsize=26,
        ha="center",
        va="center",
    )
    labels = [
        (W / 2 - 330, "plateau (2)", INK2),
        (W / 2 - 100, "coverage (3)", INK2),
        (W / 2 + 125, "rr (4)", RR_EDGE),
        (W / 2 + 300, "rate (5)", ROULETTE_C),
    ]
    for x, t, c in labels:
        bg.text(x, 770, t, fontsize=12.5, color=c, ha="center", va="center")
    bg.text(
        W / 2,
        830,
        r"$E_{\mathrm{gene,\,consequence}} \;=\; \sum_{v\,\in\,\mathrm{gene,\,consequence}} E_v$",
        fontsize=22,
        ha="center",
        va="center",
    )
    bg.text(
        W / 2,
        900,
        "μ is only an input to f; rr and m multiply afterwards; each factor is fit on"
        " different data, so none of them knows what the others already explain",
        fontsize=12.5,
        color=INK2,
        ha="center",
        va="center",
    )
    fig.savefig(
        os.path.join(OUT_DIR, "one_model_stack.png"),
        dpi=170,
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)


FIT_PREFIX = "data/one_model/one_model_full_v2"
FIT_LABELS = {
    "current stack, no rr": "current stack, no rr",
    "current stack + ridge rr": "current stack + ridge rr",
    "current stack x Roulette": "current stack × Roulette",
    "current stack x Roulette + ridge rr": "current stack × Roulette + ridge rr",
    "one model: groups + coverage": "groups + coverage",
    "one model: + window MQ": "+ window MQ",
    "one model: + DNM regional": "+ DNM regional",
    "one model: + DNM regional + window MQ": "+ DNM regional + window MQ",
    "one model: + gnomAD regional + window MQ": "+ gnomAD regional + window MQ",
    "one model: + window MQ + Roulette": "+ window MQ + Roulette",
    "one model: + DNM regional + window MQ + Roulette": "+ DNM regional + window MQ + Roulette",
}
FIT_SHOWN = [
    ("current stack, no rr", "current stack, no rr", NONE_C, True),
    ("current stack + ridge rr", "current stack + ridge rr", RR_EDGE, True),
    ("current stack x Roulette", "current stack × Roulette", ROULETTE_C, False),
    ("one model: groups + coverage", "one model: groups + coverage", LEGACY_C, True),
    (
        "one model: + DNM regional + window MQ + Roulette",
        "one model: + DNM regional + MQ + Roulette",
        ROULETTE_C,
        True,
    ),
]
"""Models in panels b and c: (name, label, color, filled marker)."""


def draw_fit_results(prefix: str = FIT_PREFIX) -> None:
    """Draw how the single cloglog models compare with the current stack."""
    m = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.metrics.tsv"), sep="\t")
    cov = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.coverage.tsv"), sep="\t")
    m = m.set_index("model").loc[[k for k in FIT_LABELS if k in set(m.model)]]

    fig = plt.figure(figsize=(18, 7.6))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.35, 1],
        hspace=0.55,
        wspace=0.32,
        left=0.2,
        right=0.98,
        top=0.9,
        bottom=0.1,
    )

    # (a) Per-gene noise.
    a = fig.add_subplot(gs[:, 0])
    y = np.arange(len(m))[::-1]
    stack = m.index.str.startswith("current")
    for yi, (name, r), st in zip(y, m.iterrows(), stack):
        c = NONE_C if st else LEGACY_C
        a.plot(
            [r.core_noise, r.outcome_free_noise], [yi, yi], color=c, lw=1.2, alpha=0.6
        )
        a.scatter(r.core_noise, yi, s=70, color=c, zorder=3, edgecolor=SURFACE, lw=1.5)
        a.scatter(
            r.outcome_free_noise,
            yi,
            s=70,
            facecolor=SURFACE,
            edgecolor=c,
            lw=1.8,
            zorder=3,
        )
        a.text(
            min(r.core_noise, r.outcome_free_noise) - 0.06,
            yi,
            f"{r.core_noise:.2f}",
            ha="right",
            va="center",
            fontsize=11,
            color=INK2,
        )
    ref = m.loc["current stack + ridge rr", "core_noise"]
    a.axvline(ref, color=RR_EDGE, lw=1, ls=(0, (3, 3)))
    a.axhline(y[stack.sum() - 1] - 0.5, color=AXIS, lw=1)
    a.set_yticks(y, [FIT_LABELS[k] for k in m.index], fontsize=12)
    for lab, yy in (
        ("current stack", y[0] + 0.55),
        (
            "one cloglog model, fit on the other half of the chromosomes",
            y[stack.sum()] + 0.55,
        ),
    ):
        a.text(
            0.01,
            yy,
            lab,
            transform=a.get_yaxis_transform(),
            fontsize=11.5,
            color=INK2,
            style="italic",
            va="center",
        )
    a.set_ylim(-0.6, y[0] + 0.9)
    a.set_xlim(1.2, 3.8)
    a.set_xlabel("per-gene synonymous noise (mean (obs − exp)² / exp; 1 = Poisson)")
    a.scatter([], [], s=70, color=INK2, label="core genes (unflagged, exp ≥ 20)")
    a.scatter(
        [],
        [],
        s=70,
        facecolor=SURFACE,
        edgecolor=INK2,
        lw=1.8,
        label="outcome-free set (no gene QC flag)",
    )
    a.legend(loc="lower right", fontsize=11.5, frameon=False)
    a.set_title("a  Per-gene noise, held-out chromosome halves", loc="left")
    style_axis(a, grid_axis="x")

    # (b) obs/exp by coverage.
    a = fig.add_subplot(gs[0, 1])
    bins = ["[20, 50)", "[50, 80)", "[80, 90)", "[90, 101)"]
    labels = ["20–50", "50–80", "80–90", "≥ 90"]
    for k, lab, c, filled in FIT_SHOWN:
        d = cov[cov.model == k].set_index("bin").reindex(bins)
        a.plot(
            range(len(bins)),
            d.oe,
            color=c,
            lw=2,
            marker="o",
            ms=8,
            label=lab,
            ls="-" if filled else (0, (4, 2)),
            mfc=c if filled else SURFACE,
            mew=1.8,
        )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.set_xticks(range(len(bins)), labels)
    a.set_xlabel("exome AN %")
    a.set_ylabel("synonymous obs/exp")
    a.set_title("b  By coverage", loc="left")
    a.legend(loc="upper right", fontsize=10.5, frameon=False)
    style_axis(a, grid_axis="y")

    # (c) obs/exp by context class.
    a = fig.add_subplot(gs[1, 1])
    cols = ["obs_over_exp", "oe_A/T", "oe_C non-CpG", "oe_CpG"]
    for j, (k, lab, c, filled) in enumerate(FIT_SHOWN):
        a.scatter(
            np.arange(len(cols)) + (j - 2) * 0.14,
            m.loc[k, cols],
            s=60,
            zorder=3,
            facecolor=c if filled else SURFACE,
            edgecolor=c if not filled else SURFACE,
            lw=1.8 if not filled else 1,
        )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.set_xticks(range(len(cols)), ["all", "A/T", "C non-CpG", "CpG"])
    a.set_ylabel("synonymous obs/exp")
    a.set_title("c  By context class (colors as in b)", loc="left")
    style_axis(a, grid_axis="y")

    fig.savefig(
        os.path.join(OUT_DIR, "one_model_fit.png"),
        dpi=150,
        bbox_inches="tight",
        facecolor=SURFACE,
    )
    plt.close(fig)


def draw_fit_checks(prefix: str = "data/one_model/one_model_eval") -> None:
    """Draw the report's other synonymous checks for the single cloglog models."""
    out = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.outliers.tsv"), sep="\t")
    sc = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.scores.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.deciles.tsv"), sep="\t")
    hot = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.hotspots.tsv"), sep="\t")
    order = [k for k in FIT_LABELS if k in set(sc.model)]

    fig = plt.figure(figsize=(20, 13))
    gs = fig.add_gridspec(
        2,
        6,
        hspace=0.42,
        wspace=0.9,
        left=0.15,
        right=0.98,
        top=0.93,
        bottom=0.07,
        width_ratios=[1, 1, 1, 0.55, 1, 1],
    )

    # (a) Outliers by category.
    a = fig.add_subplot(gs[0, :3])
    cats = [
        ("unflagged", "excess", "no QC flag, excess", INK),
        ("unflagged", "deficit", "no QC flag, deficit", INK2),
        ("histone genes", None, "histone genes", "#b0aa9c"),
        ("low_exome_coverage", None, "low exome coverage flag", "#d6d2c4"),
        (
            "low_exome_mapping_quality",
            "deficit",
            "low mapping quality flag, deficit",
            NONE_C,
        ),
        (
            "low_exome_mapping_quality",
            "excess",
            "low mapping quality flag, excess",
            "#c9c4b6",
        ),
    ]
    y = np.arange(len(order))[::-1]
    left = np.zeros(len(order))
    for cat, direction, label, color in cats:
        m = out.category.eq(cat) & (out.direction.eq(direction) if direction else True)
        n = out[m].groupby("model").n.sum().reindex(order).fillna(0).to_numpy()
        a.barh(
            y,
            n,
            left=left,
            color=color,
            edgecolor=SURFACE,
            lw=1,
            height=0.72,
            label=label,
            hatch="//" if direction == "excess" and cat != "unflagged" else None,
        )
        left += n
    for yi, k in zip(y, order):
        n = out[out.model.eq(k) & out.category.eq("unflagged")].n.sum()
        a.text(
            left[order.index(k)] + 3,
            yi,
            f"{int(left[order.index(k)])} ({int(n)} unflagged)",
            va="center",
            fontsize=11,
            color=INK2,
        )
    a.set_yticks(y, [FIT_LABELS[k] for k in order], fontsize=11.5)
    a.axhline(y[3] - 0.5, color=AXIS, lw=1)
    a.set_xlim(0, 330)
    a.set_xlabel("synonymous outliers (|z| > 8)")
    a.set_title("a  Outliers by release flag", loc="left")
    a.legend(loc="lower right", fontsize=10.5, frameon=False)
    style_axis(a)

    # (b) Calibration slope against the current stack.
    b = fig.add_subplot(gs[0, 4:])
    shown = [
        k
        for k in order
        if k not in ("current stack, no rr", "one model: groups + coverage")
    ]
    yb = np.arange(len(shown))[::-1]
    v = sc.set_index("model").slope_vs_current.reindex(shown)
    for yi, k in zip(yb, shown):
        c = NONE_C if k.startswith("current") else LEGACY_C
        b.plot([1, v[k]], [yi, yi], color=c, lw=1.2, alpha=0.6)
        b.scatter(v[k], yi, s=70, color=c, zorder=3, edgecolor=SURFACE, lw=1.2)
        b.text(
            v[k] + (0.03 if v[k] >= 1 else -0.03),
            yi,
            f"{v[k]:.2f}",
            va="center",
            ha="left" if v[k] >= 1 else "right",
            fontsize=11,
            color=INK2,
        )
    b.axvline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    b.set_yticks(yb, [FIT_LABELS[k] for k in shown], fontsize=11)
    b.set_xlim(0.4, 1.4)
    b.set_xlabel(
        "calibration slope against the current stack\n(below 1: too strong; above 1: too weak)"
    )
    b.set_title("b  How strong the changes are", loc="left")
    style_axis(b)

    # (c, d) Synonymous obs/exp by LOEUF and s_het decile.
    x = np.arange(10)
    for col, binning, title, xlabel in (
        (
            0,
            "decile",
            "c  By LOEUF decile",
            "LOEUF decile (v4.1.2; 1 = most constrained)",
        ),
        (2, "shet_decile", "d  By s_het decile", "s_het decile (1 = highest s_het)"),
    ):
        a = fig.add_subplot(gs[1, col : col + 2] if col == 0 else gs[1, 2:4])
        for k, lab, c, filled in FIT_SHOWN:
            d = dec[dec.model.eq(k) & dec.binning.eq(binning)].sort_values("decile")
            a.plot(
                x,
                d.oe,
                color=c,
                lw=2,
                marker="o",
                ms=6,
                label=lab,
                ls="-" if filled else (0, (4, 2)),
                mfc=c if filled else SURFACE,
                mew=1.6,
            )
        a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
        a.set_xticks(x, [str(i + 1) for i in x])
        a.set_ylim(0.96, 1.14)
        a.set_xlabel(xlabel)
        a.set_ylabel("synonymous obs/exp")
        a.set_title(title, loc="left")
        style_axis(a, grid_axis="y")
        if col == 0:
            a.legend(loc="upper left", fontsize=10, frameon=False)

    # (e) Maternal hotspot regions: the control genes there, as in the report.
    a = fig.add_subplot(gs[1, 4:])
    regions = list(dict.fromkeys(hot.region))
    for j, (k, lab, c, filled) in enumerate(FIT_SHOWN):
        h = (
            hot[hot.model.eq(k) & hot.set.eq("control genes")]
            .set_index("region")
            .reindex(regions)
        )
        a.scatter(
            np.arange(len(regions)) + (j - 2) * 0.13,
            h.oe,
            s=55,
            zorder=3,
            facecolor=c if filled else SURFACE,
            edgecolor=c if not filled else SURFACE,
            lw=1.6 if not filled else 1,
        )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.set_xticks(range(len(regions)), [r.replace(", first 20 Mb", "") for r in regions])
    a.set_ylabel("obs/exp of the control genes there\n(relative to the core set)")
    a.set_title("e  Maternal hotspot regions (colors as in c)", loc="left")
    style_axis(a, grid_axis="y")

    fig.savefig(
        os.path.join(OUT_DIR, "one_model_checks.png"),
        dpi=150,
        bbox_inches="tight",
        facecolor=SURFACE,
    )
    plt.close(fig)


def draw_release_checks(prefix: str = "data/one_model/one_model") -> None:
    """Draw the single cloglog models' pLoF and missense results, and NDD/HI PR."""
    sc = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}_release.scores.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}_release.deciles.tsv"), sep="\t")
    sc = sc.set_index("model")
    others = FIT_SHOWN[1:]

    fig, axes = plt.subplots(
        1, 4, figsize=(24, 6.2), gridspec_kw={"width_ratios": [1, 1, 0.9, 1.1]}
    )
    fig.subplots_adjust(left=0.04, right=0.99, top=0.8, bottom=0.17, wspace=0.62)
    x = np.arange(10)
    for ax, csq, what in ((axes[0], "lof", "a  pLoF"), (axes[1], "mis", "b  Missense")):
        base = (
            dec[dec.model.eq(FIT_SHOWN[0][0]) & dec.csq.eq(csq)].set_index("decile").oe
        )
        for k, lab, c, filled in others:
            d = dec[dec.model.eq(k) & dec.csq.eq(csq)].set_index("decile").oe
            ax.plot(
                x,
                100 * (d / base - 1),
                color=c,
                lw=2,
                marker="o",
                ms=6,
                label=lab,
                ls="-" if filled else (0, (4, 2)),
                mfc=c if filled else SURFACE,
                mew=1.6,
            )
        ax.axhline(0, color=NONE_C, lw=2.6)
        ax.set_xticks(x, [str(i + 1) for i in x])
        ax.set_xlabel("LOEUF decile (v4.1.2; 1 = most constrained)")
        ax.set_ylabel("obs/exp, % change from the current stack")
        ax.set_title(f"{what} obs/exp, against the current stack", loc="left")
        style_axis(ax, grid_axis="y")
    lo = min(a.get_ylim()[0] for a in axes[:2])
    hi = max(a.get_ylim()[1] for a in axes[:2])
    for a in axes[:2]:
        a.set_ylim(lo, hi)
    axes[0].legend(loc="best", fontsize=10.5, frameon=False)

    # (c) pLoF and missense relative to synonymous.
    ax = axes[2]
    ref = sc.loc[FIT_SHOWN[0][0]]
    for gi, (col, glabel) in enumerate(
        (
            ("lof_over_syn", "pLoF ÷ synonymous"),
            ("mis_over_syn", "missense ÷ synonymous"),
        )
    ):
        for j, (k, lab, c, filled) in enumerate(others):
            v = 100 * (sc.loc[k, col] / ref[col] - 1)
            y = gi * 5 + j
            ax.barh(
                y, v, height=0.8, color=c if filled else SURFACE, edgecolor=c, lw=1.6
            )
            ax.text(
                v + (0.15 if v >= 0 else -0.15),
                y,
                f"{v:+.1f}%",
                va="center",
                ha="left" if v >= 0 else "right",
                fontsize=11.5,
            )
        ax.text(
            0.02,
            gi * 5 - 0.95,
            glabel,
            fontsize=12,
            color=INK2,
            va="center",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color=NONE_C, lw=2.6)
    ax.set_yticks([])
    ax.set_ylim(9.0, -1.8)
    lim = max(abs(np.array(ax.get_xlim()))) + 1.5
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("% change from the current stack")
    ax.set_title("c  Levels relative to synonymous", loc="left")
    style_axis(ax)

    # (d) Genes crossing LOEUF 0.6, and changing decile.
    ax = axes[3]
    y = np.arange(len(others))[::-1]
    for yi, (k, lab, c, filled) in zip(y, others):
        out, into = int(sc.loc[k, "out_of_loeuf_lt_0.6"]), int(
            sc.loc[k, "into_loeuf_lt_0.6"]
        )
        kw = dict(height=0.6, color=c if filled else SURFACE, edgecolor=c, lw=1.6)
        ax.barh(yi, -out, **kw)
        ax.barh(yi, into, **kw)
        ax.text(
            -out - 10,
            yi,
            f"−{out}" if out else "0",
            ha="right",
            va="center",
            fontsize=12,
        )
        ax.text(into + 10, yi, f"+{into}" if into else "0", va="center", fontsize=12)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(
        y,
        [
            lab.replace("one model: + ", "one model:\n+ ")
            + f"\n({100 * sc.loc[k, 'decile_changed']:.1f}% change LOEUF decile)"
            for k, lab, c, filled in others
        ],
        fontsize=10.5,
    )
    m = max(
        sc.loc[
            [v[0] for v in others], ["out_of_loeuf_lt_0.6", "into_loeuf_lt_0.6"]
        ].max()
    )
    ax.set_xlim(-m * 1.35 - 40, m * 1.35 + 40)
    ax.set_xlabel("← no longer below 0.6      newly below 0.6 →")
    ax.set_title("d  Genes crossing LOEUF 0.6", loc="left")
    style_axis(ax)

    fig.suptitle(
        "The single cloglog models on pLoF and missense (fit on synonymous variants,"
        " applied to the other chromosome half)",
        x=0.04,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    fig.text(
        0.04,
        0.87,
        f"{int(sc.n.iloc[0]):,} autosomal Ensembl MANE Select transcripts with a pLoF"
        " expected count; LOEUF from the 90% gamma interval, as in the pipeline."
        f" ClinGen HI AUC of LOEUF: {sc.hi_auc.min():.3f}–{sc.hi_auc.max():.3f}.",
        fontsize=13,
        color=INK2,
    )
    fig.savefig(
        os.path.join(OUT_DIR, "one_model_release.png"),
        dpi=150,
        bbox_inches="tight",
        facecolor=SURFACE,
    )
    plt.close(fig)

    draw_disease_gene_pr(
        f"{prefix}_pr",
        versions=[
            (KEYS[k], lab, c, "-" if filled else (0, (4, 2)))
            for k, lab, c, filled in FIT_SHOWN
        ],
        shown=[KEYS[k] for k, *_ in FIT_SHOWN],
        emph=KEYS[FIT_SHOWN[-1][0]],
        baseline="the current stack",
        out="one_model_pr.png",
        suptitle="Classifying NDD and haploinsufficient genes by LOEUF, for the single"
        " cloglog models",
    )


ISSUE_FAMILY_C = {
    "current": NONE_C,
    "rr": RR_EDGE,
    "roulette": ROULETTE_C,
    "plateau": "#1baf7a",
    "one model": LEGACY_C,
}
"""Colors of the method families in the issue figures."""
ISSUE_LINES = [
    ("current stack, no rr", NONE_C, "-"),
    ("current stack + ridge rr", RR_EDGE, "-"),
    ("current stack x Roulette", ROULETTE_C, (0, (4, 2))),
    ("plateau: saturating", "#1baf7a", "-"),
    ("one model: groups + coverage", LEGACY_C, "-"),
    ("one model: + DNM regional + window MQ + Roulette", "#0b2f6b", "-"),
]
"""Methods drawn in the issue figures' line plots: (name, color, line style)."""


def issue_family(m: str) -> str:
    """Get a method's family for its color."""
    if m == "current stack, no rr":
        return "current"
    if m.startswith("plateau"):
        return "plateau"
    if m.startswith("one model"):
        return "one model"
    return "roulette" if "Roulette" in m else "rr"


def draw_issue_figures(prefix: str = "data/one_model/issues") -> None:
    """Draw one figure per problem of the current expected counts, with each fix."""
    sm = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.summary.tsv"), sep="\t")
    ctx = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.contexts.tsv"), sep="\t")
    cov = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.coverage.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.deciles.tsv"), sep="\t")
    hot = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.hotspots.tsv"), sep="\t")
    sm = sm.set_index("method")
    order = list(sm.index)
    labels = sm.label.to_dict()
    y = np.arange(len(order))[::-1]
    colors = [ISSUE_FAMILY_C[issue_family(m)] for m in order]
    line_labels = {m: labels[m] for m, *_ in ISSUE_LINES}

    def method_axis(ax):
        ax.set_yticks(y, [labels[m] for m in order], fontsize=11)
        ax.get_yticklabels()[0].set_fontweight("bold")
        style_axis(ax)

    def save(fig, name):
        fig.savefig(
            os.path.join(OUT_DIR, f"one_model_issue_{name}.png"),
            dpi=150,
            bbox_inches="tight",
            facecolor=SURFACE,
        )
        plt.close(fig)

    def dots(ax, col, fmt, ref=None):
        v = sm[col].reindex(order).to_numpy(float)
        if ref is not None:
            ax.axvline(ref, color=INK2, lw=1, ls=(0, (3, 3)))
            ax.hlines(y, ref, v, color=colors, lw=1.2, alpha=0.6)
        ax.scatter(v, y, s=70, color=colors, zorder=3, edgecolor=SURFACE, lw=1)
        for yi, vi in zip(y, v):
            ax.text(vi, yi + 0.3, fmt.format(vi), ha="center", fontsize=10, color=INK2)

    # 1. Synonymous level.
    fig, ax = plt.subplots(figsize=(10, 6))
    dots(ax, "syn_oe", "{:.3f}", ref=1)
    method_axis(ax)
    ax.set_xlabel("synonymous obs/exp, all exported MANE Select synonymous sites")
    ax.set_title("Synonymous level", loc="left")
    save(fig, "level")

    # 2. Offsets between contexts.
    fig, ax = plt.subplots(figsize=(11, 6.5))
    rng = np.random.default_rng(0)
    for yi, m, c in zip(y, order, colors):
        v = ctx[ctx.method == m].oe.to_numpy()
        ax.scatter(
            v, yi + rng.uniform(-0.18, 0.18, len(v)), s=18, color=c, alpha=0.7, lw=0
        )
        ax.text(
            1.31,
            yi,
            f"SD {100 * sm.loc[m, 'context_sd']:.1f}%",
            va="center",
            fontsize=10.5,
        )
    ax.axvline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.set_xlim(0.72, 1.38)
    method_axis(ax)
    ax.set_xlabel("synonymous obs/exp of each of the 32 contexts (total scaled to 1)")
    ax.set_title("Offsets between contexts", loc="left")
    save(fig, "contexts")

    # 3. Low-coverage calibration.
    fig, ax = plt.subplots(figsize=(12, 5.5))
    for m, c, ls in ISSUE_LINES:
        if m in ("current stack + ridge rr", "current stack x Roulette"):
            continue
        d = cov[cov.method == m].sort_values("coverage")
        ax.plot(
            d.coverage,
            d.oe,
            color=c,
            ls=ls,
            lw=2 if m == ISSUE_LINES[0][0] else 1.6,
            marker="o",
            ms=3.5,
            label=line_labels[m],
        )
    ax.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.axvline(90, color=AXIS, lw=1)
    ax.set_ylim(0.5, 2.2)
    ax.set_xlabel("exome AN %")
    ax.set_ylabel("synonymous obs/exp")
    ax.legend(loc="upper right", fontsize=10.5, frameon=False)
    ax.set_title(
        "Low-coverage calibration (rr and Roulette leave it as it is)", loc="left"
    )
    style_axis(ax, grid_axis="y")
    save(fig, "coverage")

    # 4. Per-gene noise.
    fig, ax = plt.subplots(figsize=(10, 6))
    for yi, m, c in zip(y, order, colors):
        a, b = sm.loc[m, "core_noise"], sm.loc[m, "outcome_free_noise"]
        ax.plot([a, b], [yi, yi], color=c, lw=1.2, alpha=0.6)
        ax.scatter(a, yi, s=70, color=c, zorder=3, edgecolor=SURFACE, lw=1)
        ax.scatter(b, yi, s=70, facecolor=SURFACE, edgecolor=c, lw=1.6, zorder=3)
        ax.text(
            min(a, b) - 0.05,
            yi,
            f"{a:.2f}",
            ha="right",
            va="center",
            fontsize=10,
            color=INK2,
        )
    ax.axvline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.text(1.02, y[-1] - 0.6, "Poisson", fontsize=10, color=INK2)
    ax.set_xlim(0.8, 4)
    method_axis(ax)
    ax.set_xlabel(
        "per-gene synonymous noise, mean (obs − exp)² / exp"
        "\n(filled: core genes; hollow: genes with no QC flag)"
    )
    ax.set_title("Regional variation: per-gene noise", loc="left")
    save(fig, "noise")

    # 5. Maternal hotspot regions.
    fig, ax = plt.subplots(figsize=(12, 5))
    regions = list(dict.fromkeys(hot.region))
    shown = [m for m, *_ in ISSUE_LINES] + ["current stack x Roulette + ridge rr"]
    for j, m in enumerate(shown):
        h = hot[hot.method == m].set_index("region").reindex(regions)
        c = dict((k, cc) for k, cc, _ in ISSUE_LINES).get(m, "#c4302b")
        ax.scatter(
            np.arange(len(regions)) + (j - 3) * 0.11,
            h.oe,
            s=45,
            color=c,
            zorder=3,
            label=labels[m],
            edgecolor=SURFACE,
            lw=0.8,
        )
    ax.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.set_xticks(range(len(regions)), regions)
    ax.set_ylabel("synonymous obs/exp of control genes\n(relative to the core set)")
    ax.legend(loc="upper right", fontsize=9.5, frameon=False, ncol=2)
    ax.set_title("Maternal hotspot regions", loc="left")
    style_axis(ax, grid_axis="y")
    save(fig, "hotspots")

    # 6. The rise in the least constrained LOEUF decile.
    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(10)
    for m, c, ls in ISSUE_LINES:
        d = dec[dec.method == m].sort_values("decile")
        ax.plot(x, d.oe, color=c, ls=ls, lw=2, marker="o", ms=5, label=line_labels[m])
    ax.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.set_xticks(x, [str(i + 1) for i in x])
    ax.set_xlabel("LOEUF decile (v4.1.2; 1 = most constrained)")
    ax.set_ylabel("synonymous obs/exp")
    ax.legend(loc="upper left", fontsize=10, frameon=False)
    ax.set_title("Rise in the least constrained LOEUF decile", loc="left")
    style_axis(ax, grid_axis="y")
    save(fig, "loeuf_rise")

    # 7-9. Outliers with no QC flag, low-mapping-quality genes, histone genes.
    for name, cols, title, xlabel in (
        (
            "outliers",
            ["outliers_no_flag_excess", "outliers_no_flag_deficit"],
            "Synonymous outliers with no QC flag (|z| > 8)",
            "genes (filled: excess; hollow: deficit)",
        ),
        (
            "low_mq",
            ["low_mq_deficit", "low_mq_excess"],
            "Genes flagged for low exome mapping quality: outliers",
            "genes (filled: deficit, looking depleted; hollow: excess, overcorrected)",
        ),
        ("histones", ["histone_outliers"], "Histone genes", "outliers (|z| > 8)"),
    ):
        fig, ax = plt.subplots(figsize=(10, 6))
        left = np.zeros(len(order))
        for j, col in enumerate(cols):
            v = sm[col].reindex(order).to_numpy(float)
            ax.barh(
                y,
                v,
                left=left,
                height=0.65,
                color=colors if j == 0 else SURFACE,
                edgecolor=colors,
                lw=1.4,
            )
            left += v
        for yi, m, lv in zip(y, order, left):
            extra = (
                f"  (pooled obs/exp {sm.loc[m, 'histone_oe']:.2f})"
                if name == "histones"
                else ""
            )
            ax.text(
                lv + 1, yi, f"{int(lv)}{extra}", va="center", fontsize=10.5, color=INK2
            )
        ax.set_xlim(0, max(left) * 1.4 + 5)
        method_axis(ax)
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left")
        save(fig, name)

    # 10. Side effects: pLoF and missense levels, NDD ranking.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"wspace": 0.08})
    for dy, col, hollow in ((0.17, "lof_vs_syn", False), (-0.17, "mis_vs_syn", True)):
        v = sm[col].reindex(order).to_numpy(float)
        a1.barh(
            y + dy,
            np.nan_to_num(v),
            height=0.32,
            color=SURFACE if hollow else colors,
            edgecolor=colors,
            lw=1.4,
        )
        for yi, vi in zip(y, v):
            a1.text(
                (vi if np.isfinite(vi) else 0)
                + (0.2 if not np.isfinite(vi) or vi >= 0 else -0.2),
                yi + dy,
                f"{vi:+.1f}%" if np.isfinite(vi) else "not run",
                va="center",
                ha="left" if not np.isfinite(vi) or vi >= 0 else "right",
                fontsize=9.5,
                color=INK2,
            )
    a1.axvline(0, color=INK2, lw=1)
    a1.set_xlim(-6, 10)
    method_axis(a1)
    a1.set_xlabel("% change relative to synonymous (filled: pLoF; hollow: missense)")
    a1.set_title("Side effect: pLoF and missense levels", loc="left")
    v, lo, hi = (
        sm[k].reindex(order).to_numpy(float) for k in ("ndd_change", "ndd_lo", "ndd_hi")
    )
    a2.hlines(y, lo, hi, color=colors, lw=1.6)
    a2.scatter(v, y, s=60, color=colors, zorder=3, edgecolor=SURFACE, lw=1)
    a2.axvline(0, color=INK2, lw=1)
    a2.set_yticks(y, [""] * len(y))
    a2.set_xlabel("NDD AUPRC, change from the current (95% interval)")
    a2.set_title("Outcome: ranking NDD genes", loc="left")
    style_axis(a2)
    save(fig, "side_effects")


def draw_coverage_check() -> None:
    """Draw the stored coverage model against the share seen at each AN %."""
    g = pd.read_csv(
        os.path.join(OUT_DIR, "data/one_model/coverage_check.by_an.tsv"), sep="\t"
    )
    fig, (a, b) = plt.subplots(
        2,
        1,
        figsize=(12, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12},
    )
    size = 8 + 160 * np.sqrt(g.possible / g.possible.max())
    a.scatter(
        g.coverage,
        g.y_pipeline,
        s=size,
        facecolor=SURFACE,
        edgecolor=NONE_C,
        lw=1.4,
        zorder=3,
        label="what the coverage model is fit to: obs / (scale × Σ possible·μ)",
    )
    a.scatter(
        g.coverage,
        g.y_plateau,
        s=size,
        color=LEGACY_C,
        alpha=0.8,
        lw=0,
        zorder=4,
        label="relative to the plateau it multiplies: obs / Σ possible·f(μ)",
    )
    x = np.linspace(20, 89.99, 100)
    a.plot(
        x,
        COVERAGE_MODEL[0] + COVERAGE_MODEL[1] * x,
        color=INK,
        lw=2.2,
        zorder=5,
        label="stored coverage model h (a straight line, one point per AN % value)",
    )
    a.plot([90, 100], [1, 1], color=INK, lw=2.2, zorder=5)
    a.plot(
        [89.99, 90],
        [COVERAGE_MODEL[0] + COVERAGE_MODEL[1] * 89.99, 1],
        color=INK,
        lw=1,
        ls=(0, (2, 2)),
    )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.axvline(90, color=AXIS, lw=1)
    a.set_ylabel("share seen, relative to AN ≥ 90")
    a.set_ylim(0, 1.3)
    a.legend(loc="lower right", fontsize=11, frameon=False)
    a.set_title(
        "Synonymous variants seen by exome AN %, and the coverage model", loc="left"
    )
    style_axis(a, grid_axis="both")

    b.bar(g.coverage, g.possible, width=0.8, color=NONE_C, lw=0)
    b.set_yscale("log")
    b.axvline(90, color=AXIS, lw=1)
    b.set_ylabel("possible variants")
    b.set_xlabel("exome AN % (sites below 20 are dropped by the pipeline)")
    style_axis(b, grid_axis="y")
    fig.savefig(
        os.path.join(OUT_DIR, "one_model_coverage_check.png"),
        dpi=150,
        bbox_inches="tight",
        facecolor=SURFACE,
    )
    plt.close(fig)


def draw_coverage_oe_current(prefix: str = FIT_PREFIX) -> None:
    """Draw synonymous obs/exp by exome AN % under the current models."""
    cov = pd.read_csv(os.path.join(OUT_DIR, f"{prefix}.coverage.tsv"), sep="\t")
    cov = cov[cov.model == "current stack, no rr"]
    g = pd.read_csv(
        os.path.join(OUT_DIR, "data/one_model/coverage_check.by_an.tsv"), sep="\t"
    )
    fig, a = plt.subplots(figsize=(12, 5.6))
    a.scatter(
        g.coverage,
        g.oe_stack,
        s=8 + 160 * np.sqrt(g.possible / g.possible.max()),
        color=NONE_C,
        alpha=0.85,
        lw=0,
        zorder=3,
        label="each AN % value (sized by possible variants)",
    )
    for lo, hi, b in (
        (20, 50, "[20, 50)"),
        (50, 80, "[50, 80)"),
        (80, 90, "[80, 90)"),
        (90, 101, "[90, 101)"),
    ):
        oe = cov.set_index("bin").loc[b, "oe"]
        a.plot(
            [lo, hi - 1],
            [oe, oe],
            color=INK,
            lw=3,
            solid_capstyle="butt",
            zorder=4,
            label="each bin, as in panel b of the fit figure" if lo == 20 else None,
        )
        a.text(
            (lo + hi - 1) / 2,
            oe + (0.1 if lo == 90 else 0.05),
            f"{oe:.2f}",
            ha="center",
            fontsize=12,
            weight="bold",
        )
    a.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)))
    a.axvline(90, color=AXIS, lw=1)
    a.set_xlabel("exome AN %")
    a.set_ylabel("synonymous obs/exp")
    a.set_ylim(0.5, 2.25)
    a.set_title(
        "Synonymous obs/exp by exome AN %, current plateau and coverage models (no rr)",
        loc="left",
    )
    a.legend(loc="upper right", fontsize=11, frameon=False)
    style_axis(a, grid_axis="y")
    fig.savefig(
        os.path.join(OUT_DIR, "one_model_coverage_oe.png"),
        dpi=150,
        bbox_inches="tight",
        facecolor=SURFACE,
    )
    plt.close(fig)


if __name__ == "__main__":
    draw_stack()
    draw_fit_results()
    draw_fit_checks()
    draw_release_checks()
    draw_issue_figures()
    draw_coverage_check()
    draw_coverage_oe_current()
