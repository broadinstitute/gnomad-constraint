"""Check the current coverage model against the exported synonymous cells.

For the notes ``notes_one_model_for_expected.md``, linked from the rr report.

1. The cells' current expected (``e_old``, possible * stored plateau line * coverage
   correction) summed per MANE Select transcript, against the release's ``syn.exp``.
2. Per exome AN % value, the share of possible variants seen relative to high
   coverage, two ways:

   - ``y_pipeline``: what ``build_models`` fits the coverage model to,
     obs / (scale * sum(possible * mu)), with scale = obs / sum(possible * mu) at
     AN >= 90;
   - ``y_plateau``: relative to the plateau prediction the model is applied to,
     obs / sum(possible * f(mu)), over the same ratio at AN >= 90.

   And the coverage model refit the pipeline's way (an unweighted least-squares line
   in AN % over the low-coverage values, one point per value), and weighted.

Writes ``data/one_model/coverage_check.by_an.tsv``.
"""

import numpy as np
import pandas as pd
from make_followup_figures import COVERAGE_MODEL

CELLS = "data/one_model/syn_model_cells.tsv.bgz"
RELEASE = "../gnomad.v4.1.2.constraint_metrics.tsv.bgz"
OUT = "data/one_model/coverage_check.by_an.tsv"


def main() -> None:
    """Run both checks."""
    c = pd.read_csv(
        CELLS,
        sep="\t",
        compression="gzip",
        usecols=[
            "transcript",
            "cpg",
            "coverage",
            "possible",
            "observed",
            "mu",
            "e_old",
        ],
    )
    rel = pd.read_csv(
        RELEASE,
        sep="\t",
        compression="gzip",
        usecols=["transcript", "mane_select", "syn.possible", "syn.obs", "syn.exp"],
    )
    rel = rel[rel.mane_select.astype(str).str.lower().eq("true")]
    rel = rel.drop_duplicates("transcript").set_index("transcript")
    t = (
        c.groupby("transcript")[["possible", "observed", "e_old"]]
        .sum()
        .join(rel, how="inner")
    )
    low_tx = set(c[c.coverage < 90].transcript)
    for a, b in (
        ("possible", "syn.possible"),
        ("observed", "syn.obs"),
        ("e_old", "syn.exp"),
    ):
        r = t[a] / t[b]
        rl = r[r.index.isin(low_tx)]
        print(
            f"{a} vs release {b}: {len(r)} transcripts, within 0.1%: "
            f"{(abs(r - 1) < 0.001).mean():.4f} (with sites below AN 90: "
            f"{(abs(rl - 1) < 0.001).mean():.4f} of {len(rl)})"
        )

    fits = pd.read_csv("plateau_models.fits.tsv", sep="\t")
    st = fits[fits.fit == "stored"].set_index("line")
    a = np.where(
        c.cpg,
        st.loc["CpG transition", "intercept"],
        st.loc["everything else", "intercept"],
    )
    b = np.where(
        c.cpg, st.loc["CpG transition", "slope"], st.loc["everything else", "slope"]
    )
    c["e_plat"] = c.possible * (a + b * c.mu)
    c["pmu"] = c.possible * c.mu
    hi = c.coverage >= 90
    scale = c.observed[hi].sum() / c.pmu[hi].sum()
    oe_hi = c.observed[hi].sum() / c.e_plat[hi].sum()
    g = c.groupby("coverage")[["observed", "pmu", "e_plat", "e_old", "possible"]].sum()
    g["y_pipeline"] = g.observed / (scale * g.pmu)
    g["y_plateau"] = g.observed / g.e_plat / oe_hi
    g["h"] = np.where(
        g.index >= 90, 1.0, COVERAGE_MODEL[0] + COVERAGE_MODEL[1] * g.index
    )
    g["oe_stack"] = g.observed / g.e_old
    low = g[g.index < 90]
    print("stored coverage model:", COVERAGE_MODEL)
    print("refit, unweighted:", np.polyfit(low.index, low.y_pipeline, 1)[::-1])
    print(
        "refit, weighted by possible * mu:",
        np.polyfit(low.index, low.y_pipeline, 1, w=np.sqrt(low.pmu))[::-1],
    )
    g.to_csv(OUT, sep="\t")


if __name__ == "__main__":
    main()
