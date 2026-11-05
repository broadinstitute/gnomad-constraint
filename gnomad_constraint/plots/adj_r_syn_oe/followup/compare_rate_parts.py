"""Compare the home-built per-variant rates with Roulette's, part by part.

Inputs: the exported cells with the rate versions' cell factors
(``cluster/export_rate_cells.py``), and ``roulette_combo.v2.cells.tsv.gz``, which splits
Roulette's rate per cell (window and context) into a sequence part (what the wider
sequence adds to the trinucleotide rate) and a regional part. CpG cells are left out
(both versions give them 1).

For the total, regional and sequence factors, per context class and over all cells:
the spread of each version's log factor, and how closely the two agree (correlation,
and the slope of home-built on Roulette, both weighted by expected count). The same
at the gene level (each transcript's expected-weighted factor), which is what changes
gene-level noise. Also each version's mean gene factor in the maternal hotspot
regions, and the per-gene factors (``<prefix>.genes.tsv.gz``).
"""

import argparse
import logging

import numpy as np
import pandas as pd
from outlier_drivers import MATERNAL_HOTSPOTS

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

PAIRS = {
    "total": ("rr_homebuilt_only", "roulette_multiplier"),
    "regional": ("rr_homebuilt_regional", "roulette_regional"),
    "sequence": ("rr_homebuilt_sequence", "roulette_sequence"),
}


def class_of(context: pd.Series) -> np.ndarray:
    """Get the context class of each strand-collapsed trinucleotide context."""
    return np.where(
        context.str[1] == "A",
        "A/T",
        np.where(context.str[2] == "G", "CpG", "C non-CpG"),
    )


def agree(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> dict:
    """Get weighted spreads of log x and log y, their correlation and the slope of x on y."""
    lx, ly = np.log(x), np.log(y)
    mx, my = np.average(lx, weights=w), np.average(ly, weights=w)
    vx = np.average((lx - mx) ** 2, weights=w)
    vy = np.average((ly - my) ** 2, weights=w)
    cxy = np.average((lx - mx) * (ly - my), weights=w)
    return {
        "sd_homebuilt": float(np.sqrt(vx)),
        "sd_roulette": float(np.sqrt(vy)),
        "corr": float(cxy / np.sqrt(vx * vy)),
        "slope_homebuilt_on_roulette": float(cxy / vy),
    }


def main(args: argparse.Namespace) -> None:
    """Compare the two sets of rates."""
    cells = pd.read_csv(
        args.cells,
        sep="\t",
        na_values=["NA"],
        compression="gzip" if args.cells.endswith("gz") else None,
        usecols=[
            "transcript",
            "win",
            "context",
            "exp",
            *[p[0] for p in PAIRS.values()],
        ],
    )
    rou = pd.read_csv(args.roulette_cells, sep="\t", compression="gzip").rename(
        columns={
            "multiplier": "roulette_multiplier",
            "regional_part": "roulette_regional",
            "sequence_part": "roulette_sequence",
        }
    )
    d = cells.merge(
        rou[["win", "context", *[p[1] for p in PAIRS.values()]]],
        on=["win", "context"],
    )
    d["cls"] = class_of(d.context)
    d = d[d.cls != "CpG"]
    cols = [c for p in PAIRS.values() for c in p]
    d = d.dropna(subset=cols)
    d = d[(d[cols] > 0).all(axis=1)]
    logger.info("%d non-CpG cells with both versions' factors.", len(d))

    rows = []
    for level in ("cell", "gene"):
        if level == "cell":
            g = d
        else:
            g = (
                d.assign(**{c: d[c] * d.exp for c in cols})
                .groupby(["transcript", "cls"])[[*cols, "exp"]]
                .sum()
            )
            g = g.assign(**{c: g[c] / g.exp for c in cols}).reset_index()
            g = g[g.exp >= 5]
        for cls in ("all", "A/T", "C non-CpG"):
            s = g if cls == "all" else g[g.cls == cls]
            for part, (hb, ro) in PAIRS.items():
                rows.append(
                    {
                        "level": level,
                        "class": cls,
                        "part": part,
                        "n": len(s),
                        **agree(s[hb].to_numpy(), s[ro].to_numpy(), s.exp.to_numpy()),
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}.agreement.tsv", sep="\t", index=False)

    # Mean gene factor in the maternal hotspot regions (all non-CpG sites of a gene).
    tx = pd.read_csv(args.genes, sep="\t")[["transcript", "chromosome", "mid"]]
    gene = (
        d.assign(**{c: d[c] * d.exp for c in cols})
        .groupby("transcript")[[*cols, "exp"]]
        .sum()
    )
    gene = gene.assign(**{c: gene[c] / gene.exp for c in cols}).reset_index()
    gene = gene.merge(tx, on="transcript")
    gene.to_csv(f"{args.out_prefix}.genes.tsv.gz", sep="\t", index=False)
    hot = []
    for name, (c, s, e) in MATERNAL_HOTSPOTS.items():
        r = gene[gene.chromosome.eq(c) & gene.mid.between(s, e)]
        hot.append(
            {
                "region": name,
                "n": len(r),
                **{
                    c: float(np.average(r[c], weights=r.exp))
                    for c in (
                        "rr_homebuilt_only",
                        "roulette_multiplier",
                        "rr_homebuilt_regional",
                        "roulette_regional",
                    )
                },
            }
        )
    hot = pd.DataFrame(hot)
    hot.to_csv(f"{args.out_prefix}.hotspots.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(out.round(3).to_string(index=False))
    print(hot.round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells", required=True, help="Rate cells TSV (bgzipped is fine)."
    )
    parser.add_argument("--roulette-cells", default="roulette_combo.v2.cells.tsv.gz")
    parser.add_argument("--genes", default="outlier_drivers.all_genes.tsv")
    parser.add_argument("--out-prefix", default="compare_rate_parts")
    main(parser.parse_args())
