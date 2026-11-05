"""Add the per-variant-rate versions (Roulette, home-built) to Ruchit's exported cells.

The constraint pipeline multiplies each site's expected count by its per-variant rate
multiplier and by the rr fit on top of the rates (``annotate_regional_rr_expected``).
For each synonymous Ensembl MANE Select site of the per-variant expected Table, this
takes that factor, and sums it per cell (transcript, 1-kb window, context) as
``sum(exp * factor) / sum(exp)``: the cell-level equivalent of an rr column. The
versions then score like every other rr in ``full_run_scores.py`` (``--rr-names``).

Versions: the rates alone (``<rates>_only``) and the rates with the ridge rr fit on
top (``<rates>_ridge``, held out with the out-of-fold rr). For the home-built rates,
also their regional and sequence parts alone, each normalized per context over the
synonymous sites (and 1 at CpG, like the multiplier), to see which part falls short
of Roulette. Writes the base cells with ``rr_<version>`` and ``rr_<version>_oof``
columns (versions not fit on DNMs get the same held-out column).
"""

import argparse

import hail as hl

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import VARIANT_RATE_MULTIPLIER_FIELD
from gnomad_constraint.utils.regional_rate import get_syn_mane_select_ht

RD = "gnomad_v4.1.1/constraint/regional_rate"
RATES = {
    "roulette": (
        "gs://gnomad-tmp/jg/rr_roulette/gnomad_v4.1.1/constraint/regional_rate/"
        "gnomad.v4.1.1.roulette_variant_rates.ht",
        "gs://gnomad-tmp/jg/rr_roulette",
    ),
    "homebuilt": (
        "gs://gnomad-tmp/jg/rr_homebuilt/rates/homebuilt_variant_rates.ht",
        "gs://gnomad-tmp/jg/rr_homebuilt",
    ),
}
"""Rates name -> (per-variant rate Table, root of the rr fit on top of it)."""
CELLS = (
    "gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export/"
    "syn_cells.tsv.bgz"
)


def main(args: argparse.Namespace) -> None:
    """Compute the cell factors and export the annotated cells."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    ht = constraint_res.get_per_variant_expected_dataset(version="4.1.1").ht()
    ht = get_syn_mane_select_ht(ht)
    start = (ht.locus.position - 1) // 1000 * 1000
    ht = ht.annotate(
        win=ht.locus.contig + "-" + hl.str(start) + "-" + hl.str(start + 1000)
    )
    factors = {}
    cpg = (ht.context[1] == "C") & (ht.context[2] == "G")
    parts = {}
    for name, (rate_path, root) in RATES.items():
        rate = hl.read_table(rate_path)
        row = rate[ht.locus, ht.alleles]
        m = hl.or_else(row[VARIANT_RATE_MULTIPLIER_FIELD], 1.0)
        factors[f"{name}_only"] = (m, m)
        for part in ("regional_part", "sequence_part"):
            if part in rate.row:
                parts[f"{name}_{part.split('_')[0]}"] = hl.if_else(
                    cpg, 1.0, hl.or_else(row[part], 1.0)
                )
        rr = {}
        for oof in (False, True):
            path = (
                f"{root}/{RD}/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_ridge.calibrated"
                f"{'.out_of_fold' if oof else ''}.ht"
            )
            rr_ht = hl.read_table(path).select("adj_r")
            rr[oof] = hl.or_else(rr_ht[ht.locus].adj_r.get(ht.context), 1.0)
        factors[f"{name}_ridge"] = (m * rr[False], m * rr[True])
    ht = ht.annotate(
        **{f"f_{v}": f for v, (f, _) in factors.items()},
        **{f"f_{v}_oof": f for v, (_, f) in factors.items()},
        **{f"raw_{v}": p for v, p in parts.items()},
    )
    ht = ht.checkpoint(hl.utils.new_temp_file("rate_cells.sites", "ht"))
    # Normalize each part per context over the synonymous sites.
    norms = ht.aggregate(
        hl.agg.group_by(
            ht.context,
            hl.struct(
                **{v: hl.agg.sum(ht.exp * ht[f"raw_{v}"]) for v in parts},
                _exp=hl.agg.sum(ht.exp),
            ),
        )
    )
    norm_lit = hl.literal(
        {c: {v: n[v] / n._exp for v in parts} for c, n in norms.items()}
    )
    ht = ht.annotate(
        **{f"f_{v}": ht[f"raw_{v}"] / norm_lit[ht.context][v] for v in parts},
        **{f"f_{v}_oof": ht[f"raw_{v}"] / norm_lit[ht.context][v] for v in parts},
    ).drop(*[f"raw_{v}" for v in parts])
    cols = [c for c in ht.row if c.startswith("f_")]
    cells = ht.group_by("transcript", "win", "context").aggregate(
        _exp=hl.agg.sum(ht.exp),
        **{c: hl.agg.sum(ht.exp * ht[c]) for c in cols},
    )
    cells = cells.annotate(
        **{
            f"rr_{c[2:]}": hl.or_missing(cells._exp > 0, cells[c] / cells._exp)
            for c in cols
        }
    ).drop("_exp", *cols)
    cells = cells.checkpoint(hl.utils.new_temp_file("rate_cells.cells", "ht"))

    base = hl.import_table(args.cells, force_bgz=True, min_partitions=64)
    base = base.key_by("transcript", "win", "context")
    base = base.annotate(**cells[base.key])
    print(
        "Base cells without a factor:",
        base.aggregate(hl.agg.count_where(hl.is_missing(base.rr_roulette_only))),
        "of",
        base.count(),
    )
    base.export(args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", default=CELLS)
    parser.add_argument(
        "--out", default="gs://gnomad-tmp/jg/rr_downstream/syn_cells_rr_rates.tsv.bgz"
    )
    main(parser.parse_args())
