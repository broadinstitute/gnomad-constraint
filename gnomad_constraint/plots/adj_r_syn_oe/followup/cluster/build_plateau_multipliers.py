"""Build per-variant multipliers that swap the v4.1.1 plateau lines for another plateau.

Each coding SNV's expected count comes from the stored plateau line of its model group
(CpG transitions, or everything else): ``f_line(mu)``. To see what another plateau
model would give without rebuilding the per-variant expected Table, this writes
``variant_rate_multiplier = f_new(mu) / f_line(mu)`` for every coding SNV, so that
``constraint_pipeline.py --aggregate-per-variant-expected --skip-regional-rr
--variant-rate-ht-path`` gives exactly ``possible * f_new(mu) * coverage correction``.
The expected counts of a variant are linear in the multiplier, so this is the same as
applying the new plateau directly.

Versions (``--params``: ``plateau_params.json`` from ``write_plateau_params.py``, the fits
in ``plateau_models.fits.tsv`` and ``per_variant_expected.totals.tsv``, copied to GCS
by ``run_plateau.sh``):

- ``saturating``: share seen = 1 - exp(-(a + b * mu)), mu in units of 1e-8, one curve
  per model group, fit by binomial likelihood;
- ``logistic``: share seen = 1 / (1 + exp(-(a + b * mu))), same units and groups;
- ``by_substitution``: one weighted line per substitution type (strand-collapsed
  A>C, A>G, A>T, C>A, C>G, C>T outside CpG, and CpG transitions).

Only autosomes and PARs get a multiplier (the refits used their training sites); other
sites get 1. Sites are the annotated context Table's SNVs within ``--coding-padding``
bp of a GENCODE CDS, with mu from the mutation rate Table.
"""

import argparse
import json
import logging

import hail as hl
from gnomad.utils.constraint import annotate_mutation_type

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import (
    MU_GROUPING,
    VARIANT_RATE_MULTIPLIER_FIELD,
)

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def main(args: argparse.Namespace) -> None:
    """Write one multiplier Table per version."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day", log="/tmp/build_plateau_multipliers.log")
    with hl.hadoop_open(args.params) as f:
        params = json.load(f)

    ctx = constraint_res.get_annotated_context_ht(version=args.version).ht()
    cds = constraint_res.get_gencode_ht(args.version)
    cds = cds.filter(cds.feature == "CDS")
    cds = cds.key_by(
        interval=hl.locus_interval(
            cds.interval.start.contig,
            hl.max(cds.interval.start.position - args.coding_padding, 1),
            cds.interval.end.position + args.coding_padding,
            includes_end=True,
            reference_genome="GRCh38",
            invalid_missing=True,
        )
    ).select()
    cds = cds.filter(hl.is_defined(cds.interval))
    ctx = ctx.filter(hl.is_defined(cds[ctx.locus])).select(*MU_GROUPING)
    # The mutation rate Table has a few hundred rows, so look rates up in a literal.
    mu_rows = (
        constraint_res.get_mutation_ht(version=args.version)
        .ht()
        .key_by()
        .select(*MU_GROUPING, "mu_snp")
        .collect()
    )
    mu_map = hl.literal(
        {tuple(r[f] for f in MU_GROUPING): float(r.mu_snp) for r in mu_rows}
    )
    logger.info("%d mutation rate groups.", len(mu_rows))
    ctx = ctx.annotate(mu=mu_map.get(hl.tuple([ctx[f] for f in MU_GROUPING])))
    ctx = annotate_mutation_type(ctx)
    ctx = ctx.filter(hl.is_defined(ctx.mu))

    x = ctx.mu * 1e8
    cpg = ctx.cpg
    line = params["stored"]

    def pick(d: dict, key: hl.expr.StringExpression) -> tuple:
        lit = hl.literal({k: (float(v[0]), float(v[1])) for k, v in d.items()})
        p = lit[key]
        return p[0], p[1]

    group = hl.if_else(cpg, "CpG transition", "everything else")
    a0, b0 = pick(line, group)
    f_line = a0 + b0 * x
    flip = ~hl.literal({"A", "C"}).contains(ctx.ref)
    comp = hl.literal(COMP)
    sub = hl.if_else(
        cpg,
        "CpG transition",
        hl.if_else(flip, comp[ctx.ref], ctx.ref)
        + ">"
        + hl.if_else(flip, comp[ctx.alt], ctx.alt),
    )
    news = {}
    a, b = pick(params["saturating"], group)
    news["saturating"] = 1 - hl.exp(-(a + b * x))
    a, b = pick(params["logistic"], group)
    news["logistic"] = 1 / (1 + hl.exp(-(a + b * x)))
    a, b = pick(params["by_substitution"], sub)
    news["by_substitution"] = a + b * x

    auto = ctx.locus.in_autosome_or_par()
    ctx = ctx.annotate(
        **{
            f"m_{k}": hl.if_else(auto, v / f_line, 1.0)
            for k, v in news.items()
            if k in args.versions
        }
    )
    ctx = ctx.select(*[f"m_{k}" for k in args.versions]).checkpoint(
        hl.utils.new_temp_file("plateau_multipliers", "ht")
    )
    for k in args.versions:
        out = f"{args.out}/gnomad.v{args.version}.plateau_{k}_multiplier.ht"
        if hl.hadoop_exists(f"{out}/_SUCCESS") and not args.overwrite:
            logger.info("%s exists, skipping.", out)
            continue
        ht = ctx.select(**{VARIANT_RATE_MULTIPLIER_FIELD: ctx[f"m_{k}"]})
        ht = ht.annotate_globals(
            plateau_multiplier_globals=hl.struct(
                version=k, params=json.dumps(params), coding_padding=args.coding_padding
            )
        )
        ht.write(out, overwrite=args.overwrite)
        t = hl.read_table(out)
        logger.info(
            "%s multiplier quantiles: %s",
            k,
            t.aggregate(
                hl.agg.approx_quantiles(
                    t[VARIANT_RATE_MULTIPLIER_FIELD], [0, 0.01, 0.5, 0.99, 1]
                )
            ),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="4.1.1")
    parser.add_argument("--params", required=True, help="JSON file of the fits.")
    parser.add_argument(
        "--versions", nargs="+", default=["saturating", "logistic", "by_substitution"]
    )
    parser.add_argument("--coding-padding", type=int, default=50)
    parser.add_argument("--out", default="gs://gnomad-tmp/jg/rr_plateau")
    parser.add_argument("--overwrite", action="store_true")
    main(parser.parse_args())
