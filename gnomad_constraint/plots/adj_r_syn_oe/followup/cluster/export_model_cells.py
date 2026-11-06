"""Export variant cells for trying a single expected-count model.

From the v4.1.1 preprocessed context Table: every possible autosomal SNV in an Ensembl
MANE Select transcript that the pipeline applies the models to, of the consequence sets
in ``--csqs`` (the pipeline's constraint groups: ``syn`` = synonymous_variant, ``mis`` =
missense_variant, ``lof`` = LOFTEE HC, the release's ``lof``), grouped into cells of
(consequence set, transcript, 1-kb window, context, ref, alt, methylation level,
coverage). Per cell:

- ``possible`` and ``observed`` (the calibration counts: PASS, adj AC > 0, AF below the
  pipeline's cutoff);
- ``mu`` (the mutation rate of the group) and ``cpg``;
- ``e_old``: the current expected count, possible * stored plateau line * coverage
  correction, as ``apply_models`` computes it;
- ``dnm``: synonymous DNMs at those variants (the full run's DNM Table);
- ``m_sum``: the sum over possible variants of Roulette's multiplier (1 where it has
  none), for a cell-mean offset.

Unless ``--skip-window-mq``, also writes ``window_mq.tsv.bgz`` next to the cells: per 1-kb window, the mean
allele-specific mapping quality (``info.AS_MQ``) of the gnomAD v4.1 exome SNV sites
in it, their number, and the share with AS_MQ below 50. Mapping quality is only known
where a variant was called, so it is summarized over the window rather than used at
the variant itself, which would leak whether the variant was seen.

The window id is ``chrom-start-end`` with start = position // 1000 * 1000, as in the
exported cells and window features.
"""

import argparse
import logging

import hail as hl

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import MU_GROUPING

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

ROOT = "gs://gnomad/v4.1.1/constraint"
DNM = "gs://gnomad-tmp/jg/rr_full/gnomad_v4.1.1/constraint/regional_rate/gnomad.v4.1.1.syn_dnm.ht"
ROULETTE = "gs://gnomad-tmp/jg/rr_roulette/gnomad_v4.1.1/constraint/regional_rate/gnomad.v4.1.1.roulette_variant_rates.ht"


def main(args: argparse.Namespace) -> None:
    """Write the cells TSV."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day", log="/tmp/export_model_cells.log")
    plateau = hl.eval(
        hl.experimental.read_expression(f"{ROOT}/models/gnomad.v4.1.1.plateau.he")
    )
    coverage = hl.eval(
        hl.experimental.read_expression(f"{ROOT}/models/gnomad.v4.1.1.coverage.he")
    )
    lines = {
        k.cpg: (float(v[0][0]), float(v[0][1]))
        for k, v in plateau.items()
        if k.genomic_region == "autosome_or_par"
    }
    logger.info("Plateau lines: %s; coverage model: %s", lines, coverage)

    mu_rows = (
        constraint_res.get_mutation_ht(version="4.1.1")
        .ht()
        .key_by()
        .select(*MU_GROUPING, "mu_snp")
        .collect()
    )
    mu_map = hl.literal(
        {tuple(r[f] for f in MU_GROUPING): float(r.mu_snp) for r in mu_rows}
    )

    ht = hl.read_table(
        f"{ROOT}/preprocessed_data/gnomad.v4.1.1.context.preprocessed.ht"
    )
    if args.test:
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chr22", "GRCh38")])
    ht = ht.filter(ht.locus.in_autosome())
    cal = ht.calibrate_mu
    ht = ht.filter(hl.is_defined(cal.apply_model) & (cal.possible_variants > 0))
    # x.lof is missing outside pLoF consequences, and a missing condition would make
    # the whole case missing.
    csq_of = (
        lambda x: hl.case()
        .when(hl.or_else(x.lof == "HC", False), "lof")
        .when(x.most_severe_consequence == "missense_variant", "mis")
        .when(x.most_severe_consequence == "synonymous_variant", "syn")
        .or_missing()
    )
    tc = ht.vep.transcript_consequences.filter(
        lambda x: hl.is_defined(x.mane_select)
        & x.transcript_id.startswith("ENST")
        & hl.literal(set(args.csqs)).contains(csq_of(x))
    )
    ht = ht.annotate(
        transcript=hl.or_missing(hl.len(tc) > 0, tc[0].transcript_id),
        csq=hl.or_missing(hl.len(tc) > 0, csq_of(tc[0])),
    )
    ht = ht.filter(hl.is_defined(ht.transcript))

    cal = ht.calibrate_mu
    mu = mu_map.get(hl.tuple([ht[f] for f in MU_GROUPING]))
    line = hl.literal(lines)[ht.cpg]
    c = ht.exomes_coverage
    h = (
        hl.case()
        .when(c == 0, 0.0)
        .when(cal.apply_model.high_or_low_coverage == "high", 1.0)
        .default(coverage[0] + coverage[1] * c)
    )
    start = ht.locus.position // 1000 * 1000
    win_expr = ht.locus.contig + "-" + hl.str(start) + "-" + hl.str(start + 1000)
    dnm = hl.read_table(DNM)
    rou = hl.read_table(ROULETTE)
    ht = ht.annotate(
        win=win_expr,
        cov=c,
        mu=mu,
        e_old=cal.possible_variants * (line[0] + line[1] * mu) * h,
        observed=cal.observed_variants[0],
        possible=cal.possible_variants,
        dnm=hl.int(hl.is_defined(dnm[ht.key])),
        m=hl.or_else(rou[ht.key].variant_rate_multiplier, 1.0),
    )
    ht = ht.filter(hl.is_defined(ht.mu))
    cells = ht.group_by(
        "csq",
        "transcript",
        "win",
        "context",
        "ref",
        "alt",
        "methylation_level",
        "cpg",
        coverage=ht.cov,
    ).aggregate(
        possible=hl.agg.sum(ht.possible),
        observed=hl.agg.sum(ht.observed),
        mu=hl.agg.take(ht.mu, 1)[0],
        e_old=hl.agg.sum(ht.e_old),
        dnm=hl.agg.sum(ht.dnm),
        m_sum=hl.agg.sum(ht.m * ht.possible),
    )
    cells.export(args.out)
    logger.info("Wrote %s", args.out)
    if args.skip_window_mq:
        return

    sites = constraint_res.get_sites_resource("exomes", version="4.1.1").ht()
    if args.test:
        sites = hl.filter_intervals(sites, [hl.parse_locus_interval("chr22", "GRCh38")])
    sites = sites.filter(
        sites.locus.in_autosome() & hl.is_snp(sites.alleles[0], sites.alleles[1])
    )
    w0 = sites.locus.position // 1000 * 1000
    sites = sites.select(
        win=sites.locus.contig + "-" + hl.str(w0) + "-" + hl.str(w0 + 1000),
        mq=sites.info.AS_MQ,
    )
    mq = sites.group_by("win").aggregate(
        n_snv=hl.agg.count(),
        mean_as_mq=hl.agg.mean(sites.mq),
        frac_mq_below_50=hl.agg.fraction(sites.mq < 50),
    )
    mq_out = args.out.replace(".tsv.bgz", ".window_mq.tsv.bgz")
    mq.export(mq_out)
    logger.info("Wrote %s", mq_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", default="gs://gnomad-tmp/jg/one_model/syn_model_cells.tsv.bgz"
    )
    parser.add_argument("--test", action="store_true", help="chr22 only.")
    parser.add_argument(
        "--csqs", nargs="+", default=["syn"], choices=["syn", "mis", "lof"]
    )
    parser.add_argument("--skip-window-mq", action="store_true")
    main(parser.parse_args())
