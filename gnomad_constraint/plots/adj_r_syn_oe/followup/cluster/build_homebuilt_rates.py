"""Build Roulette-style per-variant rate multipliers from gnomAD v4.1 genomes ourselves.

Roulette (Seplyarskiy et al. 2023) learns a mutation rate for every possible SNV from
gnomAD rare variants. Two of its parts carry what our trinucleotide mutation rates
and regional features miss, and this rebuilds both from gnomAD v4.1 genome rare
variants (PASS, adj AC > 0, AF < ``--max-af``) outside exons, low-complexity regions
and segmental duplications, on the autosomes:

- a regional rate: the observed share of possible variants of each trinucleotide
  change in each 50 kb window, relative to that change's genome-wide share. As in
  Roulette, it is per trinucleotide change; here it is shrunk toward its window's
  rate for the change's class, and that toward the window's overall rate, so that
  sparse cells don't add noise;
- a sequence-context part: the genome-wide observed share of each heptamer change,
  relative to its trinucleotide change's share.

A coding variant's rate is the product of the two. ``get_variant_rate_ht`` then
divides it by the mean of its mutation-rate group (context, ref, alt, methylation
level) over coding sites and leaves CpG contexts at 1, exactly as for Roulette's own
rates, so the result plugs into ``regional_rate_pipeline.py --use-variant-rates
--variant-rate-ht-path`` and ``constraint_pipeline.py --variant-rate-ht-path``.

Steps (each skipped if its output exists):

1. ``--count-reference``: possible sites per (window, trinucleotide) and per heptamer,
   counted from hg38.2bit with numpy on the driver (reference k-mer counts only).
2. ``--count-variants``: observed rare variants per (window, trinucleotide, alt) and
   per (heptamer, alt), in Hail.
3. ``--build-rates``: the multipliers, and the per-variant coding rate Table.
"""

import argparse
import logging
import os
import struct
import subprocess

import hail as hl
import numpy as np
import pandas as pd

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.resources.constants import MU_GROUPING
from gnomad_constraint.utils.regional_rate import get_variant_rate_ht

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

TWOBIT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
FASTA = "gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz"
FAI = "gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai"
GENOMES_SITES = (
    "gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht"
)
LCR = "gs://gcp-public-data--gnomad/resources/grch38/lcr_intervals/LCRFromHengHg38.ht"
SEGDUP = (
    "gs://gcp-public-data--gnomad/resources/grch38/seg_dup_intervals/GRCh38_segdups.ht"
)
AUTOSOMES = [f"chr{c}" for c in range(1, 23)]
CODE = np.frombuffer(b"TCAG", dtype=np.uint8)
LETTERS = np.array(list("ACGT"))
BASE_CODE = np.full(256, 4, dtype=np.int8)
BASE_CODE[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
FLANK = 3
CLASSES = {
    ("A", "G"): "A>G",
    ("A", "C"): "A>C/T",
    ("A", "T"): "A>C/T",
    ("C", "T"): "C>T",
    ("C", "A"): "C>A/G",
    ("C", "G"): "C>A/G",
}


def read_chromosome(path: str, name: str) -> np.ndarray:
    """
    Read one sequence of a .2bit file as an array of ASCII bases (N blocks as N).

    :param path: .2bit path.
    :param name: Sequence name, e.g. "chr1".
    :return: uint8 array of bases.
    """
    with open(path, "rb") as f:
        _, _, n_seq, _ = struct.unpack("<IIII", f.read(16))
        offsets = {}
        for _ in range(n_seq):
            nm = f.read(ord(f.read(1))).decode()
            offsets[nm] = struct.unpack("<I", f.read(4))[0]
        f.seek(offsets[name])
        size, n_n = struct.unpack("<II", f.read(8))
        n_starts = np.frombuffer(f.read(4 * n_n), dtype="<u4")
        n_sizes = np.frombuffer(f.read(4 * n_n), dtype="<u4")
        (n_mask,) = struct.unpack("<I", f.read(4))
        f.seek(8 * n_mask + 4, 1)
        packed = np.frombuffer(f.read((size + 3) // 4), dtype=np.uint8)
    codes = np.empty(len(packed) * 4, dtype=np.uint8)
    for i, shift in enumerate((6, 4, 2, 0)):
        codes[i::4] = (packed >> shift) & 3
    seq = CODE[codes[:size]].copy()
    for s, n in zip(n_starts, n_sizes):
        seq[s : s + n] = ord("N")
    return seq


def exclusion_ht(version: str, padding: int) -> hl.Table:
    """
    Get the regions left out of the counts: exons (padded), LCRs and segdups.

    :param version: Constraint resource version, for the GENCODE Table.
    :param padding: Bases added to each side of each exon.
    :return: Table keyed by ``interval``.
    """
    gc = constraint_res.get_gencode_ht(version)
    gc = gc.filter(gc.feature == "exon")
    gc = gc.key_by(
        interval=hl.locus_interval(
            gc.interval.start.contig,
            hl.max(gc.interval.start.position - padding, 1),
            gc.interval.end.position + padding,
            includes_end=True,
            reference_genome="GRCh38",
            invalid_missing=True,
        )
    ).select()
    gc = gc.filter(hl.is_defined(gc.interval))
    parts = [gc] + [hl.read_table(p).key_by("interval").select() for p in (LCR, SEGDUP)]
    return parts[0].union(*parts[1:])


def code7(codes: np.ndarray, pos: np.ndarray) -> tuple:
    """
    Get the pyrimidine-centred heptamer and trinucleotide code of positions.

    :param codes: Base codes (0-3 = A, C, G, T; 4 = N) of the chromosome.
    :param pos: 0-based positions to code (their heptamer has no N).
    :return: (heptamer codes, trinucleotide codes), heptamers as base-4 numbers of
        the seven bases, after reverse-complementing positions whose middle base is
        G or T.
    """
    win = np.stack([codes[pos + k] for k in range(-FLANK, FLANK + 1)], axis=1).astype(
        np.int32
    )
    flip = win[:, FLANK] >= 2
    win[flip] = (3 - win[flip])[:, ::-1]
    weights = 4 ** np.arange(2 * FLANK, -1, -1)
    hept = win @ weights
    tri = win[:, FLANK - 1 : FLANK + 2] @ np.array([16, 4, 1])
    return hept, tri


def count_reference(args: argparse.Namespace) -> None:
    """Count possible sites per (window, trinucleotide) and per heptamer."""
    twobit = "/tmp/hg38.2bit"
    if not os.path.exists(twobit):
        subprocess.run(["curl", "-sSfL", TWOBIT_URL, "-o", twobit], check=True)
    excl = exclusion_ht(args.version, args.exon_padding)
    rows = excl.aggregate(
        hl.agg.collect(
            (
                excl.interval.start.contig,
                excl.interval.start.position,
                excl.interval.end.position,
            )
        )
    )
    excl_df = pd.DataFrame(rows, columns=["contig", "start", "end"])
    windows, hepts = [], np.zeros(4**7, dtype=np.int64)
    for chrom in args.contigs:
        seq = read_chromosome(twobit, chrom)
        codes = BASE_CODE[seq]
        n = len(codes)
        # Mask exons, LCRs and segdups (1-based inclusive -> 0-based).
        diff = np.zeros(n + 1, dtype=np.int32)
        e = excl_df[excl_df.contig == chrom]
        np.add.at(diff, np.clip(e.start.to_numpy() - 1, 0, n), 1)
        np.add.at(diff, np.clip(e.end.to_numpy(), 0, n), -1)
        keep = np.cumsum(diff[:-1]) == 0
        valid = keep.copy()
        valid[:FLANK] = valid[n - FLANK :] = False
        for k in range(-FLANK, FLANK + 1):
            valid[FLANK : n - FLANK] &= codes[FLANK + k : n - FLANK + k] < 4
        nw = n // args.window + 1
        counts = np.zeros(nw * 64, dtype=np.int64)
        pos_all = np.flatnonzero(valid)
        # In chunks, to keep the 7-column arrays small.
        for c in range(0, len(pos_all), 20_000_000):
            pos = pos_all[c : c + 20_000_000]
            hept, tri = code7(codes, pos)
            hepts += np.bincount(hept, minlength=4**7)
            counts += np.bincount((pos // args.window) * 64 + tri, minlength=nw * 64)
        counts = counts.reshape(nw, 64)
        wi, ti = np.nonzero(counts)
        windows.append(
            pd.DataFrame(
                {"contig": chrom, "window": wi, "tri": ti, "possible": counts[wi, ti]}
            )
        )
        logger.info("%s: %d positions counted", chrom, len(pos_all))
    df = pd.concat(windows)
    tri_s = ["".join(LETTERS[[t // 16, (t // 4) % 4, t % 4]]) for t in range(64)]
    df["context"] = [tri_s[t] for t in df.tri]
    with hl.hadoop_open(f"{args.out}/reference_window_counts.tsv", "w") as f:
        df.drop(columns="tri").to_csv(f, sep="\t", index=False)
    nz = np.flatnonzero(hepts)
    hept_s = [
        "".join(LETTERS[[(h >> (2 * (6 - i))) & 3 for i in range(7)]]) for h in nz
    ]
    with hl.hadoop_open(f"{args.out}/reference_heptamer_counts.tsv", "w") as f:
        pd.DataFrame({"heptamer": hept_s, "possible": hepts[nz]}).to_csv(
            f, sep="\t", index=False
        )


def collapsed_context(ht: hl.Table) -> tuple:
    """
    Get the pyrimidine-centred heptamer and alt of each SNV row.

    :param ht: Table keyed by ``locus`` and ``alleles``.
    :return: (heptamer, alt, whether the heptamer is valid) expressions.
    """
    seq = hl.get_sequence(
        ht.locus.contig,
        ht.locus.position,
        before=FLANK,
        after=FLANK,
        reference_genome="GRCh38",
    ).upper()
    flip = hl.literal({"G", "T"}).contains(ht.alleles[0])
    hept = hl.if_else(flip, hl.reverse_complement(seq), seq)
    alt = hl.if_else(flip, hl.reverse_complement(ht.alleles[1]), ht.alleles[1])
    ok = (
        (hl.len(seq) == 2 * FLANK + 1)
        & (seq[FLANK] == ht.alleles[0])
        & seq.matches("^[ACGT]+$")
    )
    return hept, alt, ok


def count_variants(args: argparse.Namespace) -> None:
    """Count observed rare variants per (window, trinucleotide, alt) and heptamer."""
    g = hl.read_table(GENOMES_SITES)
    adj = hl.eval(g.freq_meta).index({"group": "adj"})
    g = g.select(freq=g.freq[adj], filters=g.filters)
    g = hl.filter_intervals(
        g, [hl.parse_locus_interval(c, reference_genome="GRCh38") for c in args.contigs]
    )
    g = g.filter(
        hl.is_snp(g.alleles[0], g.alleles[1])
        & (hl.len(g.filters) == 0)
        & (g.freq.AC > 0)
        & (g.freq.AF < args.max_af)
    )
    excl = exclusion_ht(args.version, args.exon_padding)
    g = g.filter(hl.is_missing(excl[g.locus]))
    hept, alt, ok = collapsed_context(g)
    g = g.annotate(hept=hept, alt=alt)
    g = g.filter(ok)
    g = g.checkpoint(hl.utils.new_temp_file("rare_genome_snvs", "ht"))
    w = g.group_by(
        contig=g.locus.contig,
        window=g.locus.position // args.window,
        context=g.hept[FLANK - 1 : FLANK + 2],
        alt=g.alt,
    ).aggregate(observed=hl.agg.count())
    w.write(f"{args.out}/rare_window_counts.ht", overwrite=args.overwrite)
    h = g.aggregate(hl.agg.counter(hl.tuple([g.hept, g.alt])))
    with hl.hadoop_open(f"{args.out}/rare_heptamer_counts.tsv", "w") as f:
        pd.DataFrame(
            [(k[0], k[1], v) for k, v in h.items()],
            columns=["heptamer", "alt", "observed"],
        ).to_csv(f, sep="\t", index=False)


def build_rates(args: argparse.Namespace) -> None:
    """Build the regional and heptamer multipliers, and the per-variant coding rates."""
    with hl.hadoop_open(f"{args.out}/reference_window_counts.tsv") as f:
        ref_w = pd.read_csv(f, sep="\t")
    with hl.hadoop_open(f"{args.out}/reference_heptamer_counts.tsv") as f:
        ref_h = pd.read_csv(f, sep="\t")
    with hl.hadoop_open(f"{args.out}/rare_heptamer_counts.tsv") as f:
        obs_h = pd.read_csv(f, sep="\t")
    obs_w = hl.read_table(f"{args.out}/rare_window_counts.ht").to_pandas()

    # Genome-wide share of possible variants seen, per trinucleotide change.
    ref_w = ref_w[ref_w.context.str[1].isin(["A", "C"])]
    alts = {"A": ["C", "G", "T"], "C": ["A", "G", "T"]}
    poss = pd.concat(
        [ref_w[ref_w.context.str[1] == r].assign(alt=a) for r in "AC" for a in alts[r]]
    )
    cells = poss.merge(obs_w, on=["contig", "window", "context", "alt"], how="left")
    cells["observed"] = cells.observed.fillna(0)
    tot = cells.groupby(["context", "alt"])[["observed", "possible"]].sum()
    share = (tot.observed / tot.possible).rename("share")
    cells = cells.join(share, on=["context", "alt"])
    cells["expected"] = cells.possible * cells.share
    cells["cls"] = [CLASSES[(c[1], a)] for c, a in zip(cells.context, cells.alt)]
    cpg = (cells.context.str[1] == "C") & (cells.context.str[2] == "G")
    cells["cls"] = np.where(cpg, "CpG " + cells.cls, cells.cls)

    # Window rate (outside CpG), shrunk class rate, shrunk trinucleotide rate.
    key = ["contig", "window"]
    nc = cells[~cpg].groupby(key)[["observed", "expected"]].sum()
    win_rate = (nc.observed / nc.expected).rename("win_rate")
    cells = cells.join(win_rate, on=key)
    cells["win_rate"] = cells.win_rate.fillna(1.0)
    cl = cells.groupby(key + ["cls"]).agg(
        o=("observed", "sum"), e=("expected", "sum"), w=("win_rate", "first")
    )
    cl["cls_rate"] = (cl.o + args.class_prior * cl.w) / (cl.e + args.class_prior)
    cells = cells.join(cl["cls_rate"], on=key + ["cls"])
    cells["regional"] = (cells.observed + args.tri_prior * cells.cls_rate) / (
        cells.expected + args.tri_prior
    )
    logger.info(
        "Regional multiplier quantiles by class:\n%s",
        cells.groupby("cls").regional.quantile([0.05, 0.5, 0.95]).unstack().round(3),
    )

    # Heptamer share relative to its trinucleotide change's share.
    ref_h = pd.concat(
        [
            ref_h[ref_h.heptamer.str[FLANK] == r].assign(alt=a)
            for r in "AC"
            for a in alts[r]
        ]
    )
    hp = ref_h.merge(obs_h, on=["heptamer", "alt"], how="left").fillna({"observed": 0})
    hp["context"] = hp.heptamer.str[FLANK - 1 : FLANK + 2]
    hp = hp.join(share, on=["context", "alt"])
    hp["sequence"] = (hp.observed + args.tri_prior) / (
        hp.possible * hp.share + args.tri_prior
    )
    with hl.hadoop_open(f"{args.out}/regional_multipliers.tsv.bgz", "w") as f:
        cells[key + ["context", "alt", "observed", "expected", "regional"]].to_csv(
            f, sep="\t", index=False
        )
    with hl.hadoop_open(f"{args.out}/heptamer_multipliers.tsv", "w") as f:
        hp[["heptamer", "alt", "observed", "possible", "sequence"]].to_csv(
            f, sep="\t", index=False
        )

    # Per-variant rates over coding sites.
    reg_ht = hl.Table.from_pandas(
        cells[key + ["context", "alt", "regional"]].astype({"window": "int32"})
    ).key_by("contig", "window", "context", "alt")
    seq_map = hl.literal(
        {(h, a): float(v) for h, a, v in zip(hp.heptamer, hp.alt, hp.sequence)}
    )
    ctx = constraint_res.get_annotated_context_ht(version=args.version).ht()
    ctx = hl.filter_intervals(
        ctx,
        [hl.parse_locus_interval(c, reference_genome="GRCh38") for c in args.contigs],
    )
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
    hept, alt, ok = collapsed_context(ctx)
    ctx = ctx.annotate(
        window=hl.int32(ctx.locus.position // args.window),
        hept=hl.or_missing(ok, hept),
    )
    ctx = ctx.annotate(
        regional=hl.or_else(
            reg_ht[ctx.locus.contig, ctx.window, ctx.context, ctx.alt].regional, 1.0
        ),
        sequence=hl.or_else(seq_map.get((ctx.hept, ctx.alt)), 1.0),
    )
    ctx = ctx.checkpoint(hl.utils.new_temp_file("homebuilt_coding_rates", "ht"))
    rates = ctx.select(
        roulette_mu=ctx.regional * ctx.sequence, roulette_filters=hl.empty_set(hl.tstr)
    )
    ht = get_variant_rate_ht(rates, ctx.select(*MU_GROUPING), cds)
    ht = ht.annotate(
        regional_part=ctx[ht.key].regional, sequence_part=ctx[ht.key].sequence
    )
    ht = ht.annotate_globals(
        variant_rate_globals=ht.variant_rate_globals.annotate(
            source="gnomAD v4.1 genomes rare variants (homebuilt, Roulette-style)",
            window=args.window,
            max_af=args.max_af,
            class_prior=args.class_prior,
            tri_prior=args.tri_prior,
        )
    )
    ht.write(f"{args.out}/homebuilt_variant_rates.ht", overwrite=args.overwrite)


def main(args: argparse.Namespace) -> None:
    """Run the requested steps."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")
    hl.get_reference("GRCh38").add_sequence(FASTA, FAI)
    if args.count_reference and not hl.hadoop_exists(
        f"{args.out}/reference_heptamer_counts.tsv"
    ):
        count_reference(args)
    if args.count_variants and not hl.hadoop_exists(
        f"{args.out}/rare_heptamer_counts.tsv"
    ):
        count_variants(args)
    if args.build_rates:
        build_rates(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="GCS folder for the outputs.")
    parser.add_argument("--version", default="4.1.1")
    parser.add_argument("--count-reference", action="store_true")
    parser.add_argument("--count-variants", action="store_true")
    parser.add_argument("--build-rates", action="store_true")
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=AUTOSOMES,
        help="Contigs to count and rate (e.g. chr22 for a quick test).",
    )
    parser.add_argument("--window", type=int, default=50_000)
    parser.add_argument("--max-af", type=float, default=0.001)
    parser.add_argument("--exon-padding", type=int, default=50)
    parser.add_argument("--coding-padding", type=int, default=10)
    parser.add_argument(
        "--class-prior",
        type=float,
        default=200.0,
        help="Pseudo-expected variants shrinking a window's class rate to its rate.",
    )
    parser.add_argument(
        "--tri-prior",
        type=float,
        default=50.0,
        help="Pseudo-expected variants shrinking a trinucleotide rate to its class.",
    )
    parser.add_argument("--overwrite", action="store_true")
    main(parser.parse_args())
