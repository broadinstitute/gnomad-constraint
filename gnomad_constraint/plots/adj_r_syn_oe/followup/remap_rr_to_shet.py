"""Remap synonymous obs/exp under each rr version onto s_het and LOEUF deciles.

Uses Ruchit Panchal's rr evaluation export
(gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export/), which has
one row per (Ensembl MANE Select transcript, 1kb window, context) cell with the
unadjusted synonymous expected and observed counts. Each cell is multiplied by the rr
of its window and context from each rr Table, then cells are summed per transcript and
per decile.

Steps:

- ``--annotate-rr`` (Hail, local): add one rr column per rr Table to the cells.
- ``--summarize`` (pandas): sum obs/exp per transcript, bin by gnomAD v4 s_het
  (Kar et al. 2026, Supplementary Table 3) and by v4.1.2 LOEUF decile.
- ``--outliers`` (pandas): count synonymous outliers by direction and release gene
  flag.
"""

import argparse
import logging
import re

import hail as hl
import numpy as np
import pandas as pd

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

RR_TABLES = {
    "adj_r": "legacy/adj_r_per_context_methyl_genome_1kb_autosome.agg.ht",
    "syn_adj_r": "legacy/adj_r_syn_dnm_per_context_methyl_genome_1kb_autosome.agg.ht",
    "ruchit_17k": "v1_ddd_asd/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_bonferroni.xy_masked.calibrated.ht",
    "ruchit_satterstrom": "satterstrom/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_bonferroni.xy_masked.calibrated.ht",
}
"""rr Tables, relative to ``--rr-dir``. Each is keyed by window interval with an
``adj_r`` dict of context to rr."""


def annotate_rr(cells_path: str, rr_dir: str, out_path: str, tmp_dir: str) -> None:
    """
    Add the rr of each cell's window and context from every rr Table.

    A window's interval in every rr Table starts at the window's BED start, so the
    locus at BED start + 1 is inside it. Cells without a window get a missing rr.

    :param cells_path: Path to ``syn_cells.tsv.bgz``.
    :param rr_dir: Directory holding the rr Tables in ``RR_TABLES``.
    :param out_path: Output TSV path.
    :param tmp_dir: Hail temporary directory.
    :return: None.
    """
    hl.init(local="local[8]", tmp_dir=tmp_dir, quiet=True)
    ht = hl.import_table(
        cells_path,
        types={"exp": hl.tfloat64, "obs": hl.tint64, "n_sites": hl.tint64},
        force_bgz=True,
        min_partitions=32,
    )
    win = ht.win.split("-")
    ht = ht.annotate(
        _locus=hl.locus(win[0], hl.int(win[1]) + 1, reference_genome="GRCh38")
    )
    for name, path in RR_TABLES.items():
        rr_ht = hl.read_table(f"{rr_dir}/{path}").select("adj_r")
        ht = ht.annotate(**{f"rr_{name}": rr_ht[ht._locus].adj_r.get(ht.context)})

    ht.drop("_locus").export(out_path)


def read_shet(path: str) -> pd.DataFrame:
    """
    Read gnomAD v4 LoF s_het (Kar et al. 2026, Supplementary Table 3).

    :param path: Path to the xlsx.
    :return: DataFrame with ``gene_id`` and ``shet``.
    """
    shet = pd.read_excel(path, sheet_name="LoF_s_het")
    return shet.rename(columns={"Posterior_mean": "shet"})[["gene_id", "shet"]]


def per_transcript(cells_rr_path: str) -> pd.DataFrame:
    """
    Sum syn obs and expected under each rr per transcript.

    Every ``rr_*`` column of the cells file is used (e.g. ``rr_recommended`` from
    ``rr_shrinkage_experiments.py``). A missing rr (no window) counts as 1.
    ``syn_adj_r_norm`` is the legacy table divided by its expected-weighted mean per
    context over all exported cells, the same normalization as Ruchit's calibrated
    Tables.

    :param cells_rr_path: Output of :func:`annotate_rr`.
    :return: DataFrame indexed by transcript with ``obs`` and one expected column per
        correction (``none`` is unadjusted).
    """
    cells = pd.read_csv(cells_rr_path, sep="\t", na_values=["NA"])
    rr_cols = [c for c in cells.columns if c.startswith("rr_")]
    cells[rr_cols] = cells[rr_cols].fillna(1.0)

    norm = (
        cells.assign(w=cells.exp * cells.rr_syn_adj_r)
        .groupby("context")
        .apply(lambda d: d.w.sum() / d.exp.sum(), include_groups=False)
    )
    cells["rr_syn_adj_r_norm"] = cells.rr_syn_adj_r / cells.context.map(norm)
    rr_cols.append("rr_syn_adj_r_norm")

    exp_cols = {"none": cells.exp}
    for c in rr_cols:
        exp_cols[c[3:]] = cells.exp * cells[c]
        logger.info(
            "%s: expected-weighted mean rr over all cells %.4f",
            c,
            exp_cols[c[3:]].sum() / cells.exp.sum(),
        )

    return (
        pd.DataFrame({"transcript": cells.transcript, "obs": cells.obs, **exp_cols})
        .groupby("transcript")
        .sum()
    )


def read_release_flags(release_tsv: str) -> pd.DataFrame:
    """
    Read gene ID and release gene flags per Ensembl transcript.

    :param release_tsv: v4.1.2 release TSV (bgz).
    :return: DataFrame with ``transcript``, ``gene``, ``gene_id`` and ``gene_flags``
        (a list).
    """
    rel = pd.read_csv(
        release_tsv,
        sep="\t",
        usecols=["gene", "gene_id", "transcript", "gene_flags"],
        compression="gzip",
        na_values=["NA"],
    )
    rel = rel[rel.transcript.str.startswith("ENST")].drop_duplicates("transcript")
    rel["gene_flags"] = (
        rel.gene_flags.fillna("")
        .str.strip("[]")
        .str.replace('"', "")
        .str.split(",")
        .map(lambda f: [x for x in f if x])
    )
    return rel


def summarize(
    cells_rr_path: str,
    transcripts_path: str,
    release_tsv: str,
    shet_path: str,
    out_path: str,
) -> None:
    """
    Sum syn obs and exp under each rr per s_het and LOEUF decile.

    Kept transcripts are autosomal ones with a LOEUF decile and an s_het, and without
    a release QC gene flag (``low_exome_mapping_quality``, ``low_exome_coverage``).
    Constraint outlier flags are not used: they come from the uncorrected z-scores
    being evaluated, so dropping them would flatten the deciles by construction. The
    rr Tables are autosome-only or chrX/Y-masked, so chrX/Y are left out.

    :param cells_rr_path: Output of :func:`annotate_rr`.
    :param transcripts_path: Path to ``transcripts.tsv.bgz``.
    :param release_tsv: v4.1.2 release TSV, used to map transcript to gene ID.
    :param shet_path: Path to the s_het xlsx.
    :param out_path: Output summary TSV path.
    :return: None.
    """
    per_tx = per_transcript(cells_rr_path)
    exp_cols = [c for c in per_tx.columns if c != "obs"]

    tx = pd.read_csv(transcripts_path, sep="\t", na_values=["NA"], compression="gzip")
    check = per_tx.join(tx.set_index("transcript")[["obs", "exp_unadj"]], rsuffix="_tx")
    logger.info(
        "Cell sums vs transcript table: max |obs diff| %s, max |exp diff| %.3g",
        (check.obs - check.obs_tx).abs().max(),
        (check["none"] - check.exp_unadj).abs().max(),
    )

    tx = tx[tx.is_xy.astype(str).str.lower().eq("false") & tx.decile.notna()]
    rel = read_release_flags(release_tsv)
    rel = rel[rel.gene_flags.map(len) == 0][["transcript", "gene_id"]]
    tx = tx.merge(rel, on="transcript").merge(read_shet(shet_path), on="gene_id")
    # Decile 0 = highest s_het (most constrained), as in plot_syn_oe_by_shet_decile.py.
    tx["shet_decile"] = pd.qcut(
        tx.shet.rank(method="first", ascending=False), 10, labels=False
    )
    tx["loeuf_decile"] = tx.decile.astype(int)
    logger.info(
        "%d transcripts kept (no QC gene flag, autosomal, with s_het).", len(tx)
    )

    df = tx[["transcript", "shet_decile", "loeuf_decile"]].merge(
        per_tx.reset_index(), on="transcript"
    )
    rows = []
    for binning in ("shet_decile", "loeuf_decile"):
        g = df.groupby(binning)
        for corr in exp_cols:
            s = g[["obs", corr]].sum()
            for dec, r in s.iterrows():
                rows.append(
                    {
                        "binning": binning.replace("_decile", ""),
                        "decile": dec,
                        "correction": corr,
                        "n": int(g.size()[dec]),
                        "obs": r.obs,
                        "exp": r[corr],
                        "oe": r.obs / r[corr],
                    }
                )
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


GENE_FAMILIES = {
    "histone genes": r"^(H1-|H2AC|H2BC|H3C|H4C|H2AW|H2BW|H3-|H4-)",
}
"""Gene families checked among synonymous outliers, by gene symbol."""


def outlier_category(gene_flags: list, gene: str) -> str:
    """
    Get the outlier category of a transcript.

    Release QC gene flags come first (``low_exome_mapping_quality``, then
    ``low_exome_coverage``, the only two in the v4.1.2 release), then the gene
    families in ``GENE_FAMILIES``, then "unflagged".

    :param gene_flags: Release gene flags of the transcript.
    :param gene: Gene symbol.
    :return: Category.
    """
    for flag in ("low_exome_mapping_quality", "low_exome_coverage"):
        if flag in gene_flags:
            return flag
    if gene_flags:
        return "other flag"
    for family, pattern in GENE_FAMILIES.items():
        if re.match(pattern, gene or ""):
            return family
    return "unflagged"


def summarize_outliers(
    cells_rr_path: str, release_tsv: str, out_path: str, z_cutoff: float = 8.0
) -> None:
    """
    Count synonymous outliers under each rr by direction and category.

    Categories are release QC gene flags, then gene families (see
    :func:`outlier_category`).

    Uses every exported transcript (autosomes and chrX/Y) and the raw z-score
    (obs - exp) / sqrt(exp), as in Ruchit's report.

    :param cells_rr_path: Output of :func:`annotate_rr`.
    :param release_tsv: v4.1.2 release TSV, for gene flags.
    :param out_path: Output TSV path.
    :param z_cutoff: Absolute raw z-score above which a transcript is an outlier.
        Default is 8.0.
    :return: None.
    """
    per_tx = per_transcript(cells_rr_path)
    rel = read_release_flags(release_tsv).set_index("transcript").reindex(per_tx.index)
    flag_group = pd.Series(
        [
            outlier_category(
                f if isinstance(f, list) else [], g if isinstance(g, str) else ""
            )
            for f, g in zip(rel.gene_flags, rel.gene)
        ],
        index=per_tx.index,
    )

    rows = []
    for corr in [c for c in per_tx.columns if c != "obs"]:
        z = (per_tx.obs - per_tx[corr]) / per_tx[corr] ** 0.5
        direction = pd.Series("", index=z.index)
        direction[z > z_cutoff] = "excess"
        direction[z < -z_cutoff] = "deficit"
        out = pd.DataFrame({"direction": direction, "flag": flag_group})
        out = out[out.direction != ""]
        for (d, f), n in out.value_counts().items():
            rows.append({"correction": corr, "direction": d, "flag": f, "n": n})
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


HISTONE_CLUSTER_6P22 = ("chr6", 25_000_000, 28_500_000)
"""Bounds of the chr6p22 replication-dependent histone gene cluster (GRCh38)."""


def write_syn_z(cells_rr_path: str, release_tsv: str, out_path: str) -> None:
    """
    Write per-transcript synonymous counts and z-scores with a histone group label.

    Groups: histone genes in the chr6p22 cluster, other histone genes, non-histone
    genes in the same chr6 region, and every other gene. Non-histone genes with a
    release QC gene flag are left out.

    :param cells_rr_path: Output of :func:`annotate_rr` (any ``rr_*`` columns).
    :param release_tsv: v4.1.2 release TSV, for gene names, positions and flags.
    :param out_path: Output TSV path.
    :return: None.
    """
    per_tx = per_transcript(cells_rr_path)
    rel = pd.read_csv(
        release_tsv,
        sep="\t",
        compression="gzip",
        na_values=["NA"],
        usecols=[
            "gene",
            "transcript",
            "chromosome",
            "start_position",
            "cds_length",
            "gene_flags",
        ],
    )
    rel = rel[rel.transcript.str.startswith("ENST")].drop_duplicates("transcript")
    d = per_tx.join(rel.set_index("transcript"))
    histone = d.gene.fillna("").str.match(GENE_FAMILIES["histone genes"])
    contig, start, end = HISTONE_CLUSTER_6P22
    in_cluster = (d.chromosome == contig) & d.start_position.between(start, end)
    qc = d.gene_flags.fillna("[]").str.strip("[]").str.len() > 0
    d["group"] = np.select(
        [histone & in_cluster, histone, qc, in_cluster],
        [
            "histone, chr6p22 cluster",
            "histone, elsewhere",
            "",
            "non-histone, same chr6 region",
        ],
        default="other genes",
    )
    out = d[d.group != ""][["gene", "group", "cds_length", "obs"]].copy()
    for corr in [c for c in per_tx.columns if c != "obs"]:
        out[f"exp_{corr}"] = d[corr]
        out[f"z_{corr}"] = (d.obs - d[corr]) / np.sqrt(d[corr])
    out.to_csv(out_path, sep="\t")


def main(args: argparse.Namespace) -> None:
    """Run the requested steps."""
    if args.annotate_rr:
        annotate_rr(args.cells, args.rr_dir, args.cells_rr, args.tmp_dir)
    if args.summarize:
        summarize(
            args.cells_rr, args.transcripts, args.release_tsv, args.shet, args.out
        )
    if args.outliers:
        summarize_outliers(args.cells_rr, args.release_tsv, args.outliers_out)
    if args.syn_z:
        write_syn_z(args.cells_rr, args.release_tsv, args.syn_z_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotate-rr", action="store_true")
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--outliers", action="store_true")
    parser.add_argument("--syn-z", action="store_true")
    parser.add_argument("--cells", help="Path to syn_cells.tsv.bgz.")
    parser.add_argument("--transcripts", help="Path to transcripts.tsv.bgz.")
    parser.add_argument("--rr-dir", help="Directory holding the rr Tables.")
    parser.add_argument("--cells-rr", help="Cells TSV with rr columns.")
    parser.add_argument("--release-tsv", help="v4.1.2 release TSV (bgz).")
    parser.add_argument("--shet", help="Kar et al. 2026 Supplementary Table 3 xlsx.")
    parser.add_argument("--tmp-dir", default="/tmp")
    parser.add_argument("--out", help="Output summary TSV.")
    parser.add_argument("--outliers-out", help="Output outlier count TSV.")
    parser.add_argument("--syn-z-out", help="Output per-transcript syn z TSV.")
    main(parser.parse_args())
