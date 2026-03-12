#!/usr/bin/env python3
"""
Create missense viewer gene HT with all 6 ProEMIS3D methods as regional constraint data.

Reads the 6 individual ProEMIS3D forward run HTs directly (which have obs/exp/oe
populated per residue), joins with gene_mapping.tsv to get gene_id, merges adjacent
residues into contiguous regions, and annotates the existing gnomAD gene HT.

The key difference from a transcript-based approach: we join by gene_id (Ensembl ID)
so that ProEMIS3D transcripts (e.g. ENST00000354534 for SCN8A) correctly map to their
gene even when the gnomAD canonical transcript differs.

Run locally with Java 11:
  JAVA_HOME=/opt/homebrew/opt/openjdk@11 python create_proemis3d_viewer_ht.py

Submit to Dataproc via hailctl:
  hailctl dataproc submit CLUSTER_NAME \\
    gnomad_constraint/experimental/proemis3d/create_proemis3d_viewer_ht.py \\
    --overwrite
"""
import argparse
import logging

import hail as hl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Input paths
GENE_MAPPING_PATH = "gs://gnomad-julia/proemis3d/gene_mapping.tsv"
GENCODE_POS_HT_PATH = "gs://gnomad/v4.1/constraint/proemis3d/preprocessed_data/gencode_positions.ht"
GENES_HT_PATH = "gs://gnomad-julia/constraint/genes_v4.add_promis3d.ht"

# Output path
OUTPUT_HT_PATH = "gs://gnomad-julia/constraint/genes_v4.proemis3d_all_methods.ht"

# Mapping from output field name → individual forward run HT path
METHODS = {
    "proemis3d_aic": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.aic.min_exp_mis_16.ht"
    ),
    "proemis3d_lrt": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.lrt.min_exp_mis_16.ht"
    ),
    "proemis3d_plddt_exclude_aic": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.aic.min_exp_mis_16"
        ".plddt_cutoff_70.0.plddt_method_exclude_low_plddt_from_stats.ht"
    ),
    "proemis3d_plddt_exclude_lrt": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.lrt.min_exp_mis_16"
        ".plddt_cutoff_70.0.plddt_method_exclude_low_plddt_from_stats.ht"
    ),
    "proemis3d_pae_filter_exclude_aic": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.aic.min_exp_mis_16"
        ".pae_cutoff_15.0.pae_method_filter_on_pairwise_pae_with_center.ht"
    ),
    "proemis3d_pae_filter_exclude_lrt": (
        "gs://gnomad/v4.1/constraint/proemis3d/output/"
        "proemis3D_forward.oe_upper_gamma.lrt.min_exp_mis_16"
        ".pae_cutoff_15.0.pae_method_filter_on_pairwise_pae_with_center.ht"
    ),
}


def prepare_pos_ht(raw_pos_ht: hl.Table) -> hl.Table:
    """
    Build a (uniprot_id, transcript_id) → locus_by_aapos lookup table.

    :param raw_pos_ht: Raw gencode positions HT with columns: uniprot_id, enst, gene,
        aalength, cds_len, strand, aapos, locus.
    :return: HT keyed by (uniprot_id, enst) with field locus_by_aapos: dict<int, locus>.
    """
    pos_ht = raw_pos_ht.key_by(
        "uniprot_id", "enst", "gene", "aalength", "cds_len", "strand", "aapos"
    )
    pos_ht = pos_ht.select("locus")
    pos_ht = pos_ht.collect_by_key("locus")
    pos_ht = pos_ht.annotate(
        locus=hl.sorted(pos_ht.locus.locus, key=lambda x: x.position)[0]
    )
    pos_ht = pos_ht.group_by(
        "uniprot_id", "enst", "gene", "aalength", "cds_len", "strand"
    ).aggregate(locus_by_aapos=hl.dict(hl.agg.collect((pos_ht.aapos, pos_ht.locus))))
    return pos_ht.key_by("uniprot_id", "enst").persist()


def _chisq(obs, exp):
    """Compute OE≠1 chi-square statistic."""
    return hl.if_else(obs == 0, 0.0, 2.0 * (obs * hl.log(obs / exp) - (obs - exp)))


def create_method_regions_ht(
    forward_ht: hl.Table, gene_map_ht: hl.Table, pos_ht: hl.Table
) -> hl.Table:
    """
    Convert a per-residue ProEMIS3D forward run HT into a regional constraint HT.

    Uses gene_id (via gene_map_ht) as the final key so that the result can be joined
    to the gnomAD gene HT regardless of canonical transcript differences.

    :param forward_ht: Per-residue forward run HT with fields: uniprot_id, transcript_id,
        residue_index, region_index, obs, exp, oe, oe_upper, is_null, and optionally more.
    :param gene_map_ht: HT keyed by transcript_id with field gene_id (Ensembl ID).
    :param pos_ht: Prepared positions HT keyed by (uniprot_id, enst) with locus_by_aapos.
    :return: HT keyed by gene_id with fields: uniprot_id, regions array.
    """
    # Compute chi-square and p-value
    chisq_expr = _chisq(forward_ht.obs, forward_ht.exp)
    ht = forward_ht.annotate(
        chisq=chisq_expr,
        p_value=hl.pchisqtail(chisq_expr, 1),
    )

    # Join with gene_map to get gene_id; filter to only ProEMIS3D-analysed transcripts
    ht = ht.annotate(gene_id=gene_map_ht[ht.transcript_id].gene_id)
    ht = ht.filter(hl.is_defined(ht.gene_id))

    # Group by (gene_id, uniprot_id, transcript_id, region_index, is_null)
    # and collect all residues in each group
    ht = ht.key_by("gene_id", "uniprot_id", "transcript_id", "region_index", "is_null")
    ht = ht.collect_by_key("pos")

    # Sort residues by residue_index within each group
    ht = ht.annotate(pos=hl.sorted(ht.pos, key=lambda x: x.residue_index))

    # Fold: merge adjacent residues into contiguous segments.
    # Each tuple: (aa_start, aa_stop, obs, exp, oe, oe_upper, chisq, p_value)
    ht = ht.annotate(
        pos=hl.fold(
            lambda acc, j: hl.if_else(
                j.residue_index > (acc[-1][1] + 1),
                # Gap: start a new contiguous segment
                acc.append((
                    j.residue_index, j.residue_index,
                    j.obs, j.exp, j.oe, j.oe_upper, j.chisq, j.p_value,
                )),
                # Adjacent: extend current segment's stop, use latest stats
                acc[:-1].append((
                    acc[-1][0], j.residue_index,
                    j.obs, j.exp, j.oe, j.oe_upper, j.chisq, j.p_value,
                )),
            ),
            [(
                ht.pos[0].residue_index, ht.pos[0].residue_index,
                ht.pos[0].obs, ht.pos[0].exp, ht.pos[0].oe, ht.pos[0].oe_upper,
                ht.pos[0].chisq, ht.pos[0].p_value,
            )],
            ht.pos[1:],
        )
    )

    # Explode so each contiguous segment is one row
    ht = ht.explode("pos")

    # Unpack the pos tuple into named fields
    ht = ht.key_by("gene_id", "uniprot_id", "transcript_id", "region_index", "is_null")
    ht = ht.transmute(
        aa_start=ht.pos[0],
        aa_stop=ht.pos[1],
        obs_mis=ht.pos[2],
        exp_mis=ht.pos[3],
        obs_exp=ht.pos[4],
        oe_upper=ht.pos[5],
        chisq=ht.pos[6],
        p_value=ht.pos[7],
    )

    # Collect contiguous segments into a "regions" array per
    # (gene_id, uniprot_id, transcript_id, region_index, is_null)
    ht = ht.collect_by_key("regions")

    # Look up genomic coordinates from pos_ht using (uniprot_id, transcript_id)
    ht = ht.annotate(**pos_ht[ht.uniprot_id, ht.transcript_id])

    # Map each region: convert aa positions to genomic start/stop
    ht = ht.annotate(
        regions=ht.regions.map(
            lambda x: x.select(
                chrom=ht.locus_by_aapos[x.aa_start].contig,
                start=hl.if_else(
                    ht.locus_by_aapos[x.aa_start].position
                    <= ht.locus_by_aapos[x.aa_stop].position,
                    ht.locus_by_aapos[x.aa_start].position,
                    ht.locus_by_aapos[x.aa_stop].position,
                ),
                stop=hl.if_else(
                    ht.locus_by_aapos[x.aa_start].position
                    <= ht.locus_by_aapos[x.aa_stop].position,
                    ht.locus_by_aapos[x.aa_stop].position + 2,
                    ht.locus_by_aapos[x.aa_start].position + 2,
                ),
                aa_start=x.aa_start,
                aa_stop=x.aa_stop,
                obs_mis=x.obs_mis,
                exp_mis=x.exp_mis,
                obs_exp=x.obs_exp,
                oe_upper=x.oe_upper,
                region_index=ht.region_index,
                is_null=ht.is_null,
                chisq_diff_null=x.chisq,
                p_value=x.p_value,
            )
        )
    )

    # Group by gene_id, flattening regions across all (transcript, region_index, is_null)
    ht = ht.group_by("gene_id").aggregate(
        uniprot_id=hl.agg.take(ht.uniprot_id, 1)[0],
        regions=hl.flatten(hl.agg.collect(ht.regions)),
    )
    return ht.key_by("gene_id")


def main(args):
    hl.init(
        log="/create_proemis3d_viewer_ht.log",
        tmp_dir="gs://gnomad-tmp-4day/proemis3d/viewer_ht",
    )

    logger.info("Reading gene_mapping from %s", GENE_MAPPING_PATH)
    gene_map_ht = hl.import_table(
        GENE_MAPPING_PATH, delimiter="\t", key="transcript_id"
    ).persist()

    logger.info("Preparing gencode positions HT from %s", GENCODE_POS_HT_PATH)
    raw_pos_ht = hl.read_table(GENCODE_POS_HT_PATH)
    pos_ht = prepare_pos_ht(raw_pos_ht)

    logger.info("Reading gene HT from %s", GENES_HT_PATH)
    gene_ht = hl.read_table(GENES_HT_PATH)

    # Process each of the 6 methods into a per-gene regional HT
    method_annotations = {}
    aic_method_ht = None
    for out_field, ht_path in METHODS.items():
        logger.info("Processing method: %s", out_field)
        forward_ht = hl.read_table(ht_path)
        method_ht = create_method_regions_ht(forward_ht, gene_map_ht, pos_ht)
        method_ht = method_ht.checkpoint(
            f"gs://gnomad-tmp-4day/proemis3d/viewer_ht/{out_field}.ht",
            overwrite=args.overwrite,
        )
        method_annotations[out_field] = hl.struct(
            has_no_rmc_evidence=False,
            passed_qc=True,
            regions=method_ht[gene_ht.gene_id].regions,
        )
        if out_field == "proemis3d_aic":
            aic_method_ht = method_ht

    # Annotate gene HT with all 6 method regional structs and uniprot_id
    gene_ht = gene_ht.annotate(
        uniprot_id=hl.coalesce(
            gene_ht.uniprot_id,
            aic_method_ht[gene_ht.gene_id].uniprot_id,
        ),
        **method_annotations,
    )

    logger.info("Writing output HT to %s", OUTPUT_HT_PATH)
    gene_ht = gene_ht.checkpoint(OUTPUT_HT_PATH, overwrite=args.overwrite)
    gene_ht.describe()
    logger.info("Done. Row count: %d", gene_ht.count())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create missense viewer gene HT with all 6 ProEMIS3D methods."
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing output files."
    )
    main(parser.parse_args())
