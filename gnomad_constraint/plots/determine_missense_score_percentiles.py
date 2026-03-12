"""
Script to compute percentiles for missense prediction scores and annotate constraint data.

This script:
1. Loads missense prediction scores from Genetics Gym tables
2. Computes percentiles for each score across all variants
3. Annotates constraint pipeline data with these percentiles
4. Aggregates observed/expected counts by transcript and score percentile
5. Exports aggregated results for analysis

The script supports multiple missense prediction scores including:
- ProteinMPNN, ESM, REVEL, RASP, AM, MisFit, PolyPhen, CPT1, popEVE, EVE, MPC, CADD, GPN-MSA
"""

import argparse
import logging
from typing import Dict, List, Optional

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("determine_missense_score_percentiles")
logger.setLevel(logging.INFO)

MU_GROUPING = ("context", "ref", "alt", "methylation_level")
"""Mutation rate grouping fields."""

DEFAULT_SCORES = [
    "proteinmpnn_llr_neg",
    "esm1b_neg",
    "score_PAI3D",
    "revel",
    "rasp_score",
    "AM",
    "MisFit_D",
    "MisFit_S",
    "polyphen_score",
    "cpt1_score",
    "popEVE_neg",
    "EVE",
    "ESM_1v_neg",
    "mpc",
    "cadd_score",
    "gpn_msa_score",
]
"""Default list of missense prediction scores."""

SEVERE_HI_GENES = [
    "ADNP", "AHDC1", "ANKRD11", "ARID1A", "ARID1B", "ARID2", "AUTS2",
    "CHAMP1", "CHD2", "CHD7", "CHD8", "CREBBP", "CTCF", "CTNNB1", "DMRT1",
    "DYRK1A", "EFTUD2", "EHMT1", "EP300", "FOXG1", "FOXP1", "GATA2",
    "GATA6", "GLI2", "GRIN2B", "HIVEP2", "HNRNPK", "KANSL1", "KMT2D",
    "MBD5", "MED13L", "MEF2C", "NFIA", "NIPBL", "OTX2", "PAFAH1B1", "PURA",
    "RAI1", "RERE", "SATB2", "SCN1A", "SCN2A", "SETBP1", "SETD5", "SHANK3",
    "SHH", "SIX3", "SLC2A1", "SOX5", "SOX9", "STXBP1", "SYNGAP1", "TCF4",
    "TGIF1", "ZEB2", "ZIC2", "CASK", "CDKL5", "HCCS", "KDM6A", "MECP2",
    "PCDH19", "WDR45",
]
"""List of severe haploinsufficient genes."""

MODERATE_HI_GENES = [
    "BMP4", "BMPR1A", "BMPR2", "CHRNA7", "COL11A2", "COL1A1", "COL2A1",
    "COL5A1", "ELN", "EYA1", "FAS", "FBN1", "FGF10", "FOXC1", "FOXC2",
    "FZD4", "GATA3", "GCH1", "HOXD13", "IRF6", "JAG1", "KCNH2", "KCNQ1",
    "KCNQ2", "KIF11", "LHX4", "MITF", "MNX1", "MYCN", "NF1", "NKX2-5",
    "NOG", "NRXN1", "NSD1", "PAX2", "PAX3", "PITX2", "PTEN", "RB1", "RET",
    "RPS19", "RPS24", "RPS26", "RUNX2", "SALL1", "SALL4", "SF3B4", "SMAD3",
    "SMAD4", "SMARCB1", "SOX10", "SOX2", "SPRED1", "STK11", "TBX1", "TBX3",
    "TBX5", "TCF12", "TCOF1", "TRPS1", "TSC1", "TSC2", "WT1", "ACSL4",
    "BCOR", "DCX", "EBP", "EFNB1", "FLNA", "FMR1", "GRIA3", "LAMP2",
    "NSDHL", "OFD1", "OTC", "PDHA1", "PHF6", "PLP1", "PORCN", "SLC6A8",
    "SLC9A6",
]
"""List of moderate haploinsufficient genes."""

NEW_HI_GENES = [
    "ADNP", "AHDC1", "ANKRD11", "ARID1A", "ARID1B", "ARID2", "AUTS2",
    "CHAMP1", "CHD2", "CHD7", "CHD8", "CREBBP", "CTNNB1", "DYRK1A",
    "EFTUD2", "EHMT1", "EP300", "FOXG1", "FOXP1", "GATA2", "GATA6",
    "GRIN2B", "HIVEP2", "KANSL1", "KMT2D", "MBD5", "MED13L", "MEF2C",
    "NFIA", "NIPBL", "PAFAH1B1", "PURA", "SATB2", "SCN1A", "SCN2A",
    "SETD5", "SHANK3", "SLC2A1", "SOX5", "SOX9", "STXBP1", "SYNGAP1",
    "TCF4", "ZEB2", "ZIC2", "CTCF", "HNRNPK", "RAI1", "RERE", "SETBP1",
    "ASH1L", "ASXL1", "ASXL3", "BCL11A", "CIC", "GATAD2B", "KAT6A",
    "KAT6B", "KMT2A", "KMT2C", "MYT1L", "NFIX", "OTX2", "PBX1", "PHIP",
    "POGZ", "SETD2", "SON", "SOX11", "TBL1XR1", "TBR1", "TCF20", "TRIP12",
    "WAC", "ZBTB18", "ZMYND11", "ZNF462",
]
"""List of newly identified haploinsufficient genes."""


def get_hi_genes() -> Dict[str, List[str]]:
    """
    Get haploinsufficient gene lists.

    :return: Dictionary with keys 'severe', 'moderate', 'new', and 'all'
    """
    hi_genes = list(set(SEVERE_HI_GENES + MODERATE_HI_GENES + NEW_HI_GENES))
    return {
        "severe": SEVERE_HI_GENES,
        "moderate": MODERATE_HI_GENES,
        "new": NEW_HI_GENES,
        "all": hi_genes,
    }


def preprocess_missense_scores(
    ht: hl.Table,
    scores: List[str],
    gene_id_field: Optional[str] = None,
    transcript_id_field: Optional[str] = None,
) -> hl.Table:
    """
    Preprocess missense scores from Genetics Gym table.

    Selects score fields and an ID field (gene_id or transcript_id), filters to
    variants with at least one defined score, adds a random number for
    deterministic ordering, and keys by locus, alleles, and the ID field.

    :param ht: Input Genetics Gym missense scores table.
    :param scores: List of score field names to include.
    :param gene_id_field: Field name for gene ID in the input table.
    :param transcript_id_field: Field name for transcript ID in the input table.
    :return: Preprocessed Hail table with missense scores.
    :raises ValueError: If both ID fields are None, if the specified ID field(s)
        are not found, or if any score fields are not found.
    """
    if gene_id_field is None and transcript_id_field is None:
        raise ValueError(
            "Either gene_id_field or transcript_id_field must be provided."
        )
    if gene_id_field is not None and gene_id_field not in ht.row:
        raise ValueError(
            f"Gene ID field '{gene_id_field}' not found in table. "
            f"Available fields: {list(ht.row)}"
        )
    if transcript_id_field is not None and transcript_id_field not in ht.row:
        raise ValueError(
            f"Transcript ID field '{transcript_id_field}' not found in table. "
            f"Available fields: {list(ht.row)}"
        )

    scores_not_found = [s for s in scores if s not in ht.row]
    if scores_not_found:
        raise ValueError(
            f"Score fields {scores_not_found} not found in table. "
            f"Available fields: {list(ht.row)}"
        )

    # Build the renamed ID field(s) and determine key.
    id_fields = {}
    if transcript_id_field:
        id_fields["transcript_id"] = ht[transcript_id_field]
        id_key = "transcript_id"
    else:
        id_fields["gene_id"] = ht[gene_id_field]
        id_key = "gene_id"

    # Select scores + renamed ID field.
    ht = ht.select(*scores, **id_fields)

    logger.info("Filtering to variants with at least one defined score")
    ht = ht.filter(hl.any([hl.is_defined(ht[s]) for s in scores]))

    logger.info("Adding random number for deterministic ordering")
    ht = ht.annotate(rand_n=hl.rand_unif(0, 1))

    key_fields = ["locus", "alleles", id_key]
    logger.info("Keying by: %s", ", ".join(key_fields))
    ht = ht.key_by(*key_fields).distinct()
    ht = ht.select(*ht.row_value)

    return ht


def preprocess_constraint_table(ht: hl.Table) -> hl.Table:
    """
    Preprocess a constraint pipeline table.

    Flattens calibrate_mu struct fields to top-level and filters to missense
    variants.

    :param ht: Raw constraint table from the constraint pipeline.
    :return: Preprocessed table filtered to missense variants.
    """
    logger.info("Flattening calibrate_mu fields")
    ht = ht.annotate(**ht.calibrate_mu)

    logger.info("Filtering to missense variants")
    ht = ht.filter(ht.annotation == "missense_variant")

    return ht


def compute_score_percentiles(
    ht: hl.Table,
    scores: List[str],
    checkpoint_prefix: Optional[str] = None,
    overwrite: bool = True,
) -> hl.Table:
    """
    Compute percentiles for each missense prediction score.

    Orders variants by score value (with random tie-breaker), assigns an index,
    and converts to a percentile (0-100).

    :param ht: Preprocessed table with missense scores.
    :param scores: List of score field names.
    :param checkpoint_prefix: If provided, checkpoint each per-score intermediate
        table to ``{checkpoint_prefix}_{score}.ht``.
    :param overwrite: Whether to overwrite existing checkpoints.
    :return: Table with percentile annotations for each score.
    """
    logger.info("Computing percentiles for %d scores", len(scores))

    hts = {}
    for s in scores:
        logger.info("  Processing score: %s", s)
        tmp_ht = ht.select("rand_n", **{s: hl.struct(score=ht[s])})
        tmp_ht = tmp_ht.order_by(tmp_ht[s].score, tmp_ht.rand_n)
        tmp_ht = tmp_ht.annotate(
            **{
                s: tmp_ht[s].annotate(
                    idx=hl.or_missing(
                        hl.is_defined(tmp_ht[s].score),
                        hl.scan.count_where(hl.is_defined(tmp_ht[s].score)),
                    )
                )
            }
        ).key_by(*ht.key)

        if checkpoint_prefix:
            tmp_ht = tmp_ht.checkpoint(
                f"{checkpoint_prefix}_{s}.ht", overwrite=overwrite
            )
        hts[s] = tmp_ht

    # Join indices back to original table.
    ht = ht.select(**{s: _ht[ht.key][s] for s, _ht in hts.items()})

    # Compute max index for each score.
    logger.info("Computing maximum indices")
    max_idx = ht.aggregate(hl.struct(**{s: hl.agg.max(ht[s].idx) for s in scores}))

    # Convert index to percentile.
    logger.info("Converting indices to percentiles")
    ht = ht.select(
        **{
            s: ht[s].annotate(percentile=hl.int((ht[s].idx / max_idx[s]) * 100))
            for s in scores
        }
    )

    return ht


def annotate_constraint_data_with_scores(
    constraint_ht: hl.Table,
    percentile_ht: hl.Table,
    additional_constraint_hts: Optional[Dict[str, hl.Table]] = None,
) -> hl.Table:
    """
    Annotate constraint data with missense score percentiles.

    Joins additional constraint tables (by locus, alleles, transcript) and
    score percentiles (by locus, alleles, gene_id) onto the main constraint
    table.

    :param constraint_ht: Main constraint table (preprocessed).
    :param percentile_ht: Table with score percentiles keyed by
        (locus, alleles, gene_id).
    :param additional_constraint_hts: Optional dict mapping names to additional
        constraint tables keyed by (locus, alleles, transcript).
    :return: Annotated constraint table.
    """
    annotations = {}

    # Join additional constraint tables by (locus, alleles, transcript).
    if additional_constraint_hts:
        for name, additional_ht in additional_constraint_hts.items():
            logger.info("Joining additional constraint table: %s", name)
            annotations[name] = additional_ht[
                constraint_ht.locus,
                constraint_ht.alleles,
                constraint_ht.transcript,
            ]

    # Join score percentiles by (locus, alleles, gene_id).
    logger.info("Annotating constraint data with score percentiles")
    percentile_row = percentile_ht[
        constraint_ht.locus, constraint_ht.alleles, constraint_ht.gene_id
    ]
    for field in percentile_ht.row_value:
        annotations[field] = percentile_row[field]

    constraint_ht = constraint_ht.annotate(**annotations)

    return constraint_ht


def aggregate_by_transcript(
    ht: hl.Table,
    scores: List[str],
    include_hi_genes: bool = False,
) -> hl.Table:
    """
    Aggregate observed/expected counts by transcript and score percentile.

    :param ht: Constraint table annotated with score percentiles.
    :param scores: List of score field names.
    :param include_hi_genes: Whether to annotate HI gene categories.
    :return: Aggregated table.
    """
    logger.info("Preparing table for aggregation")

    # Coalesce gene metadata from additional constraint tables where the main
    # table may have missing values.
    coalesce_fields = ["gene", "gene_id", "canonical", "mane_select"]
    additional_tables = ["loeuf_all_rerun_6_12_25_ht", "loeuf_all_3_10_25"]

    coalesce_exprs = {}
    for field in coalesce_fields:
        sources = [ht[field]]
        for tbl in additional_tables:
            if tbl in ht.row:
                sources.append(ht[tbl][field])
        coalesce_exprs[field] = hl.coalesce(*sources)

    ht = ht.annotate(**coalesce_exprs)

    if include_hi_genes:
        hi_genes_dict = get_hi_genes()
        ht = ht.annotate(
            hi_gene=hl.set(hi_genes_dict["all"]).contains(ht.gene),
            hi_gene_category=(
                hl.case()
                .when(
                    hl.set(hi_genes_dict["severe"]).contains(ht.gene),
                    "Severe Haploinsufficient",
                )
                .when(
                    hl.set(hi_genes_dict["moderate"]).contains(ht.gene),
                    "Moderate Haploinsufficient",
                )
                .when(
                    hl.set(hi_genes_dict["new"]).contains(ht.gene),
                    "New Haploinsufficient",
                )
                .default("NA")
            ),
        )

    logger.info("Filtering to valid transcripts")
    ht = ht.filter(hl.is_defined(ht.gene_id) & ht.transcript.startswith("ENST"))

    logger.info("Aggregating by transcript and score percentile")
    agg_ht = ht.group_by(
        "transcript",
        "gene",
        "gene_id",
        "canonical",
        "mane_select",
        "hi_gene",
        "hi_gene_category",
    ).aggregate(
        **{
            s: hl.agg.group_by(
                ht[s].percentile,
                hl.struct(
                    **hl.agg.filter(
                        hl.is_defined(ht[s].percentile)
                        & (ht.calibrate_mu.possible_variants == 1),
                        hl.struct(
                            pos=hl.agg.sum(ht.possible_variants),
                            obs=hl.agg.sum(ht.observed_variants[0]),
                            exp=hl.agg.sum(ht.expected_variants[0]),
                            exp_with_adj_r=hl.agg.sum(
                                ht.expected_variants[0] * hl.or_else(ht.adj_r, 1)
                            ),
                        ),
                    )
                ),
            )
            for s in scores
        }
    )

    logger.info("Filling in missing percentiles")
    agg_ht = agg_ht.annotate(
        **{
            s: [
                agg_ht[s]
                .get(
                    100 - i,
                    hl.struct(pos=0, obs=0, exp=0, exp_with_adj_r=0),
                )
                .annotate(percentile=100 - i, score=s)
                for i in range(0, 101)
            ]
            for s in scores
        }
    )

    return agg_ht


def compute_cumulative_counts(agg_ht: hl.Table, scores: List[str]) -> hl.Table:
    """
    Compute cumulative counts (greater than or equal to percentile).

    :param agg_ht: Aggregated table by transcript and percentile.
    :param scores: List of score field names.
    :return: Table with cumulative counts.
    """
    logger.info("Computing cumulative counts")
    tmp_ht = agg_ht.annotate(
        **{
            s: hl.array_scan(
                lambda i, j: j.annotate(
                    pos_ge_percentile=hl.or_else(i.pos_ge_percentile, 0)
                    + hl.or_else(j.pos, 0),
                    obs_ge_percentile=hl.or_else(i.obs_ge_percentile, 0)
                    + hl.or_else(j.obs, 0),
                    exp_ge_percentile=hl.or_else(i.exp_ge_percentile, 0)
                    + hl.or_else(j.exp, 0),
                    exp_with_adj_r_ge_percentile=hl.or_else(
                        i.exp_with_adj_r_ge_percentile, 0
                    )
                    + hl.or_else(j.exp_with_adj_r, 0),
                ),
                agg_ht[s][0].annotate(
                    pos_ge_percentile=hl.or_else(agg_ht[s][0].pos, 0),
                    obs_ge_percentile=hl.or_else(agg_ht[s][0].obs, 0),
                    exp_ge_percentile=hl.or_else(agg_ht[s][0].exp, 0),
                    exp_with_adj_r_ge_percentile=hl.or_else(
                        agg_ht[s][0].exp_with_adj_r, 0
                    ),
                ),
                agg_ht[s][1:],
            )
            for s in scores
        }
    )

    logger.info("Flattening score arrays")
    tmp_ht = tmp_ht.select(scores=hl.flatten([tmp_ht[s] for s in scores]))
    tmp_ht = tmp_ht.explode("scores")
    tmp_ht = tmp_ht.annotate(**tmp_ht.scores)

    logger.info("Keying by transcript and score information")
    tmp_ht = tmp_ht.key_by(
        "transcript",
        "gene",
        "gene_id",
        "canonical",
        "mane_select",
        "hi_gene",
        "hi_gene_category",
        "score",
        "percentile",
    )

    tmp_ht = tmp_ht.select(
        "pos",
        "obs",
        "exp",
        "pos_ge_percentile",
        "obs_ge_percentile",
        "exp_ge_percentile",
        "exp_with_adj_r_ge_percentile",
    )

    return tmp_ht


def export_percentile_summary(
    ht: hl.Table,
    scores: List[str],
    output_path: str,
    filter_hi_genes: bool = False,
    filter_mane_select: bool = True,
) -> None:
    """
    Export summary of observed/expected counts by percentile.

    :param ht: Constraint table annotated with score percentiles.
    :param scores: List of score field names.
    :param output_path: Path to save summary TSV.
    :param filter_hi_genes: Whether to filter to HI genes only.
    :param filter_mane_select: Whether to filter to MANE select transcripts.
    """
    logger.info("Preparing data for export")

    # Resolve transcript field name (may be "transcript" or "transcript_id").
    transcript_field = "transcript" if "transcript" in ht.row else "transcript_id"

    # Filter scores to only those present in the table.
    available_scores = [s for s in scores if s in ht.row]
    if len(available_scores) < len(scores):
        missing = set(scores) - set(available_scores)
        logger.warning("Scores not found in table (skipping): %s", missing)
    scores = available_scores

    if filter_mane_select:
        ht = ht.filter(ht.mane_select)
        logger.info("Filtered to MANE select transcripts")

    if filter_hi_genes:
        hi_genes_dict = get_hi_genes()
        ht = ht.filter(hl.set(hi_genes_dict["all"]).contains(ht.gene))
        logger.info("Filtered to HI genes")

    ht = ht.filter(ht[transcript_field].startswith("ENST"))

    logger.info("Aggregating by percentile")
    percentiles = ht.aggregate(
        {
            s: hl.agg.group_by(
                ht[s].percentile,
                hl.struct(
                    **hl.agg.filter(
                        hl.is_defined(ht[s].percentile)
                        & (ht.calibrate_mu.possible_variants == 1),
                        hl.struct(
                            pos=hl.agg.sum(ht.possible_variants),
                            obs=hl.agg.sum(ht.observed_variants[0]),
                            exp=hl.agg.sum(ht.expected_variants[0]),
                            exp_with_adj_r=hl.agg.sum(
                                ht.expected_variants[0] * hl.or_else(ht.adj_r, 1)
                            ),
                        ),
                    )
                ),
            )
            for s in scores
        }
    )

    # Match notebook: access percentile keys 0-99, labeled 1-100.
    logger.info("Creating summary table")
    summary_ht = hl.Table.parallelize(
        [
            {
                "percentile": i + 1,
                **{s: percentiles[s][i] for s in scores},
            }
            for i in range(0, 100)
        ]
    )

    logger.info("Exporting to %s", output_path)
    summary_ht.flatten().export(output_path, delimiter="\t")


def compute_plof_per_transcript(per_snv_ht: hl.Table) -> hl.Table:
    """
    Compute adj_r-corrected gene-level pLoF obs/exp from the per-variant table.

    Filters to LOFTEE high-confidence LoF variants (``modifier == "HC"``),
    then aggregates observed and expected counts per transcript with the
    regional depletion correction (``adj_r``).

    :param per_snv_ht: Per-variant expected table (all annotations, not
        pre-filtered to missense). Must have fields: ``modifier``,
        ``observed_variants``, ``expected_variants``, ``adj_r``,
        ``calibrate_mu``, ``gene``, ``transcript``, ``canonical``.
    :return: Table keyed by ``(gene, transcript, canonical)`` with a ``lof``
        struct containing ``obs`` and ``exp`` fields.
    """
    logger.info("Computing adj_r-corrected pLoF obs/exp from per-SNV table")
    logger.info("Per-SNV table key: %s", list(per_snv_ht.key))
    logger.info("Per-SNV table row fields: %s", list(per_snv_ht.row))

    # Flatten the calibrate_mu struct if present (adds gene, transcript,
    # canonical, modifier, observed_variants, expected_variants,
    # possible_variants, etc. as top-level fields).
    if "calibrate_mu" in per_snv_ht.row and "possible_variants" not in per_snv_ht.row:
        per_snv_ht = per_snv_ht.annotate(**per_snv_ht.calibrate_mu)
        logger.info(
            "Flattened calibrate_mu. Row fields now: %s", list(per_snv_ht.row)
        )

    # Filter to LOFTEE HC LoF variants with possible_variants == 1
    # (matching compute_constraint_metrics).
    lof_ht = per_snv_ht.filter(
        (per_snv_ht.modifier == "HC")
        & (per_snv_ht.possible_variants == 1)
    )

    # Coalesce after aggressive filter to reduce empty partitions and
    # avoid shuffle skew during the downstream group_by.
    lof_ht = lof_ht.naive_coalesce(200)

    # Aggregate per transcript with adj_r correction.
    plof_ht = lof_ht.group_by(
        lof_ht.gene,
        lof_ht.transcript,
        lof_ht.canonical,
    ).aggregate(
        lof=hl.struct(
            obs=hl.agg.sum(lof_ht.observed_variants[0]),
            exp=hl.agg.sum(lof_ht.expected_variants[0] * hl.or_else(lof_ht.adj_r, 1)),
        )
    )
    logger.info("pLoF per-transcript aggregation complete")
    return plof_ht


def annotate_plof_data(ht: hl.Table, plof_ht: hl.Table) -> hl.Table:
    """
    Annotate per-variant constraint table with gene-level pLoF counts.

    Joins gene-level pLoF observed and expected counts onto the per-variant
    constraint table. Handles field name differences between tables
    (e.g. ``transcript`` vs ``transcript_id``).

    :param ht: Annotated constraint table with gene/transcript/canonical
        fields (may use ``_id`` suffixed names).
    :param plof_ht: Table keyed by ``(gene, transcript, canonical)`` with a
        ``lof`` struct containing ``obs`` and ``exp`` fields.
    :return: Table with ``plof_obs`` and ``plof_exp`` annotations added.
    """
    logger.info("Annotating with gene-level pLoF observed and expected counts")

    # Log schemas for debugging.
    logger.info(
        "Annotated constraint table row fields: %s", list(ht.row)
    )
    logger.info(
        "Annotated constraint table key: %s", list(ht.key)
    )
    logger.info(
        "pLoF table key: %s, count: available", list(plof_ht.key)
    )

    # Resolve field names: the annotated constraint table may use
    # "transcript_id" / "gene_id" while the pLoF table uses
    # "transcript" / "gene" (or vice versa).
    def _resolve(table, *candidates):
        for c in candidates:
            if c in table.row:
                return table[c]
        raise ValueError(
            f"None of {candidates} found in table fields: {list(table.row)}"
        )

    # Get gene/transcript/canonical from the annotated constraint table.
    ht_gene = _resolve(ht, "gene", "gene_id")
    ht_transcript = _resolve(ht, "transcript", "transcript_id")
    ht_canonical = _resolve(ht, "canonical")

    plof_row = plof_ht[ht_gene, ht_transcript, ht_canonical]
    ht = ht.annotate(
        plof_obs=plof_row.lof.obs,
        plof_exp=plof_row.lof.exp,
    )
    return ht


def export_matched_plof_summary(
    ht: hl.Table,
    scores: List[str],
    output_path: str,
    filter_hi_genes: bool = False,
    filter_mane_select: bool = True,
) -> None:
    """
    Export matched pLoF o/e summary by missense score percentile bin.

    For each score and percentile bin, computes the weighted average pLoF
    observed and expected counts across all genes contributing variants to
    that bin. Both weight types (expected missense count and possible missense
    count) are included in the output for comparison.

    The table must have ``plof_obs`` and ``plof_exp`` annotations (from
    :func:`annotate_plof_data`).

    Observed missense counts are intentionally excluded as weights because
    that would make the pLoF baseline depend on missense depletion itself,
    creating a circular comparison.

    :param ht: Annotated constraint table with score percentiles and
        plof_obs/plof_exp annotations.
    :param scores: List of score field names.
    :param output_path: Path to save summary TSV.
    :param filter_hi_genes: Whether to filter to HI genes only.
    :param filter_mane_select: Whether to filter to MANE select transcripts.
    """
    logger.info("Preparing data for matched pLoF export")

    # Resolve transcript field name (may be "transcript" or "transcript_id").
    transcript_field = "transcript" if "transcript" in ht.row else "transcript_id"

    # Filter scores to only those present in the table.
    available_scores = [s for s in scores if s in ht.row]
    if len(available_scores) < len(scores):
        missing = set(scores) - set(available_scores)
        logger.warning("Scores not found in table (skipping): %s", missing)
    scores = available_scores

    if filter_mane_select:
        ht = ht.filter(ht.mane_select)
        logger.info("Filtered to MANE select transcripts")

    if filter_hi_genes:
        hi_genes_dict = get_hi_genes()
        ht = ht.filter(hl.set(hi_genes_dict["all"]).contains(ht.gene))
        logger.info("Filtered to HI genes")

    ht = ht.filter(ht[transcript_field].startswith("ENST"))

    # Pre-compute weight and pLoF expressions to keep the aggregation readable.
    # Use adjusted expected (expected_variants[0] * adj_r) as the "expected
    # missense" weight so that it is consistent with the missense o/e denominator
    # (exp_with_adj_r).  We must NOT weight by observed missense — that would
    # make the pLoF baseline depend on missense depletion itself (circular).
    exp_adj_weight = ht.expected_variants[0] * hl.or_else(ht.adj_r, 1)
    pos_weight = ht.possible_variants
    plof_obs = hl.or_else(ht.plof_obs, 0)
    plof_exp = hl.or_else(ht.plof_exp, 0)

    logger.info("Aggregating by percentile with matched pLoF metrics")
    percentiles = ht.aggregate(
        {
            s: hl.agg.group_by(
                ht[s].percentile,
                hl.struct(
                    **hl.agg.filter(
                        hl.is_defined(ht[s].percentile)
                        & (ht.calibrate_mu.possible_variants == 1),
                        hl.struct(
                            # Missense metrics (same as export_percentile_summary).
                            pos=hl.agg.sum(ht.possible_variants),
                            obs=hl.agg.sum(ht.observed_variants[0]),
                            exp=hl.agg.sum(ht.expected_variants[0]),
                            exp_with_adj_r=hl.agg.sum(exp_adj_weight),
                            # pLoF weighted by adjusted expected missense count.
                            plof_obs_x_exp=hl.agg.sum(exp_adj_weight * plof_obs),
                            plof_exp_x_exp=hl.agg.sum(exp_adj_weight * plof_exp),
                            # pLoF weighted by possible missense count.
                            plof_obs_x_pos=hl.agg.sum(pos_weight * plof_obs),
                            plof_exp_x_pos=hl.agg.sum(pos_weight * plof_exp),
                            # Transparency: how many variants had pLoF data.
                            n_variants_with_plof=hl.agg.count_where(
                                hl.is_defined(ht.plof_obs)
                            ),
                        ),
                    )
                ),
            )
            for s in scores
        }
    )

    # Match notebook: access percentile keys 0-99, labeled 1-100.
    logger.info("Creating summary table")
    summary_ht = hl.Table.parallelize(
        [
            {
                "percentile": i + 1,
                **{s: percentiles[s][i] for s in scores},
            }
            for i in range(0, 100)
        ]
    )

    logger.info("Exporting matched pLoF summary to %s", output_path)
    summary_ht.flatten().export(output_path, delimiter="\t")


def export_per_gene_percentile_summary(
    ht: hl.Table,
    scores: List[str],
    output_path: str,
    filter_hi_genes: bool = False,
    filter_mane_select: bool = True,
) -> None:
    """
    Export per-gene missense and pLoF stats by score percentile bin.

    For each gene and percentile bin, exports the gene's missense observed,
    expected, and possible counts alongside its gene-level pLoF observed and
    expected counts. This allows inspection of individual gene contributions
    to the matched pLoF curve.

    The table must have ``plof_obs`` and ``plof_exp`` annotations (from
    :func:`annotate_plof_data`).

    :param ht: Annotated constraint table with score percentiles and
        plof_obs/plof_exp annotations.
    :param scores: List of score field names.
    :param output_path: Path to save per-gene summary TSV.
    :param filter_hi_genes: Whether to filter to HI genes only.
    :param filter_mane_select: Whether to filter to MANE select transcripts.
    """
    logger.info("Preparing data for per-gene percentile export")

    # Resolve transcript field name (may be "transcript" or "transcript_id").
    transcript_field = "transcript" if "transcript" in ht.row else "transcript_id"

    # Filter scores to only those present in the table.
    available_scores = [s for s in scores if s in ht.row]
    if len(available_scores) < len(scores):
        missing = set(scores) - set(available_scores)
        logger.warning("Scores not found in table (skipping): %s", missing)
    scores = available_scores

    if filter_mane_select:
        ht = ht.filter(ht.mane_select)
        logger.info("Filtered to MANE select transcripts")

    if filter_hi_genes:
        hi_genes_dict = get_hi_genes()
        ht = ht.filter(hl.set(hi_genes_dict["all"]).contains(ht.gene))
        logger.info("Filtered to HI genes")

    ht = ht.filter(ht[transcript_field].startswith("ENST"))

    # Pre-compute commonly used expressions.
    exp_adj_weight = ht.expected_variants[0] * hl.or_else(ht.adj_r, 1)

    logger.info("Aggregating per gene and percentile")
    per_gene = ht.group_by(
        ht.gene,
        ht[transcript_field],
        ht.canonical,
        ht.mane_select,
    ).aggregate(
        plof_obs=hl.agg.min(hl.or_else(ht.plof_obs, 0)),
        plof_exp=hl.agg.min(hl.or_else(ht.plof_exp, 0.0)),
        **{
            s: hl.agg.group_by(
                ht[s].percentile,
                hl.struct(
                    **hl.agg.filter(
                        hl.is_defined(ht[s].percentile)
                        & (ht.calibrate_mu.possible_variants == 1),
                        hl.struct(
                            pos=hl.agg.sum(ht.possible_variants),
                            obs=hl.agg.sum(ht.observed_variants[0]),
                            exp=hl.agg.sum(ht.expected_variants[0]),
                            exp_with_adj_r=hl.agg.sum(exp_adj_weight),
                        ),
                    )
                ),
            )
            for s in scores
        },
    )

    # Convert per-score dicts to arrays of 100 entries (percentiles 0-99,
    # labeled 1-100 to match export_percentile_summary convention).
    logger.info("Converting score dicts to arrays")
    default = hl.struct(pos=0, obs=0, exp=0, exp_with_adj_r=0)
    per_gene = per_gene.annotate(
        **{
            s: hl.range(100).map(lambda i: per_gene[s].get(i, default))
            for s in scores
        }
    )

    # Create one row per (gene x percentile) by zipping score arrays.
    logger.info("Exploding to one row per gene x percentile")
    per_gene = per_gene.annotate(
        _rows=hl.range(100).map(
            lambda i: hl.struct(
                percentile=i + 1,
                **{s: per_gene[s][i] for s in scores},
            )
        )
    )
    # Drop original score array fields to avoid name conflicts with the
    # exploded struct fields.
    per_gene = per_gene.select("plof_obs", "plof_exp", "_rows")
    per_gene = per_gene.explode("_rows")
    per_gene = per_gene.transmute(**per_gene._rows)

    logger.info("Exporting per-gene summary to %s", output_path)
    per_gene.flatten().export(output_path, delimiter="\t")


def main(args):
    """Execute the missense score percentile computation pipeline."""
    logger.info("Starting missense score percentile computation pipeline")

    # Initialize Hail.
    hl.init(log=args.log_path, tmp_dir=args.tmp_dir)

    if args.use_new_shuffle:
        logger.info("Enabling new shuffle")
        hl._set_flags(use_new_shuffle="1")

    scores = args.scores or DEFAULT_SCORES

    # Step 1: Preprocess missense scores.
    if args.preprocess_scores:
        logger.info("Step 1: Preprocessing missense scores")
        logger.info("Loading missense scores from %s", args.input_scores_ht_path)
        ht = hl.read_table(args.input_scores_ht_path)
        ht = preprocess_missense_scores(
            ht,
            scores,
            gene_id_field=args.gene_id_field,
            transcript_id_field=args.transcript_id_field,
        )
        if args.preprocessed_scores_ht_path:
            logger.info(
                "Checkpointing preprocessed table to %s",
                args.preprocessed_scores_ht_path,
            )
            ht = ht.checkpoint(
                args.preprocessed_scores_ht_path, overwrite=args.overwrite
            )

    # Step 2: Compute percentiles.
    if args.compute_percentiles:
        logger.info("Step 2: Computing score percentiles")
        if not args.preprocess_scores:
            logger.info(
                "Loading preprocessed scores from %s",
                args.preprocessed_scores_ht_path,
            )
            ht = hl.read_table(args.preprocessed_scores_ht_path)
        ht = compute_score_percentiles(
            ht,
            scores,
            checkpoint_prefix=args.preprocessed_scores_ht_path,
            overwrite=args.overwrite,
        )
        if args.percentile_ht_path:
            logger.info(
                "Checkpointing percentile table to %s", args.percentile_ht_path
            )
            ht = ht.checkpoint(args.percentile_ht_path, overwrite=args.overwrite)

    # Step 3: Annotate constraint data.
    if args.annotate_constraint:
        logger.info("Step 3: Annotating constraint data with percentiles")

        # Load and preprocess the main constraint table.
        logger.info("Loading constraint data from %s", args.constraint_ht_path)
        constraint_ht = hl.read_table(args.constraint_ht_path)
        constraint_ht = preprocess_constraint_table(constraint_ht)

        # Load percentile table.
        logger.info("Loading percentile data from %s", args.percentile_ht_path)
        percentile_ht = hl.read_table(args.percentile_ht_path)

        # Load and preprocess additional constraint tables.
        additional_hts = None
        if args.additional_constraint_ht_path_1 or args.additional_constraint_ht_path_2:
            additional_hts = {}
            if args.additional_constraint_ht_path_1:
                logger.info(
                    "Loading additional constraint table 1 from %s",
                    args.additional_constraint_ht_path_1,
                )
                add_ht1 = hl.read_table(args.additional_constraint_ht_path_1)
                add_ht1 = preprocess_constraint_table(add_ht1)
                groupings1 = [
                    g
                    for g in hl.eval(add_ht1.apply_models_globals.groupings)
                    if g not in list(MU_GROUPING) + ["locus", "alleles", "transcript"]
                ]
                additional_hts["loeuf_all_rerun_6_12_25_ht"] = (
                    add_ht1.key_by("locus", "alleles", "transcript")
                    .select(
                        *groupings1,
                        *MU_GROUPING,
                        "mu_snp",
                        "mu",
                        "observed_variants",
                        "possible_variants",
                        "predicted_proportion_observed",
                        "coverage_correction",
                        "expected_variants",
                        "calibrate_mu",
                        "adj_r",
                        "sfs_bin",
                    )
                )
                if args.additional_constraint_checkpoint_1:
                    additional_hts["loeuf_all_rerun_6_12_25_ht"] = additional_hts[
                        "loeuf_all_rerun_6_12_25_ht"
                    ].checkpoint(
                        args.additional_constraint_checkpoint_1,
                        overwrite=args.overwrite,
                    )

            if args.additional_constraint_ht_path_2:
                logger.info(
                    "Loading additional constraint table 2 from %s",
                    args.additional_constraint_ht_path_2,
                )
                add_ht2 = hl.read_table(args.additional_constraint_ht_path_2)
                add_ht2 = preprocess_constraint_table(add_ht2)
                groupings2 = [
                    g
                    for g in hl.eval(add_ht2.apply_models_globals.groupings)
                    if g not in list(MU_GROUPING) + ["locus", "alleles", "transcript"]
                ]
                additional_hts["loeuf_all_3_10_25"] = (
                    add_ht2.key_by("locus", "alleles", "transcript")
                    .select(
                        *groupings2,
                        *MU_GROUPING,
                        "observed_variants",
                        "possible_variants",
                        "expected_variants",
                        "adj_r",
                    )
                )
                if args.additional_constraint_checkpoint_2:
                    additional_hts["loeuf_all_3_10_25"] = additional_hts[
                        "loeuf_all_3_10_25"
                    ].checkpoint(
                        args.additional_constraint_checkpoint_2,
                        overwrite=args.overwrite,
                    )

        # Annotate.
        constraint_ht = annotate_constraint_data_with_scores(
            constraint_ht,
            percentile_ht,
            additional_constraint_hts=additional_hts,
        )
        if args.annotated_constraint_ht_path:
            logger.info(
                "Writing annotated table to %s", args.annotated_constraint_ht_path
            )
            constraint_ht = constraint_ht.checkpoint(
                args.annotated_constraint_ht_path, overwrite=args.overwrite
            )

    # Step 4: Export percentile summaries.
    if args.export_summaries:
        logger.info("Step 4: Exporting percentile summaries")
        if not args.annotate_constraint:
            logger.info(
                "Loading annotated constraint table from %s",
                args.annotated_constraint_ht_path,
            )
            constraint_ht = hl.read_table(
                args.annotated_constraint_ht_path
            ).naive_coalesce(1000)

        if args.export_hi_genes:
            export_percentile_summary(
                constraint_ht,
                scores,
                args.hi_genes_summary_path,
                filter_hi_genes=True,
                filter_mane_select=True,
            )

        if args.export_all_genes:
            export_percentile_summary(
                constraint_ht,
                scores,
                args.all_genes_summary_path,
                filter_hi_genes=False,
                filter_mane_select=True,
            )

    # Step 5: Aggregate by transcript.
    if args.aggregate_by_transcript:
        logger.info("Step 5: Aggregating by transcript")
        if not args.annotate_constraint:
            logger.info(
                "Loading annotated constraint table from %s",
                args.annotated_constraint_ht_path,
            )
            constraint_ht = hl.read_table(
                args.annotated_constraint_ht_path
            ).naive_coalesce(1000)

        agg_ht = aggregate_by_transcript(
            constraint_ht,
            scores,
            include_hi_genes=args.include_hi_genes,
        )
        if args.aggregated_transcript_ht_path:
            logger.info(
                "Checkpointing aggregated table to %s",
                args.aggregated_transcript_ht_path,
            )
            agg_ht = agg_ht.checkpoint(
                args.aggregated_transcript_ht_path, overwrite=args.overwrite
            )

    # Step 6: Compute cumulative counts.
    if args.compute_cumulative:
        logger.info("Step 6: Computing cumulative counts")
        if not args.aggregate_by_transcript:
            logger.info(
                "Loading aggregated transcript table from %s",
                args.aggregated_transcript_ht_path,
            )
            agg_ht = hl.read_table(args.aggregated_transcript_ht_path)
        cumulative_ht = compute_cumulative_counts(agg_ht, scores)
        if args.cumulative_ht_path:
            logger.info(
                "Writing cumulative counts table to %s", args.cumulative_ht_path
            )
            cumulative_ht = cumulative_ht.checkpoint(
                args.cumulative_ht_path, overwrite=args.overwrite
            )
        if args.cumulative_tsv_path:
            logger.info(
                "Exporting cumulative counts TSV to %s", args.cumulative_tsv_path
            )
            cumulative_ht.flatten().export(args.cumulative_tsv_path)

    # Step 7: Export matched pLoF summary.
    if args.export_matched_plof_summary:
        logger.info("Step 7: Exporting matched pLoF o/e summary")

        # Load the annotated constraint table if not already in memory.
        if not args.annotate_constraint:
            logger.info(
                "Loading annotated constraint table from %s",
                args.annotated_constraint_ht_path,
            )
            constraint_ht = hl.read_table(
                args.annotated_constraint_ht_path
            ).naive_coalesce(1000)

        # Compute adj_r-corrected gene-level pLoF from the per-SNV table
        # (LOFTEE HC LoF, possible_variants == 1).
        logger.info(
            "Loading per-SNV table for pLoF computation from %s",
            args.constraint_ht_path,
        )
        per_snv_ht = hl.read_table(args.constraint_ht_path)
        plof_ht = compute_plof_per_transcript(per_snv_ht)
        constraint_ht = annotate_plof_data(constraint_ht, plof_ht)

        # Export all genes matched pLoF summary.
        if args.matched_plof_all_genes_path:
            export_matched_plof_summary(
                constraint_ht,
                scores,
                args.matched_plof_all_genes_path,
                filter_hi_genes=False,
                filter_mane_select=True,
            )

        # Export HI genes matched pLoF summary.
        if args.matched_plof_hi_genes_path:
            export_matched_plof_summary(
                constraint_ht,
                scores,
                args.matched_plof_hi_genes_path,
                filter_hi_genes=True,
                filter_mane_select=True,
            )

        # Export per-gene percentile summary.
        if args.per_gene_all_genes_path:
            export_per_gene_percentile_summary(
                constraint_ht,
                scores,
                args.per_gene_all_genes_path,
                filter_hi_genes=False,
                filter_mane_select=True,
            )
        if args.per_gene_hi_genes_path:
            export_per_gene_percentile_summary(
                constraint_ht,
                scores,
                args.per_gene_hi_genes_path,
                filter_hi_genes=True,
                filter_mane_select=True,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compute percentiles for missense prediction scores and annotate"
            " constraint data"
        )
    )

    # Input paths.
    parser.add_argument(
        "--input-scores-ht-path",
        help="Path to input Genetics Gym missense scores table",
        type=str,
        default="gs://trisha-tmp/new_VSM_temp/vsm_all_tables/vsm_all_SNP_gene_lag.ht",
    )
    parser.add_argument(
        "--constraint-ht-path",
        help="Path to constraint table",
        type=str,
        default=(
            "gs://gnomad/v4.1/constraint_coverage_corrected/apply_models"
            "/transcript_consequences/gnomad.v4.1.per_variant_expected"
            ".coverage_corrected.with_downsamplings.ht"
        ),
    )
    parser.add_argument(
        "--additional-constraint-ht-path-1",
        help="Path to first additional constraint table (loeuf_all_rerun_6_12_25)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--additional-constraint-ht-path-2",
        help="Path to second additional constraint table (loeuf_all_3_10_25)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--additional-constraint-checkpoint-1",
        help="Checkpoint path for preprocessed additional constraint table 1",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--additional-constraint-checkpoint-2",
        help="Checkpoint path for preprocessed additional constraint table 2",
        type=str,
        default=None,
    )

    # Intermediate/output paths.
    parser.add_argument(
        "--preprocessed-scores-ht-path",
        help="Path to save preprocessed scores table",
        type=str,
        default="gs://gnomad-tmp-4day/persist_Tablefm0RiyrmoQ",
    )
    parser.add_argument(
        "--percentile-ht-path",
        help="Path to save percentile table",
        type=str,
        default="gs://gnomad-tmp-4day/julia/all_missense_scores_percentile.ht",
    )
    parser.add_argument(
        "--annotated-constraint-ht-path",
        help="Path to save annotated constraint table",
        type=str,
        default=(
            "gs://gnomad/v4.1/constraint/resources"
            "/all_missense_scores_percentile.annotate_with_oe.12_23_25.ht"
        ),
    )
    parser.add_argument(
        "--aggregated-transcript-ht-path",
        help="Path to save aggregated transcript table",
        type=str,
        default="gs://gnomad-tmp-4day/persist_TableGZFTL1BCP0",
    )
    parser.add_argument(
        "--cumulative-ht-path",
        help="Path to save cumulative counts table",
        type=str,
        default=(
            "gs://gnomad/v4.1/constraint/resources/from_julia_goodrich"
            "/all_missense_scores_percentile.annotate_with_gnomad_v4.1_oe"
            ".aggregated_by_transcript.ht"
        ),
    )
    parser.add_argument(
        "--cumulative-tsv-path",
        help="Path to save cumulative counts TSV",
        type=str,
        default=(
            "gs://gnomad/v4.1/constraint/resources/from_julia_goodrich"
            "/all_missense_scores_percentile.annotate_with_gnomad_v4.1_oe"
            ".aggregated_by_transcript.tsv.gz"
        ),
    )
    parser.add_argument(
        "--all-genes-summary-path",
        help="Path to save all genes percentile summary",
        type=str,
        default=(
            "gs://gnomad-julia/loeuf_all"
            "/all_missense_scores_percentile.all_genes.mane_select.oe.tsv"
        ),
    )
    parser.add_argument(
        "--hi-genes-summary-path",
        help="Path to save HI genes percentile summary",
        type=str,
        default=(
            "gs://gnomad-julia/loeuf_all"
            "/all_missense_scores_percentile.hi_genes.oe.tsv"
        ),
    )

    # Matched pLoF output paths.
    parser.add_argument(
        "--constraint-metrics-ht-path",
        help=(
            "DEPRECATED: No longer used. Gene-level pLoF is now computed"
            " from the per-SNV table (--constraint-ht-path) with adj_r."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--matched-plof-all-genes-path",
        help="Path to save all genes matched pLoF summary TSV",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--matched-plof-hi-genes-path",
        help="Path to save HI genes matched pLoF summary TSV",
        type=str,
        default=None,
    )

    # Per-gene export paths.
    parser.add_argument(
        "--per-gene-all-genes-path",
        help="Path to save per-gene percentile summary TSV (all genes)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--per-gene-hi-genes-path",
        help="Path to save per-gene percentile summary TSV (HI genes only)",
        type=str,
        default=None,
    )

    # Field name options.
    parser.add_argument(
        "--gene-id-field",
        help="Field name for gene ID in the input scores table (e.g., 'ensg')",
        type=str,
        default="ensg",
    )
    parser.add_argument(
        "--transcript-id-field",
        help="Field name for transcript ID in the input scores table",
        type=str,
        default=None,
    )

    # General options.
    parser.add_argument(
        "--scores",
        help="List of score field names to process",
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "--tmp-dir",
        help="Temporary directory for Hail",
        type=str,
        default="gs://gnomad-tmp-4day",
    )
    parser.add_argument(
        "--log-path",
        help="Path for Hail log file",
        type=str,
        default="/determine_missense_score_percentiles.log",
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing output files",
        action="store_true",
    )
    parser.add_argument(
        "--include-hi-genes",
        help="Whether to include HI gene annotations",
        action="store_true",
    )
    parser.add_argument(
        "--use-new-shuffle",
        help="Enable Hail's new shuffle implementation",
        action="store_true",
    )

    # Pipeline steps.
    parser.add_argument(
        "--preprocess-scores",
        help="Preprocess missense scores",
        action="store_true",
    )
    parser.add_argument(
        "--compute-percentiles",
        help="Compute score percentiles",
        action="store_true",
    )
    parser.add_argument(
        "--annotate-constraint",
        help="Annotate constraint data with percentiles",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-by-transcript",
        help="Aggregate by transcript",
        action="store_true",
    )
    parser.add_argument(
        "--compute-cumulative",
        help="Compute cumulative counts",
        action="store_true",
    )
    parser.add_argument(
        "--export-summaries",
        help="Export percentile summaries",
        action="store_true",
    )
    parser.add_argument(
        "--export-all-genes",
        help="Export all genes summary",
        action="store_true",
    )
    parser.add_argument(
        "--export-hi-genes",
        help="Export HI genes summary",
        action="store_true",
    )
    parser.add_argument(
        "--export-matched-plof-summary",
        help="Export matched pLoF o/e summary by missense score percentile",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
