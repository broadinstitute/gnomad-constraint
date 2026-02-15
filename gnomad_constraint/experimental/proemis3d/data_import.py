"""Script to import data for Proemis3D.

This script imports data for Proemis3D, including COSMIS scores, Varity data, MTR3D data,
InterPro annotations, Kaplanis variants, Fu variants, ClinVar missense variants,
constraint metrics, MTR data, RMC data, context data, and Genetics Gym missense scores.


Information about the AlphaFold DB:

As per the article, https://academic.oup.com/nar/article/52/D1/D368/7337620:

"Data archiving for AlphaFold DB began with an initial release in July 2021, housing
over 360 000 structures for 20 model organism proteomes with sequences derived from
the ‘one sequence per gene’ reference proteomes provided in UniProt release 2021_02.
In December 2021, most of the reviewed sequences in UniProt, i.e. the Swiss-Prot
dataset, were incorporated from the UniProt release 2021_04. In January 2022,
proteomes relevant to global health, derived from priority lists by the World Health
Organization, were added, utilising sequences from UniProt release 2021_04 ‘one
sequence per gene’ reference proteomes. By July 2022, most of the remaining sequences
from UniProt release 2021_04 were included, featuring an additional TAR file on the
AFDB download page, EMBL-EBI’s FTP and Google Cloud Datasets, containing predictions
in MANE select (20).

A November 2022 update rectified structures affected by a temporary numerical bug
presented in the July release. This bug led to low accuracy predictions with
correspondingly low pLDDT for ∼4% of the total structure predictions in the database.
As part of this update, the coordinates for affected structures were updated (old
coordinate files remain accessible as v3 files), and minor metadata adjustments were
made in the mmCIF files for the remaining structures. We document every data version
update in our changelog at https://ftp.ebi.ac.uk/pub/databases/alphafold/CHANGELOG.txt.

As of September 2023, the EMBL-EBI’s FTP area hosts TAR files for proteomes of 48
organisms, including model organisms and WHO pathogens of interest (Supplementary
Table S1). The complete dataset is stored on Google Cloud Platform (GCP) and is
accessed through a file access API. The metadata is indexed using Apache-Solr powering
search API to facilitate data accessibility and searchability.

The database provides access to over 214 million predicted structures, although some
sequences might be outdated compared to UniProt due to less frequent data releases in
the AlphaFold DB. Predictions of UniProt sequences are outputs of a single model run.
In contrast, Swiss-Prot/proteomes entries represent the most confident prediction from
runs of five models trained with different random seeds. The following sequences are
not covered in the database: (i) those that are less than 16 amino acids, or (ii) >2700
for SwissProt or proteome sequences and 1280 for other UniProt sequences, or (iii) those
that contain non-standard amino acids, or (iv) are not in the UniProt ‘one sequence
per gene’ FASTA file, or (v) viral proteins. These limitations are under discussion.

The full dataset, housing all predictions, is accessible from Google Cloud Public
Datasets under a CC-BY-4.0 license. This dataset, approximately 23 TiB in size, is
available at the following Google Cloud Storage Bucket:
gs://public-datasets-deepmind-alphafold-v4. We suggest that most users download only
the subset of files relevant to their specific use case to optimise resources.
However, if a complete dataset is required for local processing, as might be the case
in an academic high-performance computing centre, it can be downloaded in roughly
2.5 days using a 1 Gbps internet connection. Importantly, a Google account is
necessary for the download."
"""

import argparse

import hail as hl
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.constraint import oe_confidence_interval
from gnomad.utils.liftover import default_lift_data
from hail.utils.misc import divide_null

from gnomad_constraint.experimental.proemis3d.resources import (
    CURRENT_VERSION,
    get_clinvar_missense_ht,
    get_constraint_metrics_ht,
    get_context_preprocessed_ht,
    get_cosmis_score_ht,
    get_cosmis_score_tsv,
    get_fu_variants_ht,
    get_fu_variants_tsv,
    get_genetics_gym_missense_scores_ht,
    get_gnomad_de_novo_ht,
    get_insilico_annotations_ht,
    get_interpro_annotations,
    get_interpro_annotations_ht,
    get_kaplanis_sig_variants_tsv,
    get_kaplanis_variants_ht,
    get_kaplanis_variants_tsv,
    get_mtr3d_ht,
    get_mtr3d_tsv,
    get_mtr_ht,
    get_mtr_tsv,
    get_processed_genetics_gym_missense_scores_ht,
    get_revel_csv,
    get_rmc_ht,
    get_temp_context_preprocessed_ht,
    get_temp_processed_constraint_ht,
    get_temp_processed_rmc_ht,
    get_varity_ht,
    get_varity_tsv,
)


def import_cosmis_score_data(
    model: str,
) -> hl.Table:
    """
    Import and process COSMIS TSV file for a specified structure model.

    :param model: Structure model source to process ('alphafold', 'swiss_model', or 'pdb').
    :param distance: Distance metric to use ('8a' or '10a'). Default is '8a'.
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Hail Table with COSMIS scores.
    """
    ht = hl.import_table(get_cosmis_score_tsv(model), force=True, impute=True)
    ht = ht.transmute(
        transcript_id=ht.enst_id,
        residue_index=ht.uniprot_pos - 1,
        cossyn=hl.float(ht.cossyn),
    )
    ht = ht.key_by("transcript_id", "uniprot_id", "residue_index")

    return ht


def import_varity_data() -> hl.Table:
    """
    Import Varity data.

    :return: Hail Table with Varity data.
    """
    ht = hl.import_table(
        get_varity_tsv(),
        impute=True,
        delimiter="\t",
        min_partitions=1000,
    )
    ht = ht.select(
        uniprot_id=ht.p_vid,
        residue_index=ht.aa_pos - 1,
        residue_ref=ht.aa_ref,
        residue_alt=ht.aa_alt,
        varity_r=ht.VARITY_R,
        varity_er=ht.VARITY_ER,
        varity_r_loo=ht.VARITY_R_LOO,
        varity_er_loo=ht.VARITY_ER_LOO,
    )
    ht = ht.key_by("uniprot_id", "residue_index", "residue_ref", "residue_alt")

    ht.show()

    return ht


def import_mtr3d_data() -> hl.Table:
    """
    Import MTR3D data.

    :return: Hail Table with MTR3D data.
    """
    ht = hl.import_table(
        get_mtr3d_tsv(),
        impute=True,
        delimiter=",",
        min_partitions=1000,
        missing="",
    )
    ht = ht.key_by(
        transcript_id=ht.transcript,
        uniprot_id=ht.uniprot,
        residue_index=ht.aa_num - 1,
    )
    ht = ht.select(
        "mean_pLDDT", "mtr3daf2_5a", "mtr3daf2_8a", "mtr3daf2_11a", "mtr3daf2_14a"
    )

    return ht


def import_mtr_data() -> hl.Table:
    """
    Import MTR data.

    :return: Hail Table with MTR data.
    """
    ht = hl.import_table(
        get_mtr_tsv(),
        impute=True,
        delimiter="\t",
        min_partitions=1000,
        missing="",
    )
    ht = ht.annotate(
        locus=hl.parse_locus(
            hl.format(
                "chr%s:%s", hl.if_else(ht.CHR == 23, "X", hl.str(ht.CHR)), ht.POS
            ),
            reference_genome="GRCh38",
        ),
        alleles=[ht.REF, ht.ALT],
        transcript_id=ht.TranscriptId,
        mtr=ht.MTR,
    )
    ht = ht.key_by("locus", "alleles", "transcript_id")
    ht = ht.select(
        "mtr",
        "synExp",
        "misExp",
        "expMTR",
        "synObs",
        "misObs",
        "obsMTR",
        "adj_rate",
        "pvalue",
        "qvalue",
        "proteinLength",
        "MTRpercentile_exome",
        "MTRpercentile_transcript",
    )

    return ht


def import_kaplanis_variants(
    liftover_to_grch38: bool = False,
    key_by_gene_and_transcript: bool = False,
) -> hl.Table:
    """
    Process Kaplanis de novo missense variants file and lift over loci to GRCh38.

    :param liftover_to_grch38: Whether to lift over loci to GRCh38. Default is False.
    :param key_by_gene_and_transcript: Whether to key the table by gene and transcript.
        Default is False.
    :return: Hail Table with Kaplanis de novo missense variants.
    """
    ht = hl.import_table(get_kaplanis_variants_tsv(), min_partitions=10, impute=True)

    variant_expr = ht.Variant.split(":")
    ht = ht.select(
        locus=hl.locus(
            variant_expr[0], hl.int(variant_expr[1]), reference_genome="GRCh37"
        ),
        transcript_id=ht.Transcript_ID,
        gene_id=ht.Gene_ID,
        gene_name=ht.Gene,
        alleles=[variant_expr[2], variant_expr[3]],
    )

    # Add case_control field to keep track of the number of individuals carrying the
    # variant in the Kaplanis study.
    ht = ht.annotate(case_control="DD")
    ht = ht.group_by("locus", "alleles", "gene_id", "transcript_id").aggregate(
        case_control=hl.agg.collect(ht.case_control),
    )

    # Set default select and key by fields.
    select_fields = ["case_control"]
    key_by_fields = ["locus", "alleles"]
    if liftover_to_grch38:
        # Perform liftover.
        ht = default_lift_data(ht)
        ht = ht.key_by()
        ht = ht.transmute(
            locus=ht.new_locus,
            alleles=ht.new_alleles,
            grch37_locus=ht.original_locus,
            grch37_alleles=ht.original_alleles,
        )
        select_fields = ["grch37_locus", "grch37_alleles"] + select_fields

    if key_by_gene_and_transcript:
        key_by_fields = key_by_fields + ["gene_id", "transcript_id"]

    ht = ht.key_by(*key_by_fields)
    ht = ht.select(*select_fields)

    return ht


def get_kaplanis_sig_gene_annotations(
    gene_name_expr: hl.expr.StringExpression,
) -> hl.Table:
    """
    Get Kaplanis significant gene set annotations.

    :param gene_name_expr: Gene name expression.
    :return: Struct with boolean expressions for whether the gene is in the significant,
        diagnostic consensus, or sig or diagnostic consensus gene set.
    """
    # Load significant gene set
    sig_gene_ht = hl.import_table(get_kaplanis_sig_variants_tsv())
    sig_gene_expr = sig_gene_ht.significant == "TRUE"
    diagnostic_category_consensus_expr = sig_gene_ht.diagnostic_category == "consensus"
    agg_expr = hl.agg.collect_as_set(sig_gene_ht.symbol)
    sig_genes = hl.literal(
        sig_gene_ht.aggregate(
            {
                "sig_gene": hl.agg.filter(sig_gene_expr, agg_expr),
                "consensus_gene": hl.agg.filter(
                    diagnostic_category_consensus_expr, agg_expr
                ),
                "sig_or_consensus_gene": hl.agg.filter(
                    sig_gene_expr | diagnostic_category_consensus_expr, agg_expr
                ),
            }
        )
    )
    return hl.struct(
        in_sig_gene=sig_genes["sig_gene"].contains(gene_name_expr),
        in_diagnostic_consensus_gene=(
            sig_genes["consensus_gene"].contains(gene_name_expr)
        ),
        in_sig_or_diagnostic_consensus_gene=(
            sig_genes["sig_or_consensus_gene"].contains(gene_name_expr)
        ),
    )


def import_fu_variants() -> None:
    """
    Import de novo variants from Fu et al. (2022) paper.

    Function imports variants from TSV into HT.

    :return: Hail Table with Fu de novo variants.
    """
    fu_ht = hl.import_table(
        get_fu_variants_tsv(),
        impute=True,
        # Skip blank lines at the bottom of this TSV
        missing="",
        skip_blank_lines=True,
    )

    # Remove lines from bottom of TSV that are parsed incorrectly upon import
    # These lines contain metadata about the TSV, e.g.:
    # "Supplementary Table 20. The de novo SNV/indel variants used in TADA
    # association analyses from assembled ASD cohorts"
    fu_ht = fu_ht.filter(~hl.is_missing(fu_ht.Role))
    fu_ht = fu_ht.annotate(
        locus=hl.parse_locus(
            hl.format(
                "chr%s:%s",
                fu_ht.Variant.split(":")[0],
                fu_ht.Variant.split(":")[1],
            ),
            reference_genome="GRCh38",
        ),
        alleles=[fu_ht.Variant.split(":")[2], fu_ht.Variant.split(":")[3]],
    )

    # Rename 'Proband' > 'ASD' and 'Sibling' > 'control'
    fu_ht = fu_ht.transmute(role=hl.if_else(fu_ht.Role == "Proband", "ASD", "control"))
    fu_ht = fu_ht.group_by("locus", "alleles").aggregate(
        role=hl.agg.collect(fu_ht.role),
    )

    return fu_ht


def import_interpro_annotations() -> hl.Table:
    """
    Import InterPro annotations.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Hail Table with InterPro annotations.
    """
    ht = hl.import_table(get_interpro_annotations(), impute=True, min_partitions=500)
    ht = ht.select(
        uniprot_id=ht["UniProtKB/Swiss-Prot ID"],
        transcript_id=ht["Transcript stable ID"],
        residue_index=hl.range(
            hl.int(hl.or_missing(ht["Interpro start"] != "", ht["Interpro start"])) - 1,
            hl.int(hl.or_missing(ht["Interpro end"] != "", ht["Interpro end"])),
        ),
        interpro_id=ht["Interpro ID"],
        interpro_short_description=ht["Interpro Short Description"],
        interpro_description=ht["Interpro Description"],
    )
    ht = ht.explode("residue_index")
    ht = ht.key_by("transcript_id", "uniprot_id", "residue_index")

    return ht


def process_clinvar_ht(clinvar_version: str = "20250504") -> hl.Table:
    """
    Process ClinVar HT by extracting info fields, filtering to missense variants, and rekey.

    :param clinvar_version: Version of ClinVar to use. Default is `20250504`.
    :return: Hail Table with ClinVar missense variants.
    """
    ht = clinvar.versions[clinvar_version].ht()
    ht = ht.select("rsid", **ht.info)
    ht = ht.filter(hl.any(ht.MC.map(lambda x: x.split("\\|")[1] == "missense_variant")))
    ht = ht.key_by("locus", "alleles", gene=ht.GENEINFO.split(":")[0])

    return ht


def process_gnomad_site_ht(ht) -> hl.Table:
    """
    Process gnomAD site Hail Table.

    :param ht: Hail Table to process.
    :return: Hail Table with gnomAD site variants.
    """
    return ht.select(gnomad_exomes_flags=ht.exome.flags)


def process_pext_base_ht(ht) -> hl.Table:
    """
    Process PEXT base level Hail Table.

    :param ht: Hail Table to process.
    :return: Hail Table with PEXT base level variants.
    """
    return ht.key_by("locus", "gene_id").drop("gene_symbol")


def process_pext_annotation_ht(ht) -> hl.Table:
    """
    Process PEXT annotation level Hail Table.

    :param ht: Hail Table to process.
    :return: Hail Table with PEXT annotation level variants.
    """
    return ht.key_by("locus", "alleles", "gene_id", "most_severe_consequence").drop(
        "gene_symbol"
    )


def process_gnomad_de_novo_ht(ht) -> hl.Table:
    """
    Process gnomAD de novo Hail Table.

    :return: Hail Table with gnomAD de novo variants.
    """
    ht = ht.key_by("locus", "alleles")
    ht = ht.select("de_novo_AC", "p_de_novo_stats")

    return ht


def process_rmc_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Load RMC Hail Table and annotate with p-values and per-locus expansion.

    Annotates:
       - `section_oe_upper`: Upper bound of the observed/expected confidence interval.
       - `section_p_value`: Chi-square p-value.
       - `locus`: List of loci covered by the interval.

    Explodes the `locus` array into multiple rows and keys the table by `locus` and
    `transcript`.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Hail Table with RMC data.
    """
    rmc_ht = get_rmc_ht(version).ht()

    rmc_ht = rmc_ht.annotate(
        section_oe=divide_null(rmc_ht.section_obs, rmc_ht.section_exp),
        section_oe_ci=oe_confidence_interval(rmc_ht.section_obs, rmc_ht.section_exp),
        section_p_value=hl.pchisqtail(rmc_ht.section_chisq, 1),
        locus=hl.range(
            rmc_ht.interval.start.position, rmc_ht.interval.end.position + 1
        ).map(
            lambda x: hl.locus(
                rmc_ht.interval.start.contig, x, reference_genome="GRCh38"
            )
        ),
    )

    rmc_ht = rmc_ht.explode("locus").key_by("locus", "transcript")

    return rmc_ht


def process_constraint_metrics_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Load constraint metrics Hail Table and filter to transcripts starting with "ENST".

    Only selects the syn, mis, and lof constraint groups and only the first oe_info
    struct, which contains info for the full dataset, not per-genetic ancestry group.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Hail Table with constraint metrics.
    """
    ht = get_constraint_metrics_ht(version).ht()
    ht = ht.filter(ht.transcript.startswith("ENST"))

    ht = ht.select(
        syn=ht.constraint_groups[0]
        .annotate(**ht.constraint_groups[0].oe_info[0])
        .drop("oe_info"),
        mis=ht.constraint_groups[1]
        .annotate(**ht.constraint_groups[1].oe_info[0])
        .drop("oe_info"),
        lof=ht.constraint_groups[5]
        .annotate(**ht.constraint_groups[5].oe_info[0])
        .drop("oe_info"),
    )

    ht = ht.key_by("transcript")

    return ht


def process_context_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Read, re-key, and restructure the context preprocessed Hail Table.

    This function:
    1. Reads the context table from `get_context_preprocessed_ht`.
    2. Extracts population labels from the `exomes_freq_meta` field.
    3. Re-keys the table by (`locus`, `alleles`).
    4. Selects and renames key fields:
       - Coverage statistics (`mean`, `median_approx`, `AN`, `percent_AN`)
       - Calibrated mutation frequencies per population

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    ht = get_context_preprocessed_ht(version).ht()
    pops = hl.eval(ht.exomes_freq_meta.map(lambda x: x.get("gen_anc", "total")))

    ht = ht.select(
        "context",
        "ref",
        "alt",
        "was_flipped",
        "transition",
        "cpg",
        "mutation_type",
        "methylation_level",
        transcript_consequences=ht.vep.transcript_consequences.map(
            lambda x: x.annotate(
                canonical=hl.or_else(x.canonical == 1, False),
                mane_select=hl.is_defined(x.mane_select),
                vep_domains=x.domains,
                residue_ref=x.amino_acids.split("/").first(),
                residue_alt=x.amino_acids.split("/").last(),
            ).select(
                "transcript_id",
                "gene_id",
                "gene_symbol",
                "canonical",
                "mane_select",
                "biotype",
                "most_severe_consequence",
                "sift_score",
                "polyphen_score",
                "vep_domains",
                "residue_ref",
                "residue_alt",
            )
        ),
        gnomad_exomes_filters=ht.filters.exomes,
        gnomad_exomes_coverage=hl.struct(
            mean=ht.coverage.exomes.mean,
            median_approx=ht.coverage.exomes.median_approx,
            AN=ht.AN.exomes,
            percent_AN=ht.exomes_coverage,
        ),
        gnomad_exomes_freq=hl.struct(
            **{pop: ht.calibrate_mu.exomes_freq[i] for i, pop in enumerate(pops)}
        ),
    ).explode("transcript_consequences")
    ht = ht.transmute(**ht.transcript_consequences)
    ht = ht.key_by("locus", "alleles", "transcript_id")

    return ht


def process_genetics_gym_missense_scores_ht() -> hl.Table:
    """
    Process Genetics Gym missense scores Hail Table.

    :return: Hail Table with Genetics Gym missense scores.
    """
    scores = [
        "esm_score",
        "proteinmpnn_llr",
        "am_pathogenicity",
        "rasp_score",
        "MisFit_S",
        "MisFit_D",
        "popeve",
        "eve",
        "esm1_v",
        "mpc",
        "esm_score_neg",
        "proteinmpnn_llr_neg",
        "popeve_neg",
        "esm1_v_neg",
    ]
    ht = (
        get_genetics_gym_missense_scores_ht()
        .ht()
        .select(
            "ensembl_tid",
            "uniprot_id",
            *scores,
        )
    )
    ht = ht.filter(hl.any([hl.is_defined(ht[s]) for s in scores]))
    ht = ht.annotate(
        uniprot_id=hl.or_else(ht.uniprot_id, "None"), rand_n=hl.rand_unif(0, 1)
    ).cache()

    for s in scores:
        ht = ht.order_by(ht[s], ht.rand_n)
        ht = ht.annotate(
            **{
                f"{s}_idx": hl.or_missing(
                    hl.is_defined(ht[s]), hl.scan.count_where(hl.is_defined(ht[s]))
                )
            }
        ).cache()

    ht = ht.annotate(transcript_id=ht.ensembl_tid)
    ht = ht.key_by("locus", "alleles", "transcript_id", "uniprot_id").cache()
    max_idx = ht.aggregate(hl.struct(**{s: hl.agg.max(ht[f"{s}_idx"]) for s in scores}))
    ht = ht.select(
        **{
            s: hl.struct(
                score=ht[s],
                idx=ht[f"{s}_idx"],
                percentile=hl.int((ht[f"{s}_idx"] / max_idx[s]) * 100),
            )
            for s in scores
        }
    )

    return ht


def import_revel_ht() -> hl.Table:
    """
    Import REVEL Hail Table.

    :return: Hail Table with REVEL annotations.
    """
    ht = hl.import_table(
        get_revel_csv(),
        delimiter=",",
        min_partitions=1000,
        types={"grch38_pos": hl.tstr, "REVEL": hl.tfloat64},
    )

    ht = ht.drop("hg19_pos", "aaref", "aaalt")

    # Drop variants that have no position in GRCh38 when lifted over from GRCh37.
    ht = ht.filter(ht.grch38_pos.contains("."), keep=False)
    ht = ht.transmute(chr="chr" + ht.chr)
    ht = ht.select(
        locus=hl.locus(ht.chr, hl.int(ht.grch38_pos), reference_genome="GRCh38"),
        alleles=hl.array([ht.ref, ht.alt]),
        revel=ht.REVEL,
        transcript_id=ht.Ensembl_transcriptid.strip().split(";"),
    )
    ht = ht.explode("transcript_id")
    ht = ht.key_by("locus", "alleles", "transcript_id")

    return ht


def main(args):
    """Execute the Proemis 3D pipeline."""
    hl.init(
        log="/proemis_3d_data_import.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    overwrite = args.overwrite

    if args.import_cosmis_score_data or args.import_all:
        for model in ["alphafold", "swiss_model", "pdb"]:
            import_cosmis_score_data(model).write(
                get_cosmis_score_ht(model).path,
                overwrite=overwrite,
            )

    if args.import_varity_data or args.import_all:
        import_varity_data().write(
            get_varity_ht().path,
            overwrite=overwrite,
        )

    if args.import_mtr3d_data or args.import_all:
        import_mtr3d_data().write(
            get_mtr3d_ht().path,
            overwrite=overwrite,
        )

    if args.import_interpro_annotations or args.import_all:
        import_interpro_annotations().write(
            get_interpro_annotations_ht().path,
            overwrite=overwrite,
        )

    if args.import_kaplanis_variants or args.import_all:
        import_kaplanis_variants().write(
            get_kaplanis_variants_ht().path,
            overwrite=overwrite,
        )
        import_kaplanis_variants(
            liftover_to_grch38=True,
            key_by_gene_and_transcript=True,
        ).write(
            get_kaplanis_variants_ht(
                liftover_to_grch38=True, key_by_transcript=True
            ).path,
            overwrite=overwrite,
        )
        import_kaplanis_variants(
            liftover_to_grch38=True,
            key_by_gene_and_transcript=False,
        ).write(
            get_kaplanis_variants_ht(
                liftover_to_grch38=True, key_by_transcript=False
            ).path,
            overwrite=overwrite,
        )
        import_kaplanis_variants(
            liftover_to_grch38=False,
            key_by_gene_and_transcript=True,
        ).write(
            get_kaplanis_variants_ht(
                liftover_to_grch38=False, key_by_transcript=True
            ).path,
            overwrite=overwrite,
        )

    if args.import_fu_variants or args.import_all:
        import_fu_variants().write(
            get_fu_variants_ht().path,
            overwrite=overwrite,
        )

    if args.import_revel_ht or args.import_all:
        import_revel_ht().write(
            get_insilico_annotations_ht("revel").path,
            overwrite=overwrite,
        )

    if args.process_clinvar_ht or args.import_all:
        process_clinvar_ht(args.clinvar_version).write(
            get_clinvar_missense_ht().path,
            overwrite=overwrite,
        )

    if args.process_constraint_metrics_ht or args.import_all:
        process_constraint_metrics_ht(args.constraint_metrics_ht_version).write(
            get_temp_processed_constraint_ht(args.constraint_metrics_ht_version).path,
            overwrite=overwrite,
        )

    if args.import_mtr_data or args.import_all:
        import_mtr_data().write(
            get_mtr_ht().path,
            overwrite=overwrite,
        )

    if args.process_rmc_ht or args.import_all:
        process_rmc_ht(args.rmc_ht_version).write(
            get_temp_processed_rmc_ht(args.rmc_ht_version).path,
            overwrite=overwrite,
        )

    if args.process_context_ht or args.import_all:
        process_context_ht(args.context_ht_version).write(
            get_temp_context_preprocessed_ht(args.context_ht_version).path,
            overwrite=overwrite,
        )

    if args.process_genetics_gym_missense_scores_ht or args.import_all:
        process_genetics_gym_missense_scores_ht().write(
            get_processed_genetics_gym_missense_scores_ht().path,
            overwrite=overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--import-all",
        help="Import all data.",
        action="store_true",
    )
    parser.add_argument(
        "--import-cosmis-score-data",
        help="Whether to import COSMIS score data.",
        action="store_true",
    )
    parser.add_argument(
        "--import-varity-data",
        help="Whether to import Varity data.",
        action="store_true",
    )
    parser.add_argument(
        "--import-mtr3d-data",
        help="Whether to import MTR3D data.",
        action="store_true",
    )
    parser.add_argument(
        "--import-interpro-annotations",
        help="Whether to import InterPro annotations.",
        action="store_true",
    )
    parser.add_argument(
        "--import-kaplanis-variants",
        help="Whether to import Kaplanis variants.",
        action="store_true",
    )
    parser.add_argument(
        "--import-fu-variants",
        help="Whether to import Fu variants.",
        action="store_true",
    )
    parser.add_argument(
        "--process-clinvar-ht",
        help="Whether to process ClinVar HT.",
        action="store_true",
    )
    parser.add_argument(
        "--import-revel-ht",
        help="Whether to import REVEL HT.",
        action="store_true",
    )
    parser.add_argument(
        "--clinvar-version",
        help="Version of ClinVar HT to process.",
        default="20250504",
    )
    parser.add_argument(
        "--process-constraint-metrics-ht",
        help="Whether to process constraint metrics HT.",
        action="store_true",
    )
    parser.add_argument(
        "--constraint-metrics-ht-version",
        help="Version of constraint metrics HT to process.",
        default=CURRENT_VERSION,
    )
    parser.add_argument(
        "--import-mtr-data",
        help="Whether to import MTR data.",
        action="store_true",
    )
    parser.add_argument(
        "--process-rmc-ht",
        help="Whether to process the RMC HT.",
        action="store_true",
    )
    parser.add_argument(
        "--rmc-ht-version",
        help="Version of RMC HT to process.",
        default=CURRENT_VERSION,
    )
    parser.add_argument(
        "--process-context-ht",
        help="Whether to process the context HT.",
        action="store_true",
    )
    parser.add_argument(
        "--context-ht-version",
        help="Version of context HT to process.",
        default=CURRENT_VERSION,
    )
    parser.add_argument(
        "--process-genetics-gym-missense-scores-ht",
        help="Whether to process the Genetics Gym missense scores HT.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
