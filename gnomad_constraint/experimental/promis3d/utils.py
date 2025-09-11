"""Script with utility functions for the Promis3D pipeline."""

import json
import logging
import os
from io import StringIO
from typing import Dict, Iterator, List, Optional, Union

import hail as hl
import numpy as np
import pandas as pd
from Bio.PDB import PPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from gnomad.resources.grch38.gnomad import browser_gene, browser_variant, pext
from gnomad.utils.constraint import oe_confidence_interval
from gnomad.utils.filtering import filter_gencode_ht
from gnomad.utils.reference_genome import get_reference_genome
from hail.utils.misc import divide_null
from pyspark.sql.functions import col, explode, pandas_udf, rtrim, split
from pyspark.sql.types import StringType, StructField, StructType

from gnomad_constraint.experimental.promis3d.constants import (
    HI_GENE_CATEGORIES,
    HI_GENES,
    MIN_EXP_MIS,
)
from gnomad_constraint.experimental.promis3d.data_import import (
    get_kaplanis_sig_gene_annotations,
    process_gnomad_de_novo_ht,
    process_gnomad_site_ht,
    process_pext_annotation_ht,
    process_pext_base_ht,
)
from gnomad_constraint.experimental.promis3d.resources import (
    get_clinvar_missense_ht,
    get_cosmis_score_ht,
    get_fu_variants_ht,
    get_gnomad_de_novo_ht,
    get_insilico_annotations_ht,
    get_interpro_annotations_ht,
    get_kaplanis_variants_ht,
    get_mtr3d_ht,
    get_mtr_ht,
    get_phaplo,
    get_processed_genetics_gym_missense_scores_ht,
    get_ptriplo,
    get_rmc_browser_ht,
    get_temp_context_preprocessed_ht,
    get_temp_processed_constraint_ht,
    get_temp_processed_rmc_ht,
    get_varity_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("promis3d_utils")
logger.setLevel(logging.INFO)

########################################################################################
# Functions to perform tasks from convert_gencode_fastn_to_dt.R and
# convert_gencode_fasta_to_dt.R
########################################################################################
COLNAMES_TRANSCRIPTS = [
    "enst",
    "ensg",
    "havana_g",
    "havana_t",
    "transcript",
    "gene",
    "ntlength",
    "index1",
    "index2",
    "index3",
]
"""
Column names for the GENCODE transcripts Hail Table.
"""

COLNAMES_TRANSLATIONS = {
    "2.1.1": [
        "enst",
        "ensg",
        "havana_g",
        "havana_t",
        "transcript",
        "gene",
        "aalength",
    ],
    "4.1": [
        "ensp",
        "enst",
        "ensg",
        "havana_g",
        "havana_t",
        "transcript",
        "gene",
        "aalength",
    ],
}
"""
Column names for the GENCODE translations Hail Table.
"""

VARIANT_LEVEL_ANNOTATION_CONFIG = {
    "context": {
        "ht": get_temp_context_preprocessed_ht(),
        "keys": ["locus", "alleles", "transcript_id"],
    },
    "gnomad_site": {
        "ht": browser_variant(),
        "keys": ["locus", "alleles"],
        "custom_select": process_gnomad_site_ht,
    },
    "revel": {
        "ht": get_insilico_annotations_ht("revel"),
        "keys": ["locus", "alleles", "transcript_id"],
    },
    "cadd": {
        "ht": get_insilico_annotations_ht("cadd"),
        "keys": ["locus", "alleles"],
    },
    "phylop": {"ht": get_insilico_annotations_ht("phylop"), "keys": ["locus"]},
    "genetics_gym": {
        "ht": get_processed_genetics_gym_missense_scores_ht(),
        "keys": ["locus", "alleles", "transcript_id", "uniprot_id"],
        "annotation_name": "genetics_gym_missense_scores",
    },
    "autism": {
        "ht": get_fu_variants_ht(),
        "keys": ["locus", "alleles"],
        "annotation_name": "autism",
    },
    "dd_denovo": {
        "ht": get_kaplanis_variants_ht(liftover_to_grch38=True, key_by_transcript=True),
        "keys": ["locus", "alleles", "gene_id", "transcript_id"],
        "annotation_name": "dd_denovo",
    },
    "dd_denovo_no_transcript": {
        "ht": get_kaplanis_variants_ht(liftover_to_grch38=True),
        "keys": ["locus", "alleles"],
        "annotation_name": "dd_denovo_no_transcript_match",
    },
    "gnomad_de_novo": {
        "ht": get_gnomad_de_novo_ht(),
        "keys": ["locus", "alleles"],
        "annotation_name": "gnomad_de_novo",
        "custom_select": process_gnomad_de_novo_ht,
    },
    "clinvar": {
        "ht": get_clinvar_missense_ht(),
        "keys": ["locus", "alleles", "gene_symbol"],
        "annotation_name": "clinvar",
    },
    "pext_base": {
        "ht": pext("base_level"),
        "keys": ["locus", "gene_id"],
        "annotation_name": "base_level_pext",
        "custom_select": process_pext_base_ht,
    },
    "mtr": {
        "ht": get_mtr_ht(),
        "keys": ["locus", "alleles", "transcript_id"],
        "annotation_name": "mtr",
    },
    "rmc": {
        "ht": get_temp_processed_rmc_ht(),
        "keys": ["locus", "transcript_id"],
        "annotation_name": "rmc",
    },
}
"""
Configuration for variant level annotations.
"""

RESIDUE_LEVEL_ANNOTATION_CONFIG = {
    "interpro": {
        "ht": get_interpro_annotations_ht(),
        "keys": ["transcript_id", "uniprot_id", "residue_index"],
        "annotation_name": "interpro",
    },
    "varity": {
        "ht": get_varity_ht(),
        "keys": ["uniprot_id", "residue_index", "residue_ref", "residue_alt"],
        "annotation_name": "varity",
    },
    "mtr3d": {
        "ht": get_mtr3d_ht(),
        "keys": ["transcript_id", "uniprot_id", "residue_index"],
        "annotation_name": "mtr3d",
    },
    "cosmis_alphafold": {
        "ht": get_cosmis_score_ht("alphafold"),
        "keys": ["transcript_id", "uniprot_id", "residue_index"],
        "annotation_name": "cosmis_alphafold",
    },
    "cosmis_pdb": {
        "ht": get_cosmis_score_ht("pdb"),
        "keys": ["transcript_id", "uniprot_id", "residue_index"],
        "annotation_name": "cosmis_pdb",
    },
    "cosmis_swiss_model": {
        "ht": get_cosmis_score_ht("swiss_model"),
        "keys": ["transcript_id", "uniprot_id", "residue_index"],
        "annotation_name": "cosmis_swiss_model",
    },
}
"""
Configuration for residue level annotations.
"""

BASE_LEVEL_ANNOTATION_FIELDS = [
    "gene_symbol",
    "canonical",
    "mane_select",
    "transcript_biotype",
    "most_severe_consequence",
]
"""
Fields to keep at the base level.
"""


def convert_fasta_to_table(fasta_file: str, colnames: List[str]) -> hl.Table:
    """
    Convert a FASTA file to a Hail Table.

    :param fasta_file: Path to the FASTA file.
    :param colnames: Column names for the Hail Table.
    :return: Hail Table with the FASTA file contents.
    """
    spark = hl.utils.java.Env.spark_session()
    df = spark.read.format("text").load(fasta_file, wholetext=True)
    df = df.select(explode(split(df["value"], ">")).alias("sequence"))

    # Convert the Spark DataFrame to a Hail Table.
    ht = hl.Table.from_spark(df)
    ht = ht.filter(ht.sequence != "")

    split_expr = ht.sequence.split("\n")
    split_info_expr = split_expr[0].split("\\|")
    ht = ht.select(
        **{
            n: hl.or_missing(
                (split_info_expr.length() > i)
                & ~hl.array(["", "-"]).contains(split_info_expr[i]),
                split_info_expr[i],
            )
            for i, n in enumerate(colnames)
        },
        sequence=split_expr[1].upper(),
    )

    # Remove version numbers from ENST and ENSG.
    ht = ht.annotate(enst=ht.enst.split("\\.")[0], ensg=ht.ensg.split("\\.")[0])

    return ht


def convert_gencode_transcripts_fasta_to_table(fasta_file: str) -> hl.Table:
    """
    Convert GENCODE transcripts FASTA file to a Hail Table.

    :param fasta_file: Path to the GENCODE transcripts FASTA file.
    :return: Hail Table with the GENCODE transcripts FASTA file contents parsed.
    """
    ht = convert_fasta_to_table(fasta_file, COLNAMES_TRANSCRIPTS)

    # Organize the UTR5/CDS/UTR3 indices.
    # If the UTR5/CDS/UTR3 index is not present, it is set to missing.
    ht = ht.annotate(
        index1=hl.or_missing(ht.index1.startswith("UTR5"), ht.index1),
        index2=(
            hl.case()
            .when(ht.index1.startswith("CDS"), ht.index1)
            .when(ht.index2.startswith("CDS"), ht.index2)
            .or_missing()
        ),
        index3=(
            hl.case()
            .when(ht.index1.startswith("UTR3"), ht.index1)
            .when(ht.index2.startswith("UTR3"), ht.index2)
            .when(ht.index3.startswith("UTR3"), ht.index3)
            .or_missing()
        ),
    )

    # Rename the indices to 'utr5', 'cds', and 'utr3'.
    names = ["utr5", "cds", "utr3"]
    ht = ht.rename({f"index{i + 1}": n for i, n in enumerate(names)})

    # Split the 'utr5', 'cds', and 'utr3' annotations into start and end positions.
    ht = ht.annotate(
        **{
            n: hl.bind(lambda x: x.map(hl.int), ht[n].split(":")[1].split("-"))
            for n in names
        }
    )

    # Trim the sequence to the CDS range.
    ht = ht.annotate(cds_sequence=ht.sequence[ht.cds[0] - 1 : ht.cds[1]])

    return ht


########################################################################################
# Note that the functionality in split_context_obs_exp.R is not implemented here
# because it just splits the tsv file by transcript, we will be directly using the
# Hail Table instead.
########################################################################################


########################################################################################
# Functions to perform tasks from read_af2_sequences.R
########################################################################################
def get_plddt_from_confidence_json(plddt_content: str) -> List[float]:
    """
    Get the pLDDT from a confidence JSON file.

    :param plddt_content: Content of the pLDDT JSON file.
    :return: List of pLDDT scores.
    """
    data = json.loads(plddt_content)
    return data["confidenceScore"]


def get_pae_from_json(pae_content: str) -> List[List[int]]:
    """
    Get the PAE from a PAE JSON file.

    :param pae_content: Content of the PAE JSON file.
    :return: List of PAE scores.
    """
    data = json.loads(pae_content)
    return data[0]["predicted_aligned_error"]


def get_structure_peptide(structure) -> str:
    """
    Get the sequence from a structure.

    :param structure: Structure object.
    :return: Sequence as a string.
    """
    ppb = PPBuilder()

    # Return the sequence as a string.
    return "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(structure)])


def get_structure_dist_matrix(structure: MMCIFParser) -> np.ndarray:
    """
    Calculate the "calpha" distance matrix from a structure.

    :param structure: Structure object.
    :return: Distance matrix as a NumPy array.
    """
    calpha_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True) and "CA" in residue:
                    calpha_atoms.append(residue["CA"].get_coord())

    def _calc_dist_matrix(calpha_atoms: List[np.ndarray]) -> np.ndarray:
        """
        Calculate the pairwise distance matrix between Calpha atoms.

        :param calpha_atoms: List of Calpha atoms.
        :return: Distance matrix as a NumPy array.
        """
        num_atoms = len(calpha_atoms)
        dist_matrix = np.zeros((num_atoms, num_atoms))
        for i, atom1 in enumerate(calpha_atoms):
            for j, atom2 in enumerate(calpha_atoms):
                dist_matrix[i, j] = np.linalg.norm(atom1 - atom2)

        return dist_matrix

    # Calculate the distance matrix
    return _calc_dist_matrix(calpha_atoms)


def process_af2_mmcif(
    uniprot_id: str,
    mmcif_content: str,
    distance_matrix: bool = False,
) -> Union[str, np.ndarray, List[float]]:
    """
    Process AlphaFold2 MMCIF content.

    :param uniprot_id: UniProt ID.
    :param mmcif_content: MMCIF content as a string.
    :param distance_matrix: Whether to return the distance matrix. Default is False.
    :return: Sequence or distance matrix.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(uniprot_id, StringIO(mmcif_content))

    if distance_matrix:
        return get_structure_dist_matrix(structure)
    else:
        return get_structure_peptide(structure)


def process_af2_file_by_mode(
    uniprot_id: str,
    file_content: str,
    mode: str,
) -> Union[str, np.ndarray, List[float]]:
    """
    Dispatcher to handle different AF2 modes based on filename suffix and mode.

    :param uniprot_id: UniProt ID.
    :param file_content: File content.
    :param mode: Mode for processing files. Options are 'sequence', 'distance_matrix',
        'plddt', or 'pae'. Default is 'sequence'.
    :return: Sequence or distance matrix.
    """
    if mode in {"sequence", "distance_matrix"}:
        return process_af2_mmcif(
            uniprot_id, file_content, distance_matrix=(mode == "distance_matrix")
        )
    if mode == "plddt":
        return get_plddt_from_confidence_json(file_content)

    if mode == "pae":
        return get_pae_from_json(file_content)

    raise ValueError(f"Unsupported mode: {mode}")


def process_af2_structures(
    gcs_bucket_glob: str,
    mode: str = "sequence",
) -> hl.Table:
    """
    Process AlphaFold2 structures from a GCS bucket.

    .. note::

        All files in the bucket must be in CIF format with a '.cif.gz' extension.

    :param gcs_bucket_glob: GCS bucket glob pattern.
    :param mode: Mode for processing files. Options are 'sequence', 'distance_matrix',
        'plddt', or 'pae'. Default is 'sequence'.
    :return: Hail Table with UniProt IDs and sequences or distance matrices.
    """
    # Get Spark session for file distribution and processing.
    spark = hl.utils.java.Env.spark_session()
    spark.conf.set(
        "spark.sql.execution.arrow.maxRecordsPerBatch",
        1000,
    )

    # Define schema for loading the files.
    schema = StructType(
        [
            StructField("file_content", StringType(), True),
            StructField("af2_file", StringType(), True),
        ]
    )

    # Use Spark to read files in parallel.
    # This reads the entire content of each file as a (filename, content) pair.
    af2_files_df = (
        spark.read.format("text")
        .load(gcs_bucket_glob, schema=schema, wholetext=True)
        .withColumn("af2_file", col("_metadata.file_path"))
    )
    if mode == "distance_matrix":
        col_name = "dist_mat"
        rtype = "array<array<float>>"
    elif mode == "pae":
        col_name = "pae"
        rtype = "array<array<int>>"
    elif mode == "plddt":
        col_name = "plddt"
        rtype = "array<float>"
    else:
        col_name = "sequence"
        rtype = "string"

    @pandas_udf(f"uniprot_id string, {col_name} {rtype}")
    def process_file(
        file_path_series: pd.Series, file_content_series: pd.Series
    ) -> pd.DataFrame:
        """
        Process a list of files in parallel using a Pandas UDF.

        :param file_path_series: File paths.
        :param file_content_series: File contents.
        :return: Pandas DataFrame with UniProt IDs and sequences.
        """
        result = []
        for file_path, file_content in zip(file_path_series, file_content_series):
            # Extract UniProt ID from the file path.
            uniprot_id = os.path.basename(file_path).split("-")[1]

            # Process the file content.
            # af2_data = process_af2_file_by_mode(uniprot_id, file_content, mode=mode)
            af2_data = None
            result.append((uniprot_id, af2_data))

        return pd.DataFrame(result, columns=["uniprot_id", col_name])

    # Apply the Pandas UDF to process the files.
    result_df = af2_files_df.withColumn(
        "uniprot_id_sequence", process_file(col("af2_file"), col("file_content"))
    )

    # Split the resulting column into separate columns.
    result_df = result_df.select(
        "af2_file",
        col("uniprot_id_sequence.uniprot_id"),
        col(f"uniprot_id_sequence.{col_name}"),
    )
    if mode in {"distance_matrix", "pae"}:
        from pyspark.sql.functions import posexplode

        result_df = result_df.select(
            "af2_file",
            "uniprot_id",
            posexplode(col(col_name)).alias("aa_index", col_name),
        )
    tmp_path = hl.utils.new_temp_file("process_af2_structures", "parquet")
    logger.info(f"Writing processed AlphaFold2 structures to {tmp_path}")
    result_df.write.mode("overwrite").save(tmp_path)
    logger.info(f"Finished writing.")
    result_df = spark.read.parquet(tmp_path)

    # Convert the Spark DataFrame to a Hail Table.
    key = ["af2_file", "uniprot_id"]
    if mode in {"distance_matrix", "pae"}:
        key.append("aa_index")

    ht = hl.Table.from_spark(result_df, key=key)

    return ht


def remove_multi_frag_uniprots(ht: hl.Table) -> hl.Table:
    """
    Remove UniProt IDs with multiple fragments (F2).

    .. note::

        With the current release, there is only one structure per UniProt.

    :param ht: Hail Table with structures keyed by 'af2_file' and 'uniprot_id'.
    :return: Hail Table with UniProt IDs with multiple fragments removed and
        keyed by 'uniprot_id'.
    """
    uniprots_with_multifrags = (
        ht.filter(ht.af2_file.contains("-F2-")).select().distinct()
    )

    return ht.anti_join(uniprots_with_multifrags).key_by("uniprot_id").drop("af2_file")


########################################################################################
# Functions to perform tasks from gencode_alignment.R
########################################################################################
def join_by_sequence(ht1: hl.Table, ht2: hl.Table) -> hl.Table:
    """
    Join two Hail Tables based on the 'sequence' field.

    .. note::

        The join is an inner join.

    :param ht1: First Hail Table.
    :param ht2: Second Hail Table.
    :return: Hail Table with the two input tables joined based on 'sequence'.
    """
    # Overlap the tables based on sequence.
    ht1 = ht1.key_by("sequence")
    ht2 = ht2.key_by("sequence")

    return ht1.join(ht2, how="inner")


########################################################################################
# Functions to perform tasks from make_gencode_positions_files.R
########################################################################################
def get_gencode_positions(
    transcripts_ht: hl.Table,
    translations_ht: hl.Table,
    gencode_gtf_ht: hl.Table,
    no_filter: bool = False,
) -> hl.Table:
    """
    Get GENCODE positions for the given transcripts and translations.

    :param transcripts_ht: Hail Table with GENCODE transcripts.
    :param translations_ht: Hail Table with GENCODE translations.
    :param gencode_gtf_ht: Hail Table with GENCODE GTF data.
    :param no_filter: Whether to filter the data if the CDS length is not divisible
        by 3 or if the sequence length in `transcripts_ht` is not equal to the sequence
        length in `gencode_gtf_ht`. If `no_filter` is True, then two additional fields
        are added to the table: `cds_len_mismatch` and `cds_len_not_div_by_3` to
        facilitate filtering. Default is False.
    :return: Hail Table with GENCODE positions.
    """
    build = get_reference_genome(gencode_gtf_ht.interval.start).name

    # Filter GTF data to keep only CDS features.
    gencode_gtf_ht = gencode_gtf_ht.filter(gencode_gtf_ht.feature == "CDS")

    # Get list of intervals for each transcript, and keep strand information.
    gencode_gtf_ht = gencode_gtf_ht.group_by("transcript_id", "strand").aggregate(
        intervals=hl.agg.collect(gencode_gtf_ht.interval)
    )

    # Get CDS positions and lengths for each transcript in the GTF data.
    positions = gencode_gtf_ht.intervals.flatmap(
        lambda x: hl.range(x.start.position, x.end.position + 1)
    )
    gencode_gtf_ht = gencode_gtf_ht.transmute(
        chrom=gencode_gtf_ht.intervals[0].start.contig,
        gtf_cds_pos=hl.sorted(positions),
        gtf_cds_len=hl.len(positions),
    ).key_by("transcript_id")

    # Get CDS sequences and lengths for each transcript in the transcript data.
    transcripts_ht = transcripts_ht.select(
        "enst", "cds_sequence", cds_len=hl.len(transcripts_ht.cds_sequence)
    ).key_by("enst")

    # Join the CDS data from the GTF and transcript data.
    cds_len_ht = transcripts_ht.join(gencode_gtf_ht, how="inner")

    # Filter CDS data to keep only transcripts with matching CDS lengths after removing
    # stop codons.
    cds_sequence_expr = hl.if_else(
        cds_len_ht.gtf_cds_len == (cds_len_ht.cds_len - 3),
        cds_len_ht.cds_sequence[:-3],
        cds_len_ht.cds_sequence,
    )
    cds_len_ht = cds_len_ht.annotate(
        cds_sequence=cds_sequence_expr,
        cds_len=hl.len(cds_sequence_expr),
    )
    cds_len_mismatch_expr = cds_len_ht.gtf_cds_len != cds_len_ht.cds_len
    cds_len_not_div_by_3_expr = cds_len_ht.cds_len % 3 != 0
    if no_filter:
        cds_len_ht = cds_len_ht.annotate(
            cds_len_mismatch=cds_len_mismatch_expr,
            cds_len_not_div_by_3=cds_len_not_div_by_3_expr,
        )
    else:
        cds_len_ht = cds_len_ht.filter(
            ~cds_len_mismatch_expr & ~cds_len_not_div_by_3_expr
        )

    # Get the reference sequence and amino acid positions for each position in the CDS.
    # If the strand is negative, reverse the reference sequence and amino acid
    # positions.
    aapos_expr = hl.flatten(hl.repeat(hl.range(hl.int(cds_len_ht.cds_len / 3)), 3))
    cds_len_ht = cds_len_ht.annotate(
        ref=hl.if_else(
            cds_len_ht.strand == "+",
            cds_len_ht.cds_sequence,
            hl.expr.functions.reverse_complement(cds_len_ht.cds_sequence, rna=False),
        ),
        aapos=hl.if_else(
            cds_len_ht.strand == "+",
            hl.sorted(aapos_expr),
            hl.sorted(aapos_expr, reverse=True),
        ),
    )

    # Annotate the translations data with the genomic positions, reference sequence, and
    # amino acid positions.
    # Explode the positions to get one row per position.
    ht = translations_ht.annotate(**cds_len_ht[translations_ht.enst])
    ht = ht.transmute(
        positions=hl.map(
            lambda gp, r, ap: hl.struct(
                locus=hl.locus(ht.chrom, gp, reference_genome=build), ref=r, aapos=ap
            ),
            ht.gtf_cds_pos,
            hl.range(hl.len(ht.ref)).map(lambda i: ht.ref[i]),
            ht.aapos,
        )
    ).explode("positions")

    return ht.transmute(**ht.positions)


########################################################################################
# Functions to perform tasks from run_greedy.R and run_forward.R
########################################################################################
def generate_codon_oe_table(obs_exp_ht: hl.Table, pos_ht: hl.Table) -> hl.Table:
    """
    Generate a Table with observed and expected values for codons.

    :param obs_exp_ht: Hail Table with observed and expected values.
    :param pos_ht: Hail Table with positions.
    :return: Hail Table with observed and expected values for codons.
    """
    oe_keyed = obs_exp_ht[pos_ht.locus, pos_ht.enst]
    pos_ht = pos_ht.annotate(obs=oe_keyed.obs, exp=oe_keyed.exp)

    # Get aggregate sum of observed and expected values for each codon.
    oe_codon_ht = pos_ht.group_by("enst", "uniprot_id", "aapos").aggregate(
        obs=hl.agg.sum(pos_ht.obs),
        exp=hl.agg.sum(pos_ht.exp),
    )

    # Get a list of observed and expected codon values for each transcript and UniProt
    # ID sorted by amino acid position.
    oe_codon_ht = oe_codon_ht.group_by("enst", "uniprot_id").aggregate(
        oe=hl.agg.collect(
            (
                oe_codon_ht.aapos,
                hl.struct(obs=oe_codon_ht.obs, exp=oe_codon_ht.exp),
            )
        )
    )
    oe_codon_ht = oe_codon_ht.annotate(
        oe=hl.sorted(oe_codon_ht.oe, key=lambda x: x[0]).map(lambda x: x[1])
    ).key_by("uniprot_id")

    return oe_codon_ht.collect_by_key("oe_by_transcript")


def add_idx_to_array(
    expr: hl.expr.ArrayExpression, idx_name: str, element_name: Optional[str] = None
) -> hl.expr.ArrayExpression:
    """
    Add an index to each element in an array expression.

    If the elements are structs, the index is added as a field with the name `idx_name`.
    If the elements are not structs, a new struct is created with the index as a field
    with the name `idx_name` and the element as a field with the name `element_name`.
    If `element_name` is not provided, then only hl.enumerate is used to add the index.

    :param expr: Array expression to add index to.
    :param idx_name: Name of the index field.
    :param element_name: Name of the element field. Default is None.
    :return: Array expression with index added to each element.
    """
    element_type = expr.dtype.element_type
    expr = hl.enumerate(expr)

    if isinstance(element_type, hl.tstruct):
        return expr.map(lambda x: x[1].annotate(**{idx_name: x[0]}))
    elif element_name is not None:
        return expr.map(lambda x: hl.struct(**{element_name: x[1], idx_name: x[0]}))
    else:
        return expr


def get_cumulative_oe(oe_expr):
    """
    Get the cumulative OE.

    :param oe_expr: Array expression with observed and expected values.
    :return: Array expression with cumulative OE.
    """
    oe_expr = hl.array_scan(
        lambda i, j: j.annotate(obs=i.obs + j.obs, exp=i.exp + j.exp),
        oe_expr[0],
        oe_expr[1:],
    )

    return oe_expr


def calculate_oe_upper(oe_expr, alpha=0.05):
    """
    Calculate the upper bound of the OE confidence interval.

    :param oe_expr: Array expression with observed and expected values.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :return: Array expression with upper bound of the OE confidence interval.
    """
    # Calculate upper bound of oe confidence interval.
    oe_upper_expr = oe_expr.map(
        lambda x: x.annotate(
            oe=divide_null(x.obs, x.exp),
            oe_upper=(
                hl.qchisqtail(1 - alpha / 2, 2 * (x.obs + 1), lower_tail=True)
                / (2 * x.exp)
            ),
        )
    )

    return oe_upper_expr


def get_min_oe_upper(oe_expr, min_exp_mis=None):
    """
    Get the 3D residue with the lowest upper bound of the OE confidence interval.

    :param oe_expr: Array expression with observed and expected values.
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is None.
    :return: Struct expression with the 3D residue with the lowest upper bound of the OE
        confidence interval.
    """
    oe_expr = add_idx_to_array(oe_expr, "dist_index")
    if min_exp_mis is None:
        filtered_oe_expr = oe_expr
    else:
        filtered_oe_expr = oe_expr.filter(lambda x: x.exp >= min_exp_mis)
        filtered_oe_expr = hl.or_missing(
            filtered_oe_expr.length() > 0, filtered_oe_expr
        )

    min_oe_upper_expr = hl.sorted(filtered_oe_expr, key=lambda x: x.oe_upper)[0]
    dist_index_expr = min_oe_upper_expr.dist_index
    oe_expr = hl.or_missing(
        hl.is_defined(filtered_oe_expr), oe_expr[: dist_index_expr + 1]
    )
    min_oe_upper_expr = min_oe_upper_expr.drop("dist_index")
    min_oe_upper_expr = min_oe_upper_expr.annotate(
        region=oe_expr.map(lambda x: x.residue_index)
    )

    return min_oe_upper_expr


def get_3d_residue(
    dist_mat_expr: hl.expr.ArrayExpression,
    oe_expr: hl.expr.ArrayExpression,
    alpha: float = 0.05,
    min_exp_mis: int = MIN_EXP_MIS,
) -> hl.expr.StructExpression:
    """
    Get the 3D residue with the lowest upper bound of the OE confidence interval.

    :param dist_mat_expr: Array expression with distance matrix.
    :param oe_expr: Array expression with observed and expected values.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is MIN_EXP_MIS.
    :return: Struct expression with the 3D residue with the lowest upper bound of the OE
        confidence interval.
    """
    # Annotate neighbor order per residue.
    dist_mat_expr = add_idx_to_array(
        dist_mat_expr, "residue_index", element_name="dist"
    )
    dist_mat_expr = hl.sorted(dist_mat_expr, key=lambda x: x.dist)
    dist_mat_expr = dist_mat_expr.map(lambda x: x.drop("dist"))

    # Annotate neighbor observed and expected, cumulative observed and expected, and
    # upper bound of OE confidence interval.
    oe_expr = dist_mat_expr.map(lambda x: x.annotate(**oe_expr[x.residue_index]))
    oe_expr = get_cumulative_oe(oe_expr)
    oe_expr = calculate_oe_upper(oe_expr, alpha=alpha)

    # Get the 3D region with the lowest upper bound of the OE confidence interval for
    # each residue.
    min_moeuf_expr = get_min_oe_upper(oe_expr, min_exp_mis=min_exp_mis)

    return min_moeuf_expr


def determine_regions_with_min_oe_upper(
    af2_ht: hl.Table, oe_codon_ht: hl.Table, min_exp_mis: int = MIN_EXP_MIS
) -> hl.Table:
    """
    Determine the most intolerant region for each UniProt ID and residue index.

    :param af2_ht: Hail Table with AlphaFold2 data.
    :param oe_codon_ht: Hail Table with observed and expected values for codons.
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is MIN_EXP_MIS.
    :return: Hail Table with the most intolerant region for each UniProt ID and residue
        index
    """
    af2_ht = af2_ht.annotate(oe=oe_codon_ht[af2_ht.uniprot_id].oe_by_transcript)
    af2_ht = af2_ht.explode(af2_ht.oe)
    af2_ht = af2_ht.annotate(**af2_ht.oe)
    af2_ht = af2_ht.transmute(
        transcript_id=af2_ht.enst,
        oe=af2_ht.oe,
        min_oe_upper=get_3d_residue(
            af2_ht.dist_mat, af2_ht.oe, min_exp_mis=min_exp_mis
        ),
    )

    af2_ht = af2_ht.group_by("uniprot_id", "transcript_id").aggregate(
        oe=hl.agg.take(af2_ht.oe, 1)[0],
        min_oe_upper=hl.agg.collect(
            af2_ht.min_oe_upper.annotate(residue_index=af2_ht.aa_index)
        ),
    )

    af2_ht = af2_ht.annotate(
        oe=hl.enumerate(af2_ht.oe).map(lambda x: x[1].annotate(residue_index=x[0])),
        min_oe_upper=hl.sorted(af2_ht.min_oe_upper, key=lambda x: x.oe),
    )

    return af2_ht


########################################################################################
# Functions specific to the greedy algorithm.
########################################################################################
def run_greedy(ht: hl.Table) -> hl.Table:
    """
    Run the greedy algorithm to find the most intolerant region.

    :param ht
    :return: Hail Table with the most intolerant region for each UniProt ID and residue
        index
    """
    min_oe_upper_expr = add_idx_to_array(ht.min_oe_upper, "region_index")
    initial_score_expr = min_oe_upper_expr.map(
        lambda x: hl.missing(min_oe_upper_expr.dtype.element_type)
    )
    score_expr = hl.fold(
        lambda i, j: hl.enumerate(i).map(
            lambda x: hl.coalesce(x[1], hl.or_missing(j.region.contains(x[0]), j))
        ),
        initial_score_expr,
        min_oe_upper_expr,
    )
    score_expr = add_idx_to_array(score_expr, "residue_index")
    ann_keep = ["residue_index", "region_index", "obs", "exp", "oe", "oe_upper"]
    ht = ht.select(score=score_expr.map(lambda x: x.select(*ann_keep)))
    ht = ht.checkpoint(hl.utils.new_temp_file("sort_regions_by_oe", "ht"))
    ht = ht.explode("score")

    ht = ht.select(**ht.score).key_by("uniprot_id", "transcript_id", "residue_index")

    return ht


########################################################################################
# Functions specific to the forward algorithm.
########################################################################################
def annotate_region_with_oe(region_expr, oe_expr):
    """
    Annotate a region with the OE.

    :param region_expr: Region expression.
    :param oe_expr: OE expression.
    :return: Region expression annotated with the OE.
    """
    return region_expr.map(lambda x: oe_expr[x])


def get_agg_oe_for_region(region_expr):
    """
    Get the aggregate OE for a region.

    :param region_expr: Region expression.
    :return: Aggregate OE expression.
    """
    oe_agg_expr = hl.or_missing(
        hl.is_defined(region_expr),
        region_expr.aggregate(
            lambda x: hl.struct(
                obs=hl.agg.sum(x.obs),
                exp=hl.agg.sum(x.exp),
            )
        ),
    )
    gamma_expr = divide_null(oe_agg_expr.obs, oe_agg_expr.exp)

    return oe_agg_expr.annotate(oe=gamma_expr)


def calculate_neg_log_likelihood(region_expr, gamma_expr):
    """
    Calculate the negative log-likelihood of a region.

    :param region_expr: Region expression.
    :param gamma_expr: Gamma expression.
    :return: Negative log-likelihood expression.
    """
    return hl.sum(
        region_expr.map(lambda x: -hl.dpois(x.obs, gamma_expr * x.exp, log_p=True))
    )


def getAIC(region_expr, nll):
    """
    Get the AIC.

    :param region_expr: Region expression.
    :param nll: Negative log-likelihood.
    :return: AIC.
    """
    if isinstance(region_expr, hl.expr.ArrayExpression):
        region_count = region_expr.length()
    else:
        region_count = hl.int(region_expr.region.length() > 0)

    return 2 * region_count + 2 * nll


def remove_residues_from_region(region_expr, remove_region_expr):
    """
    Remove residues from a region.

    :param region_expr: Region expression.
    :param remove_region_expr: Region expression to remove.
    :return: Region expression with residues removed.
    """
    remove_region_residues = hl.set(remove_region_expr.region)
    updated_region_expr = hl.set(region_expr.region).difference(remove_region_residues)
    return hl.or_missing(
        hl.is_defined(region_expr.region) & hl.is_defined(remove_region_expr.region),
        region_expr.annotate(region=hl.array(updated_region_expr)),
    )


def get_min_region(regions_expr, min_field="oe_upper", min_exp_mis=None):
    """
    Get the minimum region.

    :param regions_expr: Regions expression.
    :param min_field: Field to use for sorting. Default is "oe_upper".
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is None.
    :return: Minimum region expression.
    """
    regions_expr = hl.agg.collect(regions_expr)
    if min_exp_mis is None:
        filtered_regions_expr = regions_expr
    else:
        filtered_regions_expr = regions_expr.filter(lambda x: x.exp >= min_exp_mis)
        filtered_regions_expr = hl.or_missing(
            filtered_regions_expr.length() > 0, filtered_regions_expr
        )

    min_region_expr = hl.sorted(filtered_regions_expr, key=lambda x: x[min_field])[0]

    return min_region_expr


def prep_region_struct(region_expr, oe_expr):
    """
    Prepare a region struct.

    :param region_expr: Region expression.
    :param oe_expr: OE expression.
    :return: Region struct expression.
    """
    oe_expr = annotate_region_with_oe(region_expr, oe_expr)
    oe_agg_expr = get_agg_oe_for_region(oe_expr)
    nll_expr = calculate_neg_log_likelihood(oe_expr, oe_agg_expr.oe)
    return hl.struct(
        region=region_expr,
        **oe_agg_expr,
        region_nll=nll_expr,
        nll=nll_expr,
    )


def calculate_oe_neq_1_chisq(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Check for significance that observed/expected values for regions are different from 1.

    Formula is: (obs - exp)^2 / exp.

    :param obs_expr: Observed variant counts.
    :param exp_expr: Expected variant counts.
    :return: Chi-squared value.
    """
    return ((obs_expr - exp_expr) ** 2) / exp_expr


def run_forward(ht, min_exp_mis=MIN_EXP_MIS):
    """
    Run the forward algorithm to find the most intolerant region.

    :param ht: Hail Table with the most intolerant region for each UniProt ID and residue
        index
    :return: Hail Table annotated with the observed and expected values for each residue.
    """
    num_residues = ht.oe.length()
    null_region = hl.range(num_residues)
    null_model = prep_region_struct(null_region, ht.oe)
    ht = ht.select(
        "oe",
        num_residues=num_residues,
        null_model=null_model,
        regions=hl.enumerate(
            ht.min_oe_upper.map(lambda x: x.select("region")).filter(
                lambda x: x.region.length() < num_residues
            )
        ),
        selected=hl.empty_array(null_model.dtype),
        selected_nll=0,
        best_aic=getAIC(null_model, null_model.nll),
        found_best=False,
    )
    ht = ht.explode("regions")
    ht = ht.transmute(idx=ht.regions[0], region=ht.regions[1].region)
    ht = ht.key_by("uniprot_id", "transcript_id", "idx")
    ht = ht.repartition(5000, shuffle=True)
    ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_explode", "ht"))
    ht.describe()
    ht.show(5)
    round_num = 1
    while ht.aggregate(hl.agg.any(hl.is_defined(ht.region))):
        # For each region in regions, update the list of selected by
        # removing the residues in the region from the "catch all remaining"
        # region, and adding the new region to the selected list.
        region_expr = prep_region_struct(ht.region, ht.oe)
        # TODO: Consider adding a checkpoint here or after the next step.
        ht = ht.annotate(_region=region_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.prep", "ht")
        )
        region_expr = ht._region
        updated_null_expr = remove_residues_from_region(ht.null_model, region_expr)
        ht = ht.annotate(_updated_null=updated_null_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.remove", "ht")
        )
        updated_null_expr = ht._updated_null
        updated_null_expr = prep_region_struct(updated_null_expr.region, ht.oe)
        ht = ht.annotate(_updated_null=updated_null_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.prep2", "ht")
        )
        updated_null_expr = ht._updated_null
        region_expr = ht._region
        region_expr = region_expr.annotate(
            null_model=updated_null_expr,
            region_nll=region_expr.nll,
            nll=updated_null_expr.nll + ht.selected_nll + region_expr.nll,
        )
        ht2 = ht.select(exp=region_expr.exp, nll=region_expr.nll)
        ht2 = ht2.filter(hl.is_defined(ht2.nll) & (ht2.exp >= min_exp_mis)).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.scan1", "ht")
        )
        ht2 = (
            ht2.group_by("uniprot_id", "transcript_id")
            .aggregate(
                **hl.agg.fold(
                    hl.missing(hl.tstruct(min_idx=hl.tint, min_nll=hl.tfloat)),
                    lambda accum: (
                        hl.case()
                        .when(
                            hl.is_missing(accum) | (accum.min_nll > ht2.nll),
                            hl.struct(min_idx=ht2.idx, min_nll=ht2.nll),
                        )
                        .when(accum.min_nll <= ht2.nll, accum)
                        .or_missing()
                    ),
                    lambda accum1, accum2: (
                        hl.case()
                        .when(hl.is_missing(accum1), accum2)
                        .when(hl.is_missing(accum2), accum1)
                        .when(accum1.min_nll <= accum2.min_nll, accum1)
                        .default(accum2)
                    ),
                )
            )
            .checkpoint(
                hl.utils.new_temp_file(f"forward_round_{round_num}.scan2", "ht")
            )
        )
        _ht = ht.select(region=region_expr)
        ht2 = ht2.annotate(
            best_region=_ht[ht2.uniprot_id, ht2.transcript_id, ht2.min_idx].region
        )
        ht2 = ht2.checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.scan3", "ht")
        )

        # Get AIC of best candidate model.
        best_region = ht2[ht.uniprot_id, ht.transcript_id].best_region
        region_expr = hl.struct(region=ht.region)

        # Update region list.
        region_expr = remove_residues_from_region(region_expr, best_region).region
        region_expr = hl.or_missing(region_expr.length() > 0, region_expr)

        updated_null_model = best_region.null_model
        candidate_model = ht.selected.append(best_region.drop("null_model"))
        updated_vals = {
            "null_model": updated_null_model,
            "selected": candidate_model,
            "selected_nll": best_region.region_nll + ht.selected_nll,
            "best_aic": getAIC(updated_null_model, 0)
            + getAIC(candidate_model, best_region.nll),
            "region": region_expr,
        }
        curr_vals = {k: ht[k] for k in updated_vals}
        curr_vals["region"] = hl.missing(region_expr.dtype)

        found_best = (
            ht.found_best
            | hl.is_missing(best_region)
            | (updated_vals["best_aic"] >= ht.best_aic)
        )
        curr_vals = hl.struct(**curr_vals)
        updated_vals = hl.struct(**updated_vals)
        ht = ht.annotate(
            **hl.if_else(found_best, curr_vals, updated_vals),
            found_best=found_best,
        )

        ht = ht.filter(hl.is_defined(ht.region) | (ht.idx == 0))
        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}", "ht"))
        round_num += 1

    selected_expr = ht.selected.map(lambda x: x.annotate(is_null=False))
    ht = ht.select(
        selected=add_idx_to_array(
            selected_expr.append(ht.null_model.annotate(is_null=True)), "region_index"
        )
    )
    ht = ht.explode("selected")
    ht = ht.select(**ht.selected)
    ht = ht.annotate(region_length=hl.len(ht.region))
    ht = ht.explode("region")
    chisq_expr = calculate_oe_neq_1_chisq(ht.obs, ht.exp)
    ht = ht.annotate(
        residue_index=ht.region,
        oe_upper=(
            hl.qchisqtail(1 - 0.05 / 2, 2 * (ht.obs + 1), lower_tail=True)
            / (2 * ht.exp)
        ),
        chisq=chisq_expr,
        p_value=hl.pchisqtail(chisq_expr, 1),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index").select(
        "region_index", "obs", "exp", "oe", "oe_upper", "chisq", "p_value", "is_null"
    )

    return ht


def create_missense_viewer_input_ht(
    pos_ht: hl.Table,
    promis3d_ht: hl.Table,
) -> hl.Table:
    """
    Create missense viewer input Hail Table.

    :param ht: Input Hail Table.
    :return: Missense viewer input Hail Table.
    """
    pos_ht = pos_ht.key_by(
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
    pos_ht = pos_ht.key_by("uniprot_id", "enst").cache()

    chisq_expr = calculate_oe_neq_1_chisq(promis3d_ht.obs, promis3d_ht.exp)
    promis3d_ht = promis3d_ht.annotate(
        chisq=chisq_expr, p_value=hl.pchisqtail(chisq_expr, 1)
    )

    # Key by all fields except 'pos' and collect by key into a field named 'pos'.
    promis3d_ht = promis3d_ht.key_by(
        "uniprot_id", "transcript_id", "region_index", "is_null"
    ).collect_by_key("pos")

    # Sort the 'pos' field in ascending order.
    promis3d_ht = promis3d_ht.annotate(
        pos=hl.sorted(promis3d_ht.pos, key=lambda x: x.residue_index)
    )

    # Annotate with 'start' and 'stop' positions for regions by merging adjacent
    # positions.
    promis3d_ht = promis3d_ht.annotate(
        pos=hl.fold(
            lambda i, j: hl.if_else(
                j.residue_index > (i[-1][1] + 1),
                i.append(
                    (
                        j.residue_index,
                        j.residue_index,
                        j.obs,
                        j.exp,
                        j.oe,
                        j.oe_upper,
                        j.chisq,
                        j.p_value,
                    )
                ),
                i[:-1].append(
                    (
                        i[-1][0],
                        j.residue_index,
                        j.obs,
                        j.exp,
                        j.oe,
                        j.oe_upper,
                        j.chisq,
                        j.p_value,
                    )
                ),
            ),
            [
                (
                    promis3d_ht.pos[0].residue_index,
                    promis3d_ht.pos[0].residue_index,
                    promis3d_ht.pos[0].obs,
                    promis3d_ht.pos[0].exp,
                    promis3d_ht.pos[0].oe,
                    promis3d_ht.pos[0].oe_upper,
                    promis3d_ht.pos[0].chisq,
                    promis3d_ht.pos[0].p_value,
                )
            ],
            promis3d_ht.pos[1:],
        )
    )
    promis3d_ht = promis3d_ht.explode("pos")

    # Key by 'gene_id' and transform 'pos' into 'start' and 'stop' fields.
    promis3d_ht = promis3d_ht.key_by(
        "uniprot_id", "transcript_id", "region_index", "is_null"
    )
    promis3d_ht = promis3d_ht.transmute(
        start=promis3d_ht.pos[0],
        stop=promis3d_ht.pos[1],
        obs_mis=promis3d_ht.pos[2],
        exp_mis=promis3d_ht.pos[3],
        obs_exp=promis3d_ht.pos[4],
        oe_upper=promis3d_ht.pos[5],
        chisq=promis3d_ht.pos[6],
        p_value=promis3d_ht.pos[7],
    )

    # Select fields in preferred order and collect by key into a field named 'regions'.
    promis3d_ht = promis3d_ht.collect_by_key("regions")
    promis3d_ht = promis3d_ht.annotate(
        **pos_ht[promis3d_ht.uniprot_id, promis3d_ht.transcript_id]
    )
    promis3d_ht = promis3d_ht.annotate(
        regions=ht.regions.map(
            lambda x: x.select(
                chrom=promis3d_ht.locus_by_aapos[x.start].contig,
                start=hl.if_else(
                    promis3d_ht.locus_by_aapos[x.start].position
                    <= promis3d_ht.locus_by_aapos[x.stop].position,
                    promis3d_ht.locus_by_aapos[x.start].position,
                    promis3d_ht.locus_by_aapos[x.stop].position,
                ),
                stop=hl.if_else(
                    promis3d_ht.locus_by_aapos[x.start].position
                    <= promis3d_ht.locus_by_aapos[x.stop].position,
                    promis3d_ht.locus_by_aapos[x.stop].position + 2,
                    promis3d_ht.locus_by_aapos[x.start].position + 2,
                ),
                aa_start=x.start,
                aa_stop=x.stop,
                obs_mis=x.obs_mis,
                exp_mis=x.exp_mis,
                obs_exp=x.obs_exp,
                oe_upper=x.oe_upper,
                region_index=promis3d_ht.region_index,
                is_null=promis3d_ht.is_null,
                chisq_diff_null=x.chisq,
                p_value=x.p_value,
            )
        )
    )
    promis3d_ht = promis3d_ht.group_by("transcript_id", "uniprot_id").aggregate(
        gnomad_promis3d_constraint=hl.struct(
            has_no_rmc_evidence=False,
            passed_qc=True,
            regions=hl.flatten(hl.agg.collect(promis3d_ht.regions)),
        )
        # regions=hl.flatten(hl.agg.collect(ht.regions))
    )
    promis3d_ht = promis3d_ht.key_by("transcript_id")
    rmc_ht = get_rmc_browser_ht().ht()
    rmc_ht = rmc_ht.annotate(
        regions=rmc_ht.regions.map(
            lambda x: x.select(
                chrom=x.start_coordinate.contig,
                start=hl.if_else(
                    x.start_coordinate.position <= x.stop_coordinate.position,
                    x.start_coordinate.position,
                    x.stop_coordinate.position,
                ),
                stop=hl.if_else(
                    x.start_coordinate.position <= x.stop_coordinate.position,
                    x.stop_coordinate.position + 2,
                    x.start_coordinate.position,
                ),
                aa_start=x.start_aa,
                aa_stop=x.stop_aa,
                obs_mis=x.obs,
                exp_mis=x.exp,
                obs_exp=x.oe,
                chisq_diff_null=x.chisq,
                p_value=x.p,
            )
        )
    )
    promis3d_ht = promis3d_ht.annotate(
        gnomad_regional_missense_constraint=hl.struct(
            has_no_rmc_evidence=False,
            passed_qc=True,
            regions=rmc_ht[promis3d_ht.transcript_id].regions,
        ),
    )

    ht = browser_gene().ht()
    ht = ht.select(
        "interval",
        "gencode_symbol",
        "chrom",
        "strand",
        "start",
        "stop",
        "xstart",
        "xstop",
        "exons",
        "transcripts",
        "reference_genome",
        "canonical_transcript_id",
        "preferred_transcript_id",
        "preferred_transcript_source",
        **ht[ht.canonical_transcript_id],
    )
    ht = ht.repartition(15, shuffle=True)

    return ht


def prioritize_transcripts_and_uniprots(
    residue_ht: hl.Table,
) -> hl.Table:
    """
    Prioritize and label transcript/uniprot combinations for each gene.

    This function:
    1. De-duplicates both tables by key and selects relevant fields.
    2. Annotates residue HT with gene-level info (canonical, mane_select, cds_length, etc.).
    3. Adds a random number to break ties.
    4. Orders transcripts by gene_id, MANE select, canonical, CDS length, and random.
    5. Assigns an index for prioritization.
    6. Aggregates per gene:
       - All transcript/uniprot pairs with their priority.
       - The lowest-index (prioritized) uniprot per transcript.
       - The lowest-index (prioritized) transcript per gene.
    7. Annotates flags indicating for each row if it’s the “one uniprot per transcript” and/or “one transcript per gene”.
    8. Returns a de-nested table keyed by (uniprot_id, transcript_id).

    :param residue_ht: Hail Table keyed by (uniprot_id, transcript_id).
    :return: Annotated and prioritized Hail Table keyed by (uniprot_id, transcript_id).
    """
    # Prepare keys and fields
    ht = residue_ht.key_by("transcript_id", "uniprot_id")
    ht = ht.select("gene_symbol", "gene_id", "canonical", "mane_select", "cds_length")
    ht = ht.distinct()

    # Annotate with a random number for tie-breaking
    ht = ht.annotate(rand_n=hl.rand_unif(0, 1))

    # Order and add index
    ht = ht.order_by(
        "gene_id",
        hl.desc(ht.mane_select),
        hl.desc(ht.canonical),
        hl.desc(ht.cds_length),
        "rand_n",
    )
    ht = ht.add_index()

    # Aggregate per gene for prioritization
    ht = ht.group_by("gene_id").aggregate(
        _all_rows=hl.agg.collect_as_set(
            hl.struct(
                **{
                    k: ht[k]
                    for k in [
                        "transcript_id",
                        "uniprot_id",
                        "gene_symbol",
                        "canonical",
                        "mane_select",
                        "idx",
                    ]
                }
            )
        ),
        one_uniprot_per_transcript=hl.agg.group_by(
            ht.transcript_id, hl.agg.min(ht.idx)
        ),
        one_transcript_per_gene=hl.agg.min(ht.idx),
    )

    # Mark priority flags
    ht = ht.select(
        _all_rows=ht._all_rows.map(
            lambda x: x.annotate(
                one_uniprot_per_transcript=ht.one_uniprot_per_transcript.get(
                    x.transcript_id
                )
                == x.idx,
                one_transcript_per_gene=ht.one_transcript_per_gene == x.idx,
            ).drop("idx")
        )
    )
    ht = ht.explode("_all_rows")

    return ht.select(**ht._all_rows).key_by("transcript_id", "uniprot_id").cache()


def explode_af2_plddt_by_residue(af2_plddt_ht: hl.Table) -> hl.Table:
    """
    Explode AlphaFold2 pLDDT array into per-residue rows.

    This function:
    1. Takes an input Hail Table with a `plddt` array field containing per-residue
       pLDDT scores.
    2. Enumerates the array to associate each score with a 0-based residue index.
    3. Explodes the table to produce one row per residue.
    4. Extracts the `residue_index` and `plddt` values.
    5. Keys the table by (`uniprot_id`, `residue_index`).

    :param af2_plddt_ht: Input Hail Table containing AlphaFold2 pLDDT scores as an array.
    :return: Transformed Hail Table with one row per residue and associated pLDDT score.
    """
    ht = af2_plddt_ht.annotate(plddt=hl.enumerate(af2_plddt_ht.plddt))
    ht = ht.explode("plddt")
    ht = ht.annotate(residue_index=ht.plddt[0], plddt=ht.plddt[1])
    ht = ht.key_by("uniprot_id", "residue_index")

    return ht


def annotate_promis3d_with_af2_metrics(
    promis3D_ht: hl.Table,
    af2_plddt_ht: hl.Table,
    af2_pae_ht: hl.Table,
    af2_dist_ht: hl.Table,
) -> hl.Table:
    """
    Annotate a PROMIS3D Hail Table with per-residue AlphaFold2 metrics (pLDDT, pAE, dist).

    This function:
    1. Explodes the pLDDT array into (residue_index, score) format.
    2. Keys the pAE and distance matrices by (uniprot_id, aa_index).
    3. Filters PROMIS3D rows to coding transcripts (ENST).
    4. Annotates each region residue with pLDDT, pAE, and dist.
    5. Aggregates residue-level and region-level metrics.
    6. Outputs both residue-level and region-level structured annotations.

    :param promis3D_ht: PROMIS3D region Hail Table keyed by (uniprot_id, transcript_id,
        residue_index).
    :param af2_plddt_ht: AlphaFold2 pLDDT Hail Table with array of scores.
    :param af2_pae_ht: AlphaFold2 predicted aligned error matrix Hail Table.
    :param af2_dist_ht: AlphaFold2 distance matrix Hail Table.
    :return: Annotated PROMIS3D Hail Table keyed by (uniprot_id, transcript_id,
        residue_index).
    """
    # Preprocess inputs.
    af2_plddt_ht = explode_af2_plddt_by_residue(af2_plddt_ht)
    af2_pae_ht = af2_pae_ht.key_by("uniprot_id", "aa_index").checkpoint(
        hl.utils.new_temp_file("af2_pae.keyed", "ht")
    )
    af2_dist_ht = af2_dist_ht.key_by("uniprot_id", "aa_index").checkpoint(
        hl.utils.new_temp_file("af2_dist.keyed", "ht")
    )

    # Filter transcripts.
    promis3D_ht = promis3D_ht.filter(promis3D_ht.transcript_id.startswith("ENST"))

    # Annotate with per-residue metrics.
    promis3D_ht = promis3D_ht.annotate(
        plddt=af2_plddt_ht[promis3D_ht.uniprot_id, promis3D_ht.residue_index].plddt,
        pae=af2_pae_ht[promis3D_ht.uniprot_id, promis3D_ht.residue_index].pae,
        dist=af2_dist_ht[promis3D_ht.uniprot_id, promis3D_ht.residue_index].dist_mat,
    ).checkpoint(hl.utils.new_temp_file("promis3D.annotated", "ht"))

    # Group by region.
    promis3D_ht = promis3D_ht.key_by("uniprot_id", "transcript_id", "region_index")
    promis3D_ht = promis3D_ht.collect_by_key("by_residue")

    # Compute per-residue references
    residues_expr = promis3D_ht.by_residue.map(lambda x: x.residue_index)
    by_residue_expr = promis3D_ht.by_residue.map(
        lambda x: x.annotate(
            res=hl.array(hl.set(residues_expr).remove(x.residue_index))
        )
    )
    by_residue_expr = by_residue_expr.map(
        lambda x: x.annotate(
            res=x.res.map(lambda y: hl.abs(x.residue_index - y)),
            pae=x.res.map(lambda y: x.pae[y]),
            dist=x.res.map(lambda y: x.dist[y]),
        )
    )

    # Region-level aggregates
    region_plddt = by_residue_expr.map(lambda x: x.plddt)
    region_res = hl.flatten(by_residue_expr.map(lambda x: x.res))
    region_pae = hl.flatten(by_residue_expr.map(lambda x: x.pae))
    region_dist = hl.flatten(by_residue_expr.map(lambda x: x.dist))

    promis3D_ht = promis3D_ht.annotate(
        by_residue=by_residue_expr,
        region_residues=residues_expr,
        region_aa_dist_stats=region_res.aggregate(lambda x: hl.agg.stats(x)).annotate(
            median=hl.median(region_res)
        ),
        alphafold2_info=hl.struct(
            region_plddt=region_plddt.aggregate(lambda x: hl.agg.stats(x)).annotate(
                median=hl.median(region_plddt)
            ),
            region_pae_stats=region_pae.aggregate(lambda x: hl.agg.stats(x)).annotate(
                median=hl.median(region_pae)
            ),
            region_dist_stats=region_dist.aggregate(lambda x: hl.agg.stats(x)).annotate(
                median=hl.median(region_dist)
            ),
        ),
    )
    promis3D_ht = promis3D_ht.checkpoint(
        hl.utils.new_temp_file("promis3D.region_annotated", "ht")
    )

    # Flatten by residue
    promis3D_ht = promis3D_ht.explode("by_residue")
    promis3D_ht = (
        promis3D_ht.transmute(**promis3D_ht.by_residue)
        .key_by("uniprot_id", "transcript_id", "residue_index")
        .checkpoint(hl.utils.new_temp_file("promis3D.by_residue.exploded", "ht"))
    )

    # Final structured output
    promis3D_ht = promis3D_ht.select(
        residue_level_annotations=hl.struct(
            residue_to_region_aa_dist_stats=promis3D_ht.res.aggregate(
                lambda x: hl.agg.stats(x)
            ).annotate(median=hl.median(promis3D_ht.res)),
            alphafold2_info=hl.struct(
                residue_plddt=promis3D_ht.plddt,
                residue_to_region_pae_stats=promis3D_ht.pae.aggregate(
                    lambda x: hl.agg.stats(x)
                ).annotate(median=hl.median(promis3D_ht.pae)),
                residue_to_region_dist_stats=promis3D_ht.dist.aggregate(
                    lambda x: hl.agg.stats(x)
                ).annotate(median=hl.median(promis3D_ht.dist)),
            ),
        ),
        region_level_annotations=hl.struct(
            region_index=promis3D_ht.region_index,
            region_residues=promis3D_ht.region_residues,
            region_length=promis3D_ht.region_length,
            obs=promis3D_ht.obs,
            exp=promis3D_ht.exp,
            oe=promis3D_ht.oe,
            oe_upper=promis3D_ht.oe_upper,
            oe_ci=oe_confidence_interval(promis3D_ht.obs, promis3D_ht.exp),
            chisq=promis3D_ht.chisq,
            p_value=promis3D_ht.p_value,
            is_null=promis3D_ht.is_null,
            region_aa_dist_stats=promis3D_ht.region_aa_dist_stats,
            alphafold2_info=promis3D_ht.alphafold2_info,
        ),
    )

    return promis3D_ht


def generate_all_possible_snvs_from_gencode_positions(
    transcripts_ht: hl.Table,
    translations_ht: hl.Table,
    gencode_gtf_ht: hl.Table,
    gencode_translations_matched_ht: hl.Table,
) -> hl.Table:
    """
    Generate all possible single nucleotide variants (SNVs) from the GENCODE positions Hail Table.

    :param transcripts_ht: GENCODE transcripts Hail Table.
    :param translations_ht: GENCODE translations Hail Table.
    :param gencode_gtf_ht: GENCODE GTF Hail Table.
    :param gencode_translations_matched_ht: GENCODE translations matched Hail Table.
    :return: Hail Table with all possible SNVs at each GENCODE position.
    """
    gencode_translations_matched_ht = gencode_translations_matched_ht.group_by(
        "enst"
    ).aggregate(
        uniprot_id=hl.agg.collect_as_set(gencode_translations_matched_ht.uniprot_id)
    )
    translations_ht = translations_ht.annotate(
        uniprot_id=hl.or_else(
            gencode_translations_matched_ht[translations_ht.enst].uniprot_id, {"None"}
        )
    )
    translations_ht = translations_ht.explode("uniprot_id")

    ht = get_gencode_positions(
        transcripts_ht,
        translations_ht,
        gencode_gtf_ht,
        no_filter=True,
    )
    ht = ht.annotate(
        uniprot_id=hl.if_else(
            ht.cds_len_mismatch | ht.cds_len_not_div_by_3, "None", ht.uniprot_id
        )
    )
    ht = ht.filter(ht.locus.contig != "chrM")

    nucleotides = hl.set({"A", "T", "C", "G"})
    has_aa_info_expr = ~ht.cds_len_mismatch & ~ht.cds_len_not_div_by_3
    ht = ht.annotate(
        aminoacid_ref=hl.or_missing(has_aa_info_expr, ht.sequence[ht.aapos]),
        alleles=nucleotides.remove(ht.ref).map(lambda alt: [ht.ref, alt]),
        residue_index=hl.or_missing(has_aa_info_expr, ht.aapos),
        aminoacid_length=hl.int(ht.aalength),
        cds_length=hl.int(ht.cds_len),
        transcript_id=ht.enst,
        gene_id=ht.ensg,
        gene_symbol=ht.gene,
    )
    ht = ht.explode("alleles")
    ht = ht.key_by("locus", "alleles", "transcript_id", "uniprot_id", "gene_id")
    ht = ht.select(
        "gene_symbol",
        "strand",
        "cds_length",
        "aminoacid_length",
        "residue_index",
        "aminoacid_ref",
        "cds_len_mismatch",
        "cds_len_not_div_by_3",
    )

    return ht


def make_temp_annotation_ht(
    base_ht: hl.Table,
    annotation_ht: hl.Table,
    keys: List[str] = ["locus", "alleles"],
    temp_path_prefix: str = "tmp_annotation_ht",
    annotation_name: Optional[str] = None,
) -> hl.Table:
    """
    Make a temporary Hail Table with annotations from another Hail Table.

    :param base_ht: Base Hail Table to annotate.
    :param annotation_ht: Annotation Hail Table to index.
    :param keys: List of keys to index the annotation Hail Table with.
    :param temp_path_prefix: Prefix for the temporary file path.
    :param annotation_name: Name of the annotation to annotate the base Hail Table with.
    :return: Annotated Hail Table.
    """
    base_ht = base_ht.select(*[k for k in keys if k not in base_ht.key])
    keys_expr = [base_ht[k] for k in keys]
    fields = [f for f in annotation_ht.row_value if f not in base_ht.key]
    annotation_expr = annotation_ht.select(*fields)
    annotation_expr = annotation_expr.index(*keys_expr)
    if annotation_name is not None:
        annotation_expr = {annotation_name: annotation_expr}

    base_ht = base_ht.annotate(**annotation_expr)
    base_ht = base_ht.checkpoint(hl.utils.new_temp_file(temp_path_prefix, "ht"))

    return base_ht


def annotate_snvs_with_variant_level_data(ht: hl.Table) -> hl.Table:
    """
    Annotate a per-SNV Hail Table with variant-level annotations from multiple sources.

    Adds the following annotations:

        - context
        - gnomad_site
        - revel
        - cadd
        - phylop
        - genetics_gym
        - autism
        - dd_denovo
        - dd_denovo_no_transcript
        - gnomad_de_novo
        - clinvar
        - pext_base
        - pext_annotation
        - rmc

    :param ht: Input Hail Table.
    :return: Annotated Hail Table.
    """
    annotation_hts = {
        n: make_temp_annotation_ht(
            ht,
            (
                c["ht"].ht()
                if "custom_select" not in c
                else c["custom_select"](c["ht"].ht())
            ),
            keys=c["keys"],
            temp_path_prefix=n,
            annotation_name=c["annotation_name"],
        )[ht.key]
        for n, c in VARIANT_LEVEL_ANNOTATION_CONFIG.items()
    }
    annotation_expr = hl.struct()
    for t in annotation_hts.values():
        annotation_expr = annotation_expr.annotate(**t)

    ht = ht.annotate(variant_level_annotations=annotation_expr).checkpoint(
        hl.utils.new_temp_file("snvs_with_variant_level_data", "ht")
    )
    pext_annotation_ht = process_pext_annotation_ht(pext("annotation_level").ht())

    var_expr = ht.variant_level_annotations.rename({"biotype": "transcript_biotype"})
    var_update_expr = var_expr.annotate(
        dd_denovo=var_expr.dd_denovo.annotate(
            **get_kaplanis_sig_gene_annotations(ht.gene_symbol)
        ),
        dd_denovo_no_transcript_match=(
            var_expr.dd_denovo_no_transcript_match.annotate(
                **get_kaplanis_sig_gene_annotations(ht.gene_symbol)
            ),
        ),
        annotation_level_pext=pext_annotation_ht[
            ht.locus, ht.alleles, ht.gene_id, var_expr.most_severe_consequence
        ],
    )
    rearrange_fields = BASE_LEVEL_ANNOTATION_FIELDS + ["residue_alt"]
    var_update_expr = var_update_expr.drop(*rearrange_fields)
    ht = ht.annotate(
        **{k: var_expr[k] for k in rearrange_fields},
        variant_level_annotations=var_update_expr,
    )

    return ht


def combine_residue_level_annotations(
    ht: hl.Table,
    promis3d_ht: hl.Table,
) -> hl.Table:
    """
    Combine residue-level annotations by joining PROMIS3D regions with COSMIS (multiple sources) and InterPro data.

    Adds the following annotations:

        - interpro
        - mtr3d
        - cosmis_alphafold
        - cosmis_pdb
        - cosmis_swiss_model
        - promis3d

    :param ht: Input Hail Table.
    :param promis3d_ht: PROMIS3D Hail Table.
    :return: Annotated Hail Table.
    """
    annotation_hts = {
        **RESIDUE_LEVEL_ANNOTATION_CONFIG,
        "promis3d": {
            "ht": promis3d_ht,
            "keys": ["transcript_id", "uniprot_id", "residue_index"],
            "annotation_name": "promis3d",
        },
    }
    annotation_hts = {
        n: make_temp_annotation_ht(
            ht,
            (
                c["ht"]
                if isinstance(c["ht"], hl.Table)
                else (
                    c["ht"].ht()
                    if "custom_select" not in c
                    else c["custom_select"](c["ht"].ht())
                )
            ),
            keys=c["keys"],
            temp_path_prefix=n,
            annotation_name=c["annotation_name"],
        )[ht.key]
        for n, c in annotation_hts.items()
    }
    annotation_expr = hl.struct()
    for t in annotation_hts.values():
        annotation_expr = annotation_expr.annotate(**t)

    ht = ht.annotate(**annotation_expr)
    ht = ht.transmute(
        cosmis=hl.struct(
            alphafold=ht.row_value.cosmis_alphafold,
            pdb=ht.row_value.cosmis_pdb,
            swiss_model=ht.row_value.cosmis_swiss_model,
        )
    )

    return ht


def create_per_snv_combined_ht(
    ht: hl.Table,
    promis3d_ht: hl.Table,
    af2_plddt_ht: hl.Table,
    af2_pae_ht: hl.Table,
    af2_dist_ht: hl.Table,
) -> hl.Table:
    """
    Create a fully annotated per-SNV Hail Table with structured variant-, residue-, and gene-level annotations.

    :param ht: All possible SNVs Hail Table.
    :param promis3d_ht: PROMIS3D Hail Table.
    :param af2_plddt_ht: AlphaFold2 pLDDT Hail Table.
    :param af2_pae_ht: AlphaFold2 PAE Hail Table.
    :param af2_dist_ht: AlphaFold2 distance matrix Hail Table.
    :param partition_intervals: Partition intervals to read annotation Hail Tables with.
    :return: Annotated and checkpointed Hail Table.
    """
    hl._set_flags(use_new_shuffle="1")
    ht = annotate_snvs_with_variant_level_data(ht).naive_coalesce(5000).cache()
    hl._set_flags(use_new_shuffle=None)

    ht = ht.annotate(
        residue_ref=ht.aminoacid_ref,
        residue_alt=ht.variant_level_annotations.residue_alt,
    )

    base_residue_ht = (
        ht.key_by("transcript_id", "uniprot_id", "residue_index")
        .select(
            "gene_id",
            "gene_symbol",
            "canonical",
            "mane_select",
            "cds_length",
            "residue_ref",
            "residue_alt",
        )
        .distinct()
    ).cache()
    residue_ht = combine_residue_level_annotations(
        base_residue_ht,
        annotate_promis3d_with_af2_metrics(
            promis3d_ht, af2_plddt_ht, af2_pae_ht, af2_dist_ht
        ).cache(),
    ).cache()
    select_uniprot_transcript_ht = prioritize_transcripts_and_uniprots(
        base_residue_ht
    ).select("one_uniprot_per_transcript", "one_transcript_per_gene")

    hl._set_flags(use_new_shuffle="1")
    residue_ht = make_temp_annotation_ht(
        ht,
        residue_ht,
        keys=["transcript_id", "uniprot_id", "residue_index"],
        temp_path_prefix="residue",
    ).drop("residue_index")
    select_uniprot_transcript_ht = make_temp_annotation_ht(
        ht,
        select_uniprot_transcript_ht,
        keys=["transcript_id", "uniprot_id"],
        temp_path_prefix="select_uniprot_transcript",
    )
    gene_constraint_ht = make_temp_annotation_ht(
        ht,
        get_temp_processed_constraint_ht().ht(),
        keys=["transcript_id"],
        temp_path_prefix="gene_constraint",
    )
    hl._set_flags(use_new_shuffle=None)

    hi_expr = hl.case()
    for n, g in HI_GENE_CATEGORIES.items():
        hi_expr = hi_expr.when(hl.set(g).contains(ht.gene_symbol), n)
    hi_expr = hi_expr.or_missing()

    # Structure final output.
    ht = ht.select(
        "gene_symbol",
        "canonical",
        "mane_select",
        "transcript_biotype",
        "most_severe_consequence",
        "cds_len_mismatch",
        "cds_len_not_div_by_3",
        is_phaplo_gene=hl.set(get_phaplo().he()).contains(ht.gene_symbol),
        is_ptriplo_gene=hl.set(get_ptriplo().he()).contains(ht.gene_symbol),
        is_hi_gene=hl.set(HI_GENES).contains(ht.gene_symbol),
        hi_gene_category=hi_expr,
        **select_uniprot_transcript_ht[ht.key],
        variant_level_annotations=ht.variant_level_annotations,
        residue_level_annotations=hl.struct(**residue_ht[ht.key]),
        gene_level_annotations=hl.struct(
            strand=ht.strand,
            cds_length=ht.cds_length,
            cds_len_mismatch=ht.cds_len_mismatch,
            cds_len_not_div_by_3=ht.cds_len_not_div_by_3,
            aminoacid_length=ht.aminoacid_length,
            **gene_constraint_ht[ht.key],
        ),
    )

    return ht


def create_per_residue_ht_from_snv_ht(per_snv_ht: hl.Table) -> hl.Table:
    """
    Create a per-residue Hail Table from a fully annotated per-SNV Hail Table.

    This function:
    1. Extracts residue-relevant fields and drops alleles.
    2. Deduplicates rows at the residue level.
    3. Aggregates mean coverage and RMC sets per residue.
    4. Extracts flattened annotations from residue and gene level.

    :param per_snv_ht: Annotated per-SNV Hail Table from `create_per_snv_combined_ht`.
    :return: Final aggregated per-residue Hail Table.
    """
    keep_fields = [
        *BASE_LEVEL_ANNOTATION_FIELDS,
        "is_phaplo_gene",
        "is_ptriplo_gene",
        "is_hi_gene",
        "hi_gene_category",
        "cds_len_mismatch",
        "cds_len_not_div_by_3",
        "one_uniprot_per_transcript",
        "one_transcript_per_gene",
    ]

    # Extract and deduplicate.
    ht = (
        per_snv_ht.select(
            *keep_fields,
            residue_index=per_snv_ht.residue_level_annotations.residue_index,
            exomes_coverage=per_snv_ht.variant_level_annotations.exomes_coverage,
            rmc=per_snv_ht.variant_level_annotations.rmc,
            residue_level_annotations=per_snv_ht.residue_level_annotations,
            gene_level_annotations=per_snv_ht.gene_level_annotations,
        )
        .key_by("locus", "transcript_id", "uniprot_id")
        .drop("alleles")
    )

    ht = ht.distinct()
    ht = ht.checkpoint(hl.utils.new_temp_file("per_residue_dedup", "ht"))

    # Group and aggregate.
    ht = ht.group_by("transcript_id", "uniprot_id", "residue_index").aggregate(
        **{k: hl.agg.take(ht[k], 1)[0] for k in keep_fields},
        residue_mean_exomes_coverage=hl.struct(
            mean=hl.agg.mean(ht.exomes_coverage.mean),
            median_approx=hl.agg.mean(ht.exomes_coverage.median_approx),
            AN=hl.agg.mean(ht.exomes_coverage.AN),
            percent_AN=hl.agg.mean(ht.exomes_coverage.percent_AN),
        ),
        rmc=hl.agg.collect_as_set(ht.rmc),
        residue_level_annotations=hl.agg.take(ht.residue_level_annotations, 1)[0],
        gene_level_annotations=hl.agg.take(ht.gene_level_annotations, 1)[0],
    )
    ht = ht.checkpoint(hl.utils.new_temp_file("per_residue_agg", "ht"))

    # Flatten and finalize.
    ht = ht.select(
        *keep_fields,
        residue_ref=ht.residue_level_annotations.residue_ref,
        residue_mean_exomes_coverage=ht.residue_mean_exomes_coverage,
        interpro=ht.residue_level_annotations.interpro,
        rmc_regions=ht.rmc,
        cosmis=ht.residue_level_annotations.cosmis,
        promis3d=ht.residue_level_annotations.promis3d,
        gene_level_annotations=ht.gene_level_annotations,
    )

    return ht.naive_coalesce(2000)


def create_per_promis3d_region_ht_from_residue_ht(ht: hl.Table) -> hl.Table:
    """
    Create a PROMIS3D region-level Hail Table from a per-residue annotated Hail Table.

    This function:
    1. Extracts region_index from PROMIS3D residue annotations.
    2. Groups rows by (transcript_id, uniprot_id, region_index).
    3. Aggregates exomes coverage, collects RMC regions, computes per-region pLDDT stats.
    4. Annotates PROMIS3D region-level alphaFold2 info with region-level pLDDT stats.

    :param ht: Hail Table with residue-level annotations including PROMIS3D and coverage.
    :return: Aggregated region-level Hail Table.
    """
    keep_fields = [
        "gene_id",
        *BASE_LEVEL_ANNOTATION_FIELDS,
        "is_phaplo_gene",
        "is_ptriplo_gene",
        "is_hi_gene",
        "hi_gene_category",
        "cds_len_mismatch",
        "cds_len_not_div_by_3",
        "one_uniprot_per_transcript",
        "one_transcript_per_gene",
    ]

    ht = ht.annotate(region_index=ht.promis3d.region_level_annotations.region_index)

    ht = ht.group_by("transcript_id", "uniprot_id", "region_index").aggregate(
        **{k: hl.agg.take(ht[k], 1)[0] for k in keep_fields},
        region_mean_exomes_coverage=hl.struct(
            mean=hl.agg.mean(ht.residue_mean_exomes_coverage.mean),
            median_approx=hl.agg.mean(ht.residue_mean_exomes_coverage.median_approx),
            AN=hl.agg.mean(ht.residue_mean_exomes_coverage.AN),
            percent_AN=hl.agg.mean(ht.residue_mean_exomes_coverage.percent_AN),
        ),
        rmc_regions=hl.agg.explode(
            lambda x: hl.agg.collect_as_set(x), ht.rmc_regions
        ).filter(lambda x: hl.is_defined(x)),
        promis3d=hl.agg.take(ht.promis3d.region_level_annotations, 1)[0].annotate(
            region_plddt_stats=hl.agg.stats(
                ht.promis3d.residue_level_annotations.alphafold2_info.residue_plddt
            ).annotate(
                median=hl.median(
                    hl.agg.collect(
                        ht.promis3d.residue_level_annotations.alphafold2_info.residue_plddt
                    )
                )
            )
        ),
        gene_level_annotations=hl.agg.take(ht.gene_level_annotations, 1)[0],
    )

    # Move pLDDT stats into alphaFold2 info.
    ht = ht.annotate(
        promis3d=ht.promis3d.annotate(
            alphafold2_info=ht.promis3d.alphafold2_info.annotate(
                region_plddt_stats=ht.promis3d.region_plddt_stats
            )
        ).drop("region_plddt_stats")
    )

    return ht.naive_coalesce(100)
