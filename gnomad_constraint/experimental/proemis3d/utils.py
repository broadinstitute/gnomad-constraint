"""Script with utility functions for the Proemis3D pipeline."""

import io
import json
import logging
import os
from typing import Dict, Iterator, List, Optional, Union

import hail as hl
import numpy as np
import pandas as pd

# from Bio.PDB import PPBuilder
# from Bio.PDB.MMCIFParser import MMCIFParser
# from Bio.PDB.Polypeptide import is_aa
from gnomad.resources.grch38.gnomad import browser_gene, browser_variant, pext
from gnomad.utils.constraint import oe_confidence_interval
from gnomad.utils.filtering import filter_gencode_ht
from gnomad.utils.reference_genome import get_reference_genome
from hail.utils.misc import divide_null
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.functions import col, explode, pandas_udf, rtrim, split
from pyspark.sql.types import StringType, StructField, StructType

from gnomad_constraint.experimental.proemis3d.constants import (
    HI_GENE_CATEGORIES,
    HI_GENES,
    MIN_EXP_MIS,
)
from gnomad_constraint.experimental.proemis3d.data_import import (
    get_kaplanis_sig_gene_annotations,
    process_gnomad_de_novo_ht,
    process_gnomad_site_ht,
    process_pext_annotation_ht,
    process_pext_base_ht,
)
from gnomad_constraint.experimental.proemis3d.resources import (
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
logger = logging.getLogger("proemis3d_utils")
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


# def get_structure_peptide(structure) -> str:
#    """
#    Get the sequence from a structure.
#
#    :param structure: Structure object.
#    :return: Sequence as a string.
#    """
#    ppb = PPBuilder()

#    # Return the sequence as a string.
#    return "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(structure)])


# def get_structure_dist_matrix(structure: MMCIFParser) -> np.ndarray:
#    """
#    Calculate the "calpha" distance matrix from a structure.

#    :param structure: Structure object.
#    :return: Distance matrix as a NumPy array.
#    """
#    calpha_atoms = []
#    for model in structure:
#        for chain in model:
#            for residue in chain:
#                if is_aa(residue, standard=True) and "CA" in residue:
#                    calpha_atoms.append(residue["CA"].get_coord())

#    def _calc_dist_matrix(calpha_atoms: List[np.ndarray]) -> np.ndarray:
#        """
#        Calculate the pairwise distance matrix between Calpha atoms.
#
#        :param calpha_atoms: List of Calpha atoms.
#        :return: Distance matrix as a NumPy array.
#        """
#        num_atoms = len(calpha_atoms)
#        dist_matrix = np.zeros((num_atoms, num_atoms))
#        for i, atom1 in enumerate(calpha_atoms):
#            for j, atom2 in enumerate(calpha_atoms):
#                dist_matrix[i, j] = np.linalg.norm(atom1 - atom2)
#
#        return dist_matrix
#
#    # Calculate the distance matrix
#    return _calc_dist_matrix(calpha_atoms)


# def process_af2_mmcif(
#    uniprot_id: str,
#    mmcif_content: str,
#    distance_matrix: bool = False,
# ) -> Union[str, np.ndarray, List[float]]:
#    """
#    Process AlphaFold2 MMCIF content.
#
#    :param uniprot_id: UniProt ID.
#    :param mmcif_content: MMCIF content as a string.
#    :param distance_matrix: Whether to return the distance matrix. Default is False.
#    :return: Sequence or distance matrix.
#    """
#    parser = MMCIFParser(QUIET=True)
#    structure = parser.get_structure(uniprot_id, io.StringIO(mmcif_content))

#    if distance_matrix:
#        return get_structure_dist_matrix(structure)
#    else:
#        return get_structure_peptide(structure)


# def process_af2_file_by_mode(
#    uniprot_id: str,
#    file_content: str,
#    mode: str,
# ) -> Union[str, np.ndarray, List[float]]:
#    """
#    Dispatcher to handle different AF2 modes based on filename suffix and mode.

#    :param uniprot_id: UniProt ID.
#    :param file_content: File content.
#    :param mode: Mode for processing files. Options are 'sequence', 'distance_matrix',
#        'plddt', or 'pae'. Default is 'sequence'.
#    :return: Sequence or distance matrix.
#    """
#    if mode in {"sequence", "distance_matrix"}:
#        return process_af2_mmcif(
#            uniprot_id, file_content, distance_matrix=(mode == "distance_matrix")
#        )
#    if mode == "plddt":
#        return get_plddt_from_confidence_json(file_content)

#    if mode == "pae":
#        return get_pae_from_json(file_content)

#    raise ValueError(f"Unsupported mode: {mode}")


# def process_af2_structures(
#    gcs_bucket_glob: str,
#    mode: str = "sequence",
# ) -> hl.Table:
#    """
#    Process AlphaFold2 structures from a GCS bucket.
#
#    .. note::
#
#        All files in the bucket must be in CIF format with a '.cif.gz' extension.

#    :param gcs_bucket_glob: GCS bucket glob pattern.
#    :param mode: Mode for processing files. Options are 'sequence', 'distance_matrix',
#        'plddt', or 'pae'. Default is 'sequence'.
#    :return: Hail Table with UniProt IDs and sequences or distance matrices.
#    """
#    # Get Spark session for file distribution and processing.
#    spark = hl.utils.java.Env.spark_session()
#    spark.conf.set(
#        "spark.sql.execution.arrow.maxRecordsPerBatch",
#        1000,
#    )

#    # Define schema for loading the files.
#    schema = StructType(
#        [
#            StructField("file_content", StringType(), True),
#            StructField("af2_file", StringType(), True),
#        ]
#    )

#    # Use Spark to read files in parallel.
#    # This reads the entire content of each file as a (filename, content) pair.
#    af2_files_df = (
#        spark.read.format("text")
#        .load(gcs_bucket_glob, schema=schema, wholetext=True)
#        .withColumn("af2_file", col("_metadata.file_path"))
#    )
#    if mode == "distance_matrix":
#        col_name = "dist_mat"
#        rtype = "array<array<float>>"
#    elif mode == "pae":
#        col_name = "pae"
#        rtype = "array<array<int>>"
#    elif mode == "plddt":
#        col_name = "plddt"
#        rtype = "array<float>"
#    else:
#        col_name = "sequence"
#        rtype = "string"

#    @pandas_udf(f"uniprot_id string, {col_name} {rtype}")
#    def process_file(
#        file_path_series: pd.Series, file_content_series: pd.Series
#    ) -> pd.DataFrame:
#        """
#        Process a list of files in parallel using a Pandas UDF.

#        :param file_path_series: File paths.
#        :param file_content_series: File contents.
#        :return: Pandas DataFrame with UniProt IDs and sequences.
#        """
#        result = []
#        for file_path, file_content in zip(file_path_series, file_content_series):
#            # Extract UniProt ID from the file path.
#            uniprot_id = os.path.basename(file_path).split("-")[1]

#            # Process the file content.
#            af2_data = process_af2_file_by_mode(uniprot_id, file_content, mode=mode)
#            result.append((uniprot_id, af2_data))

#        return pd.DataFrame(result, columns=["uniprot_id", col_name])

#    if mode in {"distance_matrix", "pae"}:
#        from pyspark.sql.functions import posexplode

#        result_df = result_df.select(
#            "af2_file",
#            "uniprot_id",
#            posexplode(col(col_name)).alias("aa_index", col_name),
#        )

#    # Convert the Spark DataFrame to a Hail Table.
#    key = ["af2_file", "uniprot_id"]
#    if mode in {"distance_matrix", "pae"}:
#        key.append("aa_index")
#
#    ht = hl.Table.from_spark(result_df, key=key)
#
#    return ht


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
    gencode_gtf_ht = gencode_gtf_ht.annotate(chrom=gencode_gtf_ht.interval.start.contig)
    gencode_gtf_ht = gencode_gtf_ht.group_by(
        "transcript_id", "strand", "chrom"
    ).aggregate(intervals=hl.agg.collect(gencode_gtf_ht.interval))

    # Get CDS positions and lengths for each transcript in the GTF data.
    positions = gencode_gtf_ht.intervals.flatmap(
        lambda x: hl.range(x.start.position, x.end.position + 1)
    )
    gencode_gtf_ht = gencode_gtf_ht.transmute(
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
    # Handle empty arrays
    element_type = oe_expr.dtype.element_type
    oe_expr = (
        hl.case()
        .when(hl.len(oe_expr) == 0, hl.empty_array(element_type))
        .default(
            hl.array_scan(
                lambda i, j: j.annotate(
                    obs=hl.or_else(i.obs, 0) + hl.or_else(j.obs, 0),
                    exp=hl.or_else(i.exp, 0.0) + hl.or_else(j.exp, 0.0),
                ),
                oe_expr[0],
                oe_expr[1:],
            )
        )
    )

    return oe_expr


def gamma_upper_ci(
    obs: hl.expr.Int32Expression,
    exp: hl.expr.Float64Expression,
    alpha: float = 0.05,
) -> hl.expr.Float64Expression:
    """
    Calculate the upper bound of the OE confidence interval using the Gamma distribution.

    This function uses the built-in qgamma function from the custom Hail wheel.

    :param obs: Observed count
    :param exp: Expected count
    :param alpha: Significance level for the confidence interval. Default is 0.05.
    :return: Upper bound of the OE confidence interval
    """
    # Calculate shape and scale parameters for Gamma distribution
    shape = obs + hl.literal(1.0)
    scale = divide_null(
        hl.literal(1.0), exp
    )  # Use divide_null to handle division by zero
    p = hl.literal(1.0 - alpha)

    # Use the built-in qgamma function from the custom Hail wheel
    # divide_null will return null if exp is 0, making the result null as well
    return hl.qgamma(p, shape, scale)


def chisq_upper_ci(
    obs: hl.expr.Int32Expression,
    exp: hl.expr.Float64Expression,
    alpha: float = 0.05,
) -> hl.expr.Float64Expression:
    """
    Calculate the upper bound of the OE confidence interval using the chi-squared
    distribution.

    :param obs: Observed count.
    :param exp: Expected count.
    :param alpha: Significance level for the confidence interval. Default is 0.05.
    :return: Upper bound of the OE confidence interval.
    """
    return hl.qchisqtail(1 - alpha / 2, 2 * (obs + 1), lower_tail=True) / (2 * exp)


def calculate_oe_upper(oe_expr, alpha=0.05, oe_upper_method: str = "gamma"):
    """
    Calculate the upper bound of the OE confidence interval.

    :param oe_expr: Array expression with observed and expected values.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :return: Array expression with upper bound of the OE confidence interval.
    """
    # Calculate upper bound of oe confidence interval.
    if oe_upper_method not in ["gamma", "chisq"]:
        raise ValueError(f"Invalid OE upper method: {oe_upper_method}")

    oe_upper_func = gamma_upper_ci if oe_upper_method == "gamma" else chisq_upper_ci
    oe_upper_expr = oe_expr.map(
        lambda x: x.annotate(
            oe=divide_null(x.obs, x.exp),
            oe_upper=oe_upper_func(x.obs, x.exp, alpha),
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

    filtered_oe_expr = filtered_oe_expr.filter(lambda x: hl.is_defined(x.oe_upper))
    # TODO: minimize on oe instead of oe_upper?
    min_oe_upper_expr = hl.sorted(filtered_oe_expr, key=lambda x: x.oe_upper).first()
    dist_index_expr = min_oe_upper_expr.dist_index
    oe_expr = hl.or_missing(
        hl.is_defined(filtered_oe_expr), oe_expr[: dist_index_expr + 1]
    )
    min_oe_upper_expr = min_oe_upper_expr.drop("dist_index")
    min_oe_upper_expr = min_oe_upper_expr.annotate(
        region=oe_expr.map(lambda x: x.residue_index),
        dists=oe_expr.map(lambda x: x.dist),
    )

    return min_oe_upper_expr


def debug_print_pae_matrix_for_region(
    dist_mat_expr: hl.expr.ArrayExpression,
    center_residue_index_expr: hl.expr.Int32Expression,
    max_pae: float,
    title: str = "PAE matrix for region filtering",
    min_plddt: Optional[float] = None,
    plddt_cutoff_method: Optional[str] = None,
) -> Dict[int, str]:
    """
    Generate debug output showing PAE values between each residue and all subsequent residues
    in the sorted-by-distance array. This helps visualize how filter_on_pairwise_pae_in_region works.
    Returns a dictionary where keys are center residue indices and values are output strings.

    :param dist_mat_expr: Array expression with distance matrix entries (sorted by distance)
    :param center_residue_index_expr: Center residue index expression
    :param max_pae: Maximum allowed PAE value
    :param title: Title for the debug output
    :return: Dictionary mapping center residue indices to their debug output strings
    """
    # Dictionary to store output per center residue
    output_dict = {}

    # Get first key for debugging
    ht = dist_mat_expr._indices.source
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    # Filter and collect data
    _ht_debug = ht.annotate(
        dist_mat=dist_mat_expr, center_idx=center_residue_index_expr
    )
    _ht_debug = _ht_debug.filter(
        (_ht_debug.uniprot_id == uniprot_id) & (_ht_debug.enst == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return {}

    # Collect data for all center residues
    debug_data = _ht_debug.select("dist_mat", "center_idx").collect()

    if not debug_data:
        return {}

    # Group by center residue to show matrix for each
    center_groups = {}
    for row in debug_data:
        center_idx = row.center_idx
        if center_idx not in center_groups:
            center_groups[center_idx] = []
        center_groups[center_idx].append(row)

    # Show PAE matrix for each center residue
    for center_idx in sorted(center_groups.keys()):
        rows = center_groups[center_idx]
        # Use first row for this center (they should all have the same dist_mat)
        row = rows[0]
        dist_mat = row.dist_mat

        # Build output for this center residue
        center_output = f"    {BOLD}=== {title} ==={RESET}\n\n"
        center_output += f"        {BOLD}PAE matrix (row = from residue, column = to residue):{RESET}\n\n"
        center_output += f"        {BOLD}Color legend:{RESET}\n"
        center_output += f"            {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}\n"

        if min_plddt is not None:
            center_output += f"            {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}\n"

        center_output += f"            {BOLD}Filtered residues:{RESET} {ORANGE}Orange{RESET} (excluded from region)\n"

        center_output += f"\n"

        # Check if PAE data is available in dist_mat
        # Note: residue_i.pae might be a float (PAE from center) or an array (full PAE row)
        # Also check for pae_array field which contains the full PAE array for
        # each residue
        has_pae_data = (
            dist_mat
            and len(dist_mat) > 0
            and (hasattr(dist_mat[0], "pae") or hasattr(dist_mat[0], "pae_array"))
        )
        pae_is_array = False
        has_pae_array_field = False
        if has_pae_data:
            # Check if pae_array field exists (preferred for pairwise PAE)
            if hasattr(dist_mat[0], "pae_array") and dist_mat[0].pae_array is not None:
                try:
                    len(dist_mat[0].pae_array)
                    has_pae_array_field = True
                    pae_is_array = True
                except (TypeError, AttributeError):
                    pass
            # Fall back to checking pae field
            if (
                not has_pae_array_field
                and hasattr(dist_mat[0], "pae")
                and dist_mat[0].pae is not None
            ):
                try:
                    len(dist_mat[0].pae)
                    pae_is_array = True
                except (TypeError, AttributeError):
                    pae_is_array = False

        # Check if pLDDT data is available
        has_plddt_data = (
            dist_mat and len(dist_mat) > 0 and hasattr(dist_mat[0], "plddt")
        )

        # Get all residue indices for column headers
        all_residue_indices = [entry.residue_index for entry in dist_mat]

        # Determine which residues would be filtered by pLDDT (if pLDDT filtering
        # is used)
        plddt_filtered_residues = set()
        # Also track low pLDDT residues for exclude_low_plddt_from_stats (they
        # remain but are excluded from PAE checks)
        low_plddt_residues = set()
        if min_plddt is not None and has_plddt_data and plddt_cutoff_method:
            if plddt_cutoff_method == "truncate_at_first_low_plddt":
                # Find first residue with low pLDDT and mark all after it
                found_first_low = False
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if found_first_low:
                        plddt_filtered_residues.add(residue_idx)
                    elif (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        found_first_low = True
                        # Don't add the first low one, but mark all subsequent ones
            elif plddt_cutoff_method == "remove_low_plddt_residues":
                # Mark residues with low pLDDT
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        plddt_filtered_residues.add(residue_idx)
            elif plddt_cutoff_method == "exclude_low_plddt_from_stats":
                # Mark residues with low pLDDT (they remain in region but are excluded
                # from PAE checks)
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        low_plddt_residues.add(residue_idx)
            # For "mask_low_confidence_plddt", residues are not filtered from the region

        # Simulate the filtering logic to determine which residues would be filtered by PAE
        # This mimics the scan operation in filter_on_pairwise_pae_in_region
        # The logic: add residue if (pae_array is None) OR (NOT any(PAE > max_pae to residues in region))
        # Note: pLDDT filtering happens FIRST in the actual code, so we simulate pLDDT filtering first,
        # then PAE filtering on the remaining residues
        pae_filtered_residues = set()
        region = []

        if has_pae_array_field:
            # First, apply pLDDT filtering to get the residues that would remain after pLDDT filtering
            # (if pLDDT filtering removes residues)
            remaining_after_plddt = []
            if plddt_cutoff_method in [
                "truncate_at_first_low_plddt",
                "remove_low_plddt_residues",
            ]:
                # Only consider residues that pass pLDDT filtering
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if residue_idx not in plddt_filtered_residues:
                        remaining_after_plddt.append(residue)
            else:
                # For exclude_low_plddt_from_stats, all residues remain (they're just
                # marked)
                remaining_after_plddt = dist_mat

            # Simulate the scan: build region incrementally from residues that passed pLDDT filtering
            # A residue is added if its PAE to ALL residues in the current region is <= max_pae
            # A residue is filtered if its PAE to ANY residue in the current region is > max_pae
            # Note: For exclude_low_plddt_from_stats, low pLDDT residues are still checked for PAE,
            # but when checking PAE from a new residue to the region, we skip low
            # pLDDT residues in the region
            for residue in remaining_after_plddt:
                residue_idx = residue.residue_index

                # Check if this residue should be added to the region based on PAE
                # Low pLDDT residues are still checked (they're not automatically added)
                should_add = True

                if hasattr(residue, "pae_array") and residue.pae_array is not None:
                    # Check PAE from this residue to each residue already in the region
                    # Skip low pLDDT residues in the region (they don't count for PAE
                    # filtering)
                    for region_residue in region:
                        region_residue_idx = region_residue.residue_index
                        # Skip low pLDDT residues when checking PAE (they don't block
                        # other residues)
                        if region_residue_idx in low_plddt_residues:
                            continue
                        try:
                            if region_residue_idx < len(residue.pae_array):
                                pae_to_region = residue.pae_array[region_residue_idx]
                                if pae_to_region > max_pae:
                                    # PAE exceeds threshold, so this residue should be
                                    # filtered
                                    should_add = False
                                    break
                        except (TypeError, AttributeError, IndexError):
                            # If we can't get PAE value, treat as should not add
                            # (conservative)
                            should_add = False
                            break
                # else: pae_array is None or doesn't exist
                # In Hail code: (residue.pae_array is None) | ... means if None, add it
                # So should_add = True is already set, which is correct

                if should_add:
                    region.append(residue)
                else:
                    pae_filtered_residues.add(residue_idx)
        else:
            # If no PAE array field, all residues that passed pLDDT filtering are in
            # the region
            if plddt_cutoff_method in [
                "truncate_at_first_low_plddt",
                "remove_low_plddt_residues",
            ]:
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if residue_idx not in plddt_filtered_residues:
                        region.append(residue)
            else:
                # For exclude_low_plddt_from_stats or no pLDDT filtering, all residues
                # are in the region
                region = dist_mat

        # Determine which residues are in the final region (not filtered)
        region_residue_indices = {res.residue_index for res in region}

        # Create column header (use variable for backslash since f-strings can't
        # have backslashes)
        from_to_label = "From\\To"
        header = f"        {from_to_label:<8}"
        for res_idx in all_residue_indices:
            # Color orange if residue is not in the final region (filtered out)
            is_filtered = res_idx not in region_residue_indices
            if is_filtered:
                header += f"{ORANGE}{res_idx:>8}{RESET}"
            else:
                header += f"{res_idx:>8}"
        # Add pLDDT column header if pLDDT filtering is used
        if min_plddt is not None and has_plddt_data:
            header += f"{'pLDDT':>8}"
        center_output += header + "\n"
        header_sep_len = 8 + 8 * len(all_residue_indices)
        if min_plddt is not None and has_plddt_data:
            header_sep_len += 8
        center_output += "        " + "-" * header_sep_len + "\n"

        # For each residue in the sorted array, show its PAE to all residues
        for i, residue_i in enumerate(dist_mat):
            residue_i_idx = residue_i.residue_index
            row_str = f"        {residue_i_idx:>8}"

            # Check if this source residue is filtered (not in final region)
            # If filtered, all PAE values in this row should be "-" since they're not
            # used for filtering
            is_source_filtered = residue_i_idx not in region_residue_indices

            # For each column (target residue)
            for j, residue_j in enumerate(dist_mat):
                residue_j_idx = residue_j.residue_index

                if j <= i:
                    # Lower triangle or diagonal: show "-" or "0.0" for self
                    if j == i:
                        # row_str += f"{'0.0':>8}"
                        row_str += f"{'-':>8}"
                    else:
                        row_str += f"{'-':>8}"
                else:
                    # Upper triangle: show actual PAE value
                    # If source residue is filtered, show "-" (its PAE values don't
                    # count for filtering)
                    if is_source_filtered:
                        pae_str = f"{'-':>8}"
                    # For exclude_low_plddt_from_stats, show "-" if source residue
                    # (residue_i) has low pLDDT
                    elif (
                        plddt_cutoff_method == "exclude_low_plddt_from_stats"
                        and residue_i_idx in low_plddt_residues
                    ):
                        pae_str = f"{'-':>8}"
                    else:
                        pae_val = None
                        if has_pae_data:
                            # First try pae_array field (preferred for pairwise PAE)
                            if (
                                has_pae_array_field
                                and hasattr(residue_i, "pae_array")
                                and residue_i.pae_array is not None
                            ):
                                try:
                                    if residue_j_idx < len(residue_i.pae_array):
                                        pae_val = residue_i.pae_array[residue_j_idx]
                                except (TypeError, AttributeError, IndexError):
                                    pass
                            # Fall back to pae field if it's an array
                            elif (
                                pae_is_array
                                and hasattr(residue_i, "pae")
                                and residue_i.pae is not None
                            ):
                                try:
                                    if residue_j_idx < len(residue_i.pae):
                                        pae_val = residue_i.pae[residue_j_idx]
                                except (TypeError, AttributeError, IndexError):
                                    pass

                        # Format PAE value (ensure consistent width accounting for ANSI
                        # codes)
                        if pae_val is not None:
                            if pae_val > max_pae:
                                # Format with color, ensuring the visible part is 5
                                # chars wide
                                pae_str = f"{RED}{pae_val:5.1f}{RESET}"
                            else:
                                pae_str = f"{pae_val:5.1f}"
                        else:
                            # No PAE data available
                            pae_str = "  ?  "

                    # Calculate visible width (without ANSI codes) and pad accordingly
                    # ANSI codes don't count toward visible width, so we need to pad the
                    # visible content
                    visible_pae_str = pae_str.replace(RED, "").replace(RESET, "")
                    padding_needed = max(0, 8 - len(visible_pae_str))
                    row_str += " " * padding_needed + pae_str

            # Add pLDDT column if pLDDT filtering is used
            if min_plddt is not None and has_plddt_data:
                plddt_val = None
                if hasattr(residue_i, "plddt"):
                    plddt_val = residue_i.plddt

                if plddt_val is not None:
                    if plddt_val < min_plddt:
                        plddt_str = f"{RED}{plddt_val:>8.1f}{RESET}"
                    else:
                        plddt_str = f"{plddt_val:>8.1f}"
                else:
                    plddt_str = "     NA"
                row_str += plddt_str

            center_output += row_str + "\n"

        # Add pLDDT row if pLDDT filtering is used
        if min_plddt is not None and has_plddt_data:
            plddt_row = f"        {'pLDDT':<8}"
            for residue in dist_mat:
                residue_idx = residue.residue_index
                plddt_val = None
                if hasattr(residue, "plddt"):
                    plddt_val = residue.plddt

                if plddt_val is not None:
                    if plddt_val < min_plddt:
                        plddt_str = f"{RED}{plddt_val:>8.1f}{RESET}"
                    else:
                        plddt_str = f"{plddt_val:>8.1f}"
                else:
                    plddt_str = "     NA"
                plddt_row += plddt_str
            # Add pLDDT column value (same as row label, or could show average/NA)
            plddt_row += f"{'pLDDT':>8}"
            center_output += "        " + "-" * header_sep_len + "\n"
            center_output += "\n" + plddt_row + "\n"

        # Store output for this center residue
        output_dict[center_idx] = center_output

    return output_dict


# ANSI color codes
RESET = "\033[0m"
RED = "\033[91m"  # High PAE (above cutoff), moderate constraint
BOLD = "\033[1m"
ORANGE = "\033[38;5;214m"  # Filtered residues (bright orange, 256-color mode)
# Bold, underline, green for minimum OE upper, very constrained
HIGHLIGHT = "\033[1m\033[4m\033[92m"
UNDERLINE = "\033[4m"
GREEN = "\033[92m"  # Good values (low AIC, very negative NLL)
YELLOW = "\033[93m"  # Moderate values


def get_debug_oe_table_string(
    oe_list, plddt_lookup=None, min_plddt=None, pae_lookup=None, max_pae=None
):
    # Print OE array (per-residue observed/expected)
    output_string = ""

    # Determine which columns to include
    has_plddt = (
        min_plddt is not None and plddt_lookup is not None and len(plddt_lookup) > 0
    )
    has_pae = max_pae is not None and pae_lookup is not None and len(pae_lookup) > 0

    # Add color legend
    legend_parts = []
    if has_pae:
        legend_parts.append(
            f"        {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}"
        )
    if has_plddt:
        legend_parts.append(
            f"        {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}"
        )
    if legend_parts:
        legend = "\n".join(legend_parts)
        output_string += f"\n    {BOLD}Color legend:{RESET}\n{legend}\n"

    # Get all unique residue indices for PAE columns
    # First, determine the residue field name from oe_list
    residue_field_name = None
    if len(oe_list) > 0:
        if "residue_index" in oe_list[0]:
            residue_field_name = "residue_index"
        elif "aa_index" in oe_list[0]:
            residue_field_name = "aa_index"

    # Get all residue indices from pae_lookup (all residues that have PAE data)
    all_residue_indices = []
    if has_pae:
        # Get residue indices from pae_lookup keys (all residues with PAE data)
        all_residue_indices = sorted(pae_lookup.keys())
        # Also include any residue indices from oe_list that might not be in pae_lookup
        if residue_field_name:
            oe_residue_indices = set([entry[residue_field_name] for entry in oe_list])
            all_residue_indices = sorted(set(all_residue_indices) | oe_residue_indices)

    # Build header
    header = f"    {'Residue Index':<15} {'Observed':<12} {'Expected':<12} {'O/E':<10}"
    if has_plddt:
        header += f" {'pLDDT':<10}"
    if has_pae:
        header += " "  # Empty column before PAE
        for res_idx in all_residue_indices:
            header += f"{res_idx:>8}"

    # Add PAE header row if PAE columns are present
    if has_pae:
        # Calculate the width of columns before PAE section
        pre_pae_width = 15 + 12 + 12 + 10  # Residue Index + Observed + Expected + O/E
        if has_plddt:
            pre_pae_width += 10  # pLDDT
        pre_pae_width += 4  # "    " prefix
        pre_pae_width += 1  # Empty column before PAE

        # Create PAE header row with "PAE to residue:" label
        pae_header = (
            "    " + " " * (15 + 1) + " " * (12 + 1) + " " * (12 + 1) + " " * (10 + 1)
        )
        if has_plddt:
            pae_header += " " * (10 + 1)
        pae_header += " "  # Empty column before PAE
        # Add "PAE to residue:" label, centered over the PAE columns
        pae_col_width = 8 * len(all_residue_indices)
        pae_label = "PAE to residue:"
        padding = max(0, (pae_col_width - len(pae_label)) // 2)
        pae_header += " " * padding + pae_label
        output_string += pae_header + "\n"

    output_string += header + "\n"

    # Calculate separator length
    sep_len = 50
    if has_plddt:
        sep_len += 11
    if has_pae:
        sep_len += 1 + 8 * len(all_residue_indices)
    output_string += "    " + "-" * sep_len + "\n"

    residue_field_name = None
    if len(oe_list) > 0:
        if "residue_index" in oe_list[0]:
            residue_field_name = "residue_index"
        elif "aa_index" in oe_list[0]:
            residue_field_name = "aa_index"
        else:
            oe_list = [x.annotate(residue_index=i) for i, x in enumerate(oe_list)]
            residue_field_name = "residue_index"

    for oe_entry in oe_list:
        residue_idx = oe_entry[residue_field_name]
        obs = oe_entry.obs
        exp = oe_entry.exp
        oe_ratio = obs / exp if exp > 0 else 0.0

        row_str = f"    {residue_idx:<15} {obs:<12} {exp:<12.2f} {oe_ratio:<10.3f}"

        if has_plddt:
            plddt_val = plddt_lookup.get(residue_idx, None)
            if plddt_val is not None:
                # Color low pLDDT (below min_plddt) red
                if plddt_val < min_plddt:
                    plddt_str = f"{RED}{plddt_val:<10.1f}{RESET}"
                else:
                    plddt_str = f"{plddt_val:<10.1f}"
            else:
                plddt_str = "NA        "
            row_str += f" {plddt_str}"

        if has_pae:
            row_str += " "  # Empty column before PAE
            # Get PAE array for this residue
            pae_array = pae_lookup.get(residue_idx, None)
            for target_res_idx in all_residue_indices:
                pae_val = None
                if pae_array is not None:
                    try:
                        if target_res_idx < len(pae_array):
                            pae_val = pae_array[target_res_idx]
                    except (TypeError, IndexError):
                        pass

                if pae_val is not None:
                    if pae_val > max_pae:
                        pae_str = f"{RED}{pae_val:>8.2f}{RESET}"
                    else:
                        pae_str = f"{pae_val:>8.2f}"
                else:
                    pae_str = "     NA"
                row_str += pae_str

        output_string += row_str + "\n"

    return output_string


def debug_print_oe_table(
    oe_expr: hl.expr.ArrayExpression,
    min_plddt: Optional[float] = None,
    max_pae: Optional[float] = None,
    dist_mat_expr: Optional[hl.expr.ArrayExpression] = None,
) -> None:
    output_string = "\n\n\n"

    # Get first key for debugging.
    uniprot_id = None
    transcript_id = None
    ht = oe_expr._indices.source
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter((ht.uniprot_id == uniprot_id) & (ht.enst == transcript_id))

    output_string += (
        f"{BOLD}=== Per-residue observed and expected values ==={RESET}\n\n"
    )

    # Collect data - also get plddt_lookup, pae_lookup, and dist_mat if available
    select_fields = ["oe"]
    # Check if fields exist in the row schema
    row_schema = _ht_debug.row.dtype
    has_plddt_lookup = "plddt_lookup" in row_schema.fields
    has_pae_lookup = "pae_lookup" in row_schema.fields
    has_dist_mat = "dist_mat" in row_schema.fields
    if has_plddt_lookup:
        select_fields.append("plddt_lookup")
    if has_pae_lookup:
        select_fields.append("pae_lookup")
    if has_dist_mat and max_pae is not None:
        # Also get dist_mat to extract PAE arrays if pae_lookup is not available
        select_fields.append("dist_mat")

    debug_data = _ht_debug.select(*select_fields).collect()

    if not debug_data:
        return

    row = debug_data[0]

    # Build pLDDT lookup if available
    plddt_lookup = None
    if (
        has_plddt_lookup
        and hasattr(row, "plddt_lookup")
        and row.plddt_lookup is not None
    ):
        # Create dictionary mapping residue_index to pLDDT value
        plddt_lookup = {
            entry.residue_index: entry.plddt
            for entry in row.plddt_lookup
            if entry is not None
        }

    # Build PAE lookup if available
    # PAE lookup structure: array of structs with {residue_index: int,
    # pae_array: array<float>}
    pae_lookup = None
    if max_pae is not None:
        # First try to get it from dist_mat_expr parameter (most reliable since
        # it's passed directly)
        if dist_mat_expr is not None:
            # Collect dist_mat_expr to extract PAE arrays
            _ht_dist = dist_mat_expr._indices.source
            _ht_dist = _ht_dist.annotate(dist_mat=dist_mat_expr)
            _ht_dist = _ht_dist.filter(
                (_ht_dist.uniprot_id == uniprot_id) & (_ht_dist.enst == transcript_id)
            )
            dist_data = _ht_dist.select("dist_mat").head(1).collect()
            if dist_data and len(dist_data) > 0:
                dist_row = dist_data[0]
                if (
                    hasattr(dist_row, "dist_mat")
                    and dist_row.dist_mat is not None
                    and len(dist_row.dist_mat) > 0
                ):
                    pae_lookup = {}
                    for entry in dist_row.dist_mat:
                        if (
                            entry is not None
                            and hasattr(entry, "residue_index")
                            and hasattr(entry, "pae_array")
                            and entry.pae_array is not None
                        ):
                            pae_lookup[entry.residue_index] = entry.pae_array
        # If not available from dist_mat_expr, try pae_lookup field in table
        if (
            (pae_lookup is None or len(pae_lookup) == 0)
            and has_pae_lookup
            and hasattr(row, "pae_lookup")
            and row.pae_lookup is not None
        ):
            # Create dictionary mapping residue_index to PAE array
            pae_lookup = {}
            for entry in row.pae_lookup:
                if (
                    entry is not None
                    and hasattr(entry, "residue_index")
                    and hasattr(entry, "pae_array")
                ):
                    pae_lookup[entry.residue_index] = entry.pae_array
        # If still not available, try to get it from dist_mat if available
        if (
            (pae_lookup is None or len(pae_lookup) == 0)
            and has_dist_mat
            and hasattr(row, "dist_mat")
            and row.dist_mat is not None
        ):
            # Extract PAE arrays from dist_mat entries
            pae_lookup = {}
            for entry in row.dist_mat:
                if (
                    entry is not None
                    and hasattr(entry, "residue_index")
                    and hasattr(entry, "pae_array")
                    and entry.pae_array is not None
                ):
                    pae_lookup[entry.residue_index] = entry.pae_array

    # Print OE array (per-residue observed/expected)
    output_string += get_debug_oe_table_string(
        row.oe,
        plddt_lookup=plddt_lookup,
        min_plddt=min_plddt,
        pae_lookup=pae_lookup,
        max_pae=max_pae,
    )

    logger.info(output_string)


def debug_print_oe_and_regions(
    ht: hl.Table,
    title: str = "OE and Regions",
) -> None:
    """
    Print debug output showing OE values and min_oe_upper regions for a given UniProt/transcript.

    :param ht: Hail Table with schema containing 'oe' and 'min_oe_upper' arrays
    :param uniprot_id: UniProt ID to filter to
    :param transcript_id: Transcript ID to filter to
    :param title: Title for the debug output
    :param pae_cutoff: PAE threshold - values above this are colored red, below are uncolored. Default is 15.0
    """
    output_string = "\n\n\n"

    # Get first key for debugging.
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.transcript_id.collect()[0]

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter(
        (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    output_string += f"{BOLD}=== {title} ==={RESET}\nUniProt ID: {uniprot_id}, Transcript ID: {transcript_id}\n"

    # Collect data
    debug_data = _ht_debug.select("oe", "min_oe_upper").collect()

    if not debug_data:
        return

    row = debug_data[0]

    # Print min_oe_upper array (best region for each center residue)
    output_string += f"\n    {BOLD}=== Min OE Upper Array (best region for each center residue) ==={RESET}\n\n"

    if not row.min_oe_upper or len(row.min_oe_upper) == 0:
        output_string += "No min_oe_upper entries found.\n"
        logger.info(output_string)
        return

    # Find the entry with the minimum OE value (not oe_upper)
    min_oe_idx = None
    min_oe_val = None
    for idx, entry in enumerate(row.min_oe_upper):
        if hasattr(entry, "oe") and entry.oe is not None:
            if min_oe_val is None or entry.oe < min_oe_val:
                min_oe_val = entry.oe
                min_oe_idx = idx

    for idx, entry in enumerate(row.min_oe_upper):
        is_min = min_oe_idx is not None and idx == min_oe_idx

        output_string += f"\n    {BOLD}Center Residue: {entry.residue_index}{RESET}\n"
        output_string += f"      Observed: {entry.obs:7d}\n"
        output_string += f"      Expected: {entry.exp:7.2f}\n"

        # Highlight minimum OE
        if is_min:
            output_string += (
                f"      {HIGHLIGHT}O/E:      {entry.oe:7.2f} (MINIMUM){RESET}\n"
            )
        else:
            output_string += f"      O/E:      {entry.oe:7.2f}\n"

        output_string += f"      OE Upper: {entry.oe_upper:7.2f}\n"

        # Region (array of residue indices)
        if hasattr(entry, "region") and entry.region:
            region_str = ", ".join([f"{r:7d}" for r in entry.region])
            output_string += (
                f"      Region residues ({len(entry.region)}): {region_str}\n"
            )
        else:
            output_string += f"      Region: None or empty\n"

        # Distances (array of distances)
        if hasattr(entry, "dists") and entry.dists:
            dists_str = ", ".join([f"{d:7.2f}" for d in entry.dists])
            output_string += (
                f"      Distances ({len(entry.dists)}):       {dists_str}\n"
            )
        else:
            output_string += f"      Distances: None or empty\n"

    logger.info(output_string)


def debug_print_dist_mat_with_colors(
    dist_mat_expr: hl.expr.ArrayExpression,
    title: str = "Distance matrix with PAE",
    max_pae: Optional[float] = None,
    min_plddt: Optional[float] = None,
    oe_expr: Optional[hl.expr.ArrayExpression] = None,
    min_res_idx: Optional[int] = None,
    max_res_idx: Optional[int] = None,
    min_dist: Optional[float] = None,
    max_dist: Optional[float] = None,
) -> Dict[int, str]:
    """
    Generate colored debug output showing distance matrix with PAE values color-coded.
    Returns a dictionary where keys are center residue indices and values are output strings.

    :param dist_mat_expr: Array expression with distance matrix entries (must have dist, residue_index, and optionally pae, plddt)
    :param title: Title for the debug output
    :param max_pae: Maximum allowed PAE value - values above this are colored red, below are uncolored. If None, PAE won't be shown.
    :param min_plddt: Minimum allowed pLDDT value. If None, pLDDT won't be shown.
    :param oe_expr: Optional array expression with OE data (obs, exp, oe, oe_upper). If provided, will print
        observed, expected, cumulative observed, cumulative expected, cumulative O/E, and cumulative OE upper.
    :return: Dictionary mapping center residue indices to their debug output strings.
        For the "Before sorting" full matrix case, returns a dict with key -1 and the full matrix output.
    """

    # Get first key for debugging.
    ht = dist_mat_expr._indices.source
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    output_string = ""

    # For "Before sorting by nearest neighbor", if PAE and pLDDT are not shown,
    # show the full distance matrix (all residues as rows, distances to all
    # residues as columns)
    show_full_matrix = (
        "Before sorting by nearest neighbor" in title
        and max_pae is None
        and min_plddt is None
        and oe_expr is None
    )

    # Filter and collect data
    select_fields = {"dist_mat": dist_mat_expr}
    if oe_expr is not None:
        select_fields["oe_expr"] = oe_expr

    _ht_debug = ht.annotate(**select_fields)
    _ht_debug = _ht_debug.filter(
        (_ht_debug.uniprot_id == uniprot_id) & (_ht_debug.enst == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    header_string = ""
    if not show_full_matrix:
        header_string += "    "

    header_string += f"{BOLD}=== {title} ==={RESET}\n"

    # Collect data for colored output
    select_fields = ["dist_mat"]
    if oe_expr is not None:
        select_fields.append("oe_expr")

    debug_data = _ht_debug.select(*select_fields).collect()

    if not debug_data:
        return {}

    # Dictionary to store output per center residue
    output_dict = {}

    # Check if PAE and pLDDT fields exist in the struct type (only if requested)
    has_pae_field = False
    has_plddt_field = False
    if debug_data and len(debug_data[0].dist_mat) > 0:
        first_entry = debug_data[0].dist_mat[0]
        if max_pae is not None:
            has_pae_field = hasattr(first_entry, "pae")
        if min_plddt is not None:
            has_plddt_field = hasattr(first_entry, "plddt")

    # Check if OE fields exist
    has_oe_data = oe_expr is not None and debug_data and len(debug_data[0].oe_expr) > 0

    # ANSI color codes for distances (blue to white gradient)
    BLUE = "\033[94m"  # Closest (bright blue)
    CYAN = "\033[96m"  # Medium-close (cyan)
    WHITE = "\033[97m"  # Furthest (white)

    # ANSI color codes for residue indices (gradient from green to magenta)
    GREEN_IDX = "\033[92m"  # First in sorted array (closest)
    YELLOW_IDX = "\033[93m"  # Middle
    MAGENTA_IDX = "\033[95m"  # Last in sorted array (furthest)

    # ANSI color code for highlighting minimum OE upper
    HIGHLIGHT = "\033[1m\033[4m\033[92m"  # Bold and underline for minimum OE upper

    # Add color legends
    legend_parts = [
        f"            {BOLD}Residue index colors:{RESET} {GREEN_IDX}Lowest index{RESET} → {YELLOW_IDX}Middle{RESET} → {MAGENTA_IDX}Highest index{RESET}",
        f"            {BOLD}Distance colors:{RESET} {BLUE}Closest{RESET} → {CYAN}Medium{RESET} → {WHITE}Furthest{RESET}",
    ]
    if max_pae is not None and has_pae_field:
        legend_parts.append(
            f"            {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}"
        )
    if min_plddt is not None and has_plddt_field:
        legend_parts.append(
            f"            {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}"
        )
    legend = "\n".join(legend_parts)
    legend_string = f"\n        {BOLD}Color legend:{RESET}\n{legend}\n\n"

    # If showing full matrix, we need all rows to build it, but only show it once
    if show_full_matrix:
        # Build full distance matrix from all center residues
        # First, get all unique residue indices from the first row
        if not debug_data or len(debug_data) == 0:
            return {}

        first_row = debug_data[0]
        all_unique_residues = sorted(
            [entry.residue_index for entry in first_row.dist_mat]
        )

        # Build distance lookup: center_residue -> {target_residue: distance}
        dist_lookup = {}
        for center_row in debug_data:
            center_res = center_row.aa_index
            if center_res not in dist_lookup:
                dist_lookup[center_res] = {}
            for entry in center_row.dist_mat:
                target_res = entry.residue_index
                dist_lookup[center_res][target_res] = entry.dist

        # Calculate min/max distances for coloring
        all_distances = []
        for center_res in dist_lookup:
            for target_res in dist_lookup[center_res]:
                all_distances.append(dist_lookup[center_res][target_res])

        if all_distances:
            row_min_dist = min(all_distances) if min_dist is None else min_dist
            row_max_dist = max(all_distances) if max_dist is None else max_dist
            dist_range = (
                row_max_dist - row_min_dist if row_max_dist > row_min_dist else 1.0
            )
        else:
            row_min_dist = 0.0
            row_max_dist = 1.0
            dist_range = 1.0

        # Calculate residue index range for coloring
        if all_unique_residues:
            row_min_res_idx = (
                min(all_unique_residues) if min_res_idx is None else min_res_idx
            )
            row_max_res_idx = (
                max(all_unique_residues) if max_res_idx is None else max_res_idx
            )
            res_idx_range = (
                row_max_res_idx - row_min_res_idx
                if row_max_res_idx > row_min_res_idx
                else 1.0
            )
        else:
            row_min_res_idx = 0
            row_max_res_idx = 1
            res_idx_range = 1.0

        # Helper function to get color for residue index
        def get_residue_color(res_idx):
            if res_idx_range > 0:
                normalized_res_idx = (res_idx - row_min_res_idx) / res_idx_range
                if normalized_res_idx < 0.33:
                    return GREEN_IDX
                elif normalized_res_idx < 0.67:
                    return YELLOW_IDX
                else:
                    return MAGENTA_IDX
            else:
                return GREEN_IDX

        # Create matrix header with colored residue indices
        header = f"\n        {'Residue':<10}"
        for res_idx in all_unique_residues:
            color = get_residue_color(res_idx)
            header += f"{color}{res_idx:>8}{RESET}"
        output_string += header + "\n"
        output_string += "        " + "-" * (10 + 8 * len(all_unique_residues)) + "\n"

        # Create matrix rows with colored row headers
        for from_res in all_unique_residues:
            if from_res in dist_lookup:
                row_color = get_residue_color(from_res)
                row_str = f"        {row_color}{from_res:<10}{RESET}"
                for to_res in all_unique_residues:
                    if to_res in dist_lookup[from_res]:
                        dist_val = dist_lookup[from_res][to_res]
                        # Color code distance
                        if dist_range > 0:
                            normalized = (dist_val - row_min_dist) / dist_range
                            if normalized < 0.5:
                                color = BLUE if normalized < 0.25 else CYAN
                            else:
                                color = CYAN if normalized < 0.75 else WHITE
                        else:
                            color = BLUE
                        row_str += f"{color}{dist_val:>8.2f}{RESET}"
                    else:
                        row_str += f"{'?':>8}"
                output_string += row_str + "\n"

        # For full matrix, store with key -1 (special key for full matrix)
        full_output = "\n\n" + header_string + legend_string + output_string
        output_dict[-1] = full_output
        return output_dict

    # Normal processing: show data for each center residue
    for row in debug_data:
        center_idx = row.aa_index

        # Use provided min/max values if available, otherwise calculate from current row
        residue_indices = [entry.residue_index for entry in row.dist_mat]
        if not residue_indices:
            # Empty dist_mat - skip this row
            continue

        if min_res_idx is None or max_res_idx is None:
            # Calculate from current row if not provided
            row_min_res_idx = min(residue_indices)
            row_max_res_idx = max(residue_indices)
        else:
            # Use provided values
            row_min_res_idx = min_res_idx
            row_max_res_idx = max_res_idx
        res_idx_range = (
            row_max_res_idx - row_min_res_idx
            if row_max_res_idx > row_min_res_idx
            else 1
        )  # Avoid division by zero

        # Find min and max distances for gradient
        distances = [entry.dist for entry in row.dist_mat]
        if not distances:
            # Empty dist_mat - skip this row
            continue

        if min_dist is None or max_dist is None:
            # Calculate from current row if not provided
            row_min_dist = min(distances)
            row_max_dist = max(distances)
        else:
            # Use provided values
            row_min_dist = min_dist
            row_max_dist = max_dist
        dist_range = (
            row_max_dist - row_min_dist if row_max_dist > row_min_dist else 1.0
        )  # Avoid division by zero

        # Get OE data if available (should be parallel to dist_mat)
        oe_entries = row.oe_expr if has_oe_data else None

        # Find minimum OE upper value and its index (if OE data is available)
        min_oe_upper_idx = None
        min_oe_upper_val = None
        if has_oe_data and oe_entries:
            # Filter to entries that have oe_upper defined and not None
            valid_oe_entries = [
                (idx, entry)
                for idx, entry in enumerate(oe_entries)
                if hasattr(entry, "oe_upper") and entry.oe_upper is not None
            ]
            if valid_oe_entries:
                # Find the entry with minimum OE upper
                min_entry = min(valid_oe_entries, key=lambda x: x[1].oe_upper)
                min_oe_upper_idx = min_entry[0]
                min_oe_upper_val = min_entry[1].oe_upper

        # Build colored output
        residue_strs = []
        dist_strs = []
        pae_strs = []
        plddt_strs = []
        obs_strs = []
        exp_strs = []
        cum_obs_strs = []
        cum_exp_strs = []
        cum_oe_strs = []
        cum_oe_upper_strs = []

        for idx, entry in enumerate(row.dist_mat):
            # Color code residue indices by their actual index value (not position in
            # array)
            res_idx_val = entry.residue_index
            if res_idx_range > 0:
                # Normalize residue index to 0-1 range
                normalized_res_idx = (res_idx_val - row_min_res_idx) / res_idx_range
                # Interpolate between green (lowest index) and magenta (highest index)
                if normalized_res_idx < 0.33:
                    color_idx = GREEN_IDX
                elif normalized_res_idx < 0.67:
                    color_idx = YELLOW_IDX
                else:
                    color_idx = MAGENTA_IDX
            else:
                # All residue indices are the same
                color_idx = GREEN_IDX
            residue_strs.append(f"{color_idx}{res_idx_val:7d}{RESET}")

            # Color code distances (blue for closest, white for furthest)
            dist_val = entry.dist
            if dist_range > 0:
                # Normalize distance to 0-1 range
                normalized = (dist_val - row_min_dist) / dist_range
                # Interpolate between blue (0.0) and white (1.0)
                if normalized < 0.5:
                    # Blue to cyan gradient
                    color = BLUE if normalized < 0.25 else CYAN
                else:
                    # Cyan to white gradient
                    color = CYAN if normalized < 0.75 else WHITE
            else:
                # All distances are the same
                color = BLUE
            dist_strs.append(f"{color}{dist_val:7.2f}{RESET}")

            # Color code PAE values if available and requested
            if max_pae is not None and has_pae_field:
                pae_val = entry.pae
                if pae_val > max_pae:
                    # Color red if above cutoff
                    pae_strs.append(f"{RED}{pae_val:7.2f}{RESET}")
                else:
                    # No color if below cutoff
                    pae_strs.append(f"{pae_val:7.2f}")
            elif max_pae is not None:
                pae_strs.append("     NA")

            # Add pLDDT values if available and requested
            if min_plddt is not None and has_plddt_field:
                plddt_val = entry.plddt
                if plddt_val is not None:
                    # Color low pLDDT (below min_plddt) red
                    if plddt_val < min_plddt:
                        plddt_strs.append(f"{RED}{plddt_val:7.1f}{RESET}")
                    else:
                        plddt_strs.append(f"{plddt_val:7.1f}")
                else:
                    plddt_strs.append("     NA")
            elif min_plddt is not None:
                plddt_strs.append("     NA")

            # Add OE data if available (no coloring)
            # Check if this residue is excluded from stats (for
            # exclude_low_plddt_from_stats)
            is_excluded = False
            if hasattr(entry, "exclude_from_stats"):
                is_excluded = entry.exclude_from_stats is True

            if has_oe_data and oe_entries and idx < len(oe_entries):
                oe_entry = oe_entries[idx]

                # If excluded from stats, show NA for all OE values
                if is_excluded:
                    obs_strs.append("     NA")
                    exp_strs.append("     NA")
                    cum_obs_strs.append("     NA")
                    cum_exp_strs.append("     NA")
                    cum_oe_strs.append("     NA")
                    cum_oe_upper_strs.append("     NA")
                else:
                    # After get_cumulative_oe, obs and exp are cumulative
                    # Calculate per-residue values from cumulative
                    if idx == 0:
                        # First entry: per-residue = cumulative
                        per_residue_obs = oe_entry.obs
                        per_residue_exp = oe_entry.exp
                    else:
                        # Subsequent entries: per-residue = current cumulative -
                        # previous cumulative
                        prev_entry = oe_entries[idx - 1]
                        per_residue_obs = (
                            oe_entry.obs - prev_entry.obs
                            if (oe_entry.obs is not None and prev_entry.obs is not None)
                            else None
                        )
                        per_residue_exp = (
                            oe_entry.exp - prev_entry.exp
                            if (oe_entry.exp is not None and prev_entry.exp is not None)
                            else None
                        )

                    # Handle None values in formatting
                    if per_residue_obs is not None:
                        obs_strs.append(f"{per_residue_obs:7d}")
                    else:
                        obs_strs.append("     NA")

                    if per_residue_exp is not None:
                        exp_strs.append(f"{per_residue_exp:7.2f}")
                    else:
                        exp_strs.append("     NA")

                    # Cumulative values
                    if oe_entry.obs is not None:
                        cum_obs_strs.append(f"{oe_entry.obs:7d}")
                    else:
                        cum_obs_strs.append("     NA")

                    if oe_entry.exp is not None:
                        cum_exp_strs.append(f"{oe_entry.exp:7.2f}")
                    else:
                        cum_exp_strs.append("     NA")

                    if oe_entry.oe is not None:
                        cum_oe_strs.append(f"{oe_entry.oe:7.2f}")
                    else:
                        cum_oe_strs.append("     NA")

                    # Highlight minimum OE upper value
                    if oe_entry.oe_upper is not None:
                        if min_oe_upper_idx is not None and idx == min_oe_upper_idx:
                            cum_oe_upper_strs.append(
                                f"{HIGHLIGHT}{oe_entry.oe_upper:7.2f}{RESET}"
                            )
                        else:
                            cum_oe_upper_strs.append(f"{oe_entry.oe_upper:7.2f}")
                    else:
                        cum_oe_upper_strs.append("     NA")
            elif has_oe_data:
                # OE data exists but this index is out of range
                obs_strs.append("     NA")
                exp_strs.append("     NA")
                cum_obs_strs.append("     NA")
                cum_exp_strs.append("     NA")
                cum_oe_strs.append("     NA")
                cum_oe_upper_strs.append("     NA")

        # Build output for this center residue (content only, no header/legend)
        output = (
            f"        Residue indices:     {', '.join(residue_strs)}\n"
            f"        Distances:           {', '.join(dist_strs)}\n"
        )

        if max_pae is not None and has_pae_field:
            output += f"        PAE:                 {', '.join(pae_strs)}\n"

        if min_plddt is not None and has_plddt_field:
            output += f"        pLDDT:               {', '.join(plddt_strs)}\n"

        if has_oe_data:
            output += (
                f"        Observed:            {', '.join(obs_strs)}\n"
                f"        Expected:            {', '.join(exp_strs)}\n"
                f"        Cumulative Observed: {', '.join(cum_obs_strs)}\n"
                f"        Cumulative Expected: {', '.join(cum_exp_strs)}\n"
                f"        Cumulative O/E:      {', '.join(cum_oe_strs)}\n"
                f"        Cumulative OE Upper: {', '.join(cum_oe_upper_strs)}\n"
            )
            if min_oe_upper_idx is not None:
                min_oe_entry = oe_entries[min_oe_upper_idx]
                output += f"\n        {BOLD}Minimum OE Upper:{RESET} {HIGHLIGHT}{min_oe_upper_val:.2f}{RESET} at index {min_oe_upper_idx} (residue {row.dist_mat[min_oe_upper_idx].residue_index})\n\n"

                # Add Observed, Expected, O/E, and OE Upper values
                if hasattr(min_oe_entry, "obs") and min_oe_entry.obs is not None:
                    output += f"            Region observed:    {min_oe_entry.obs:7d}\n"
                else:
                    output += f"            Region observed:    {'NA':>7}\n"

                if hasattr(min_oe_entry, "exp") and min_oe_entry.exp is not None:
                    output += (
                        f"            Region expected:    {min_oe_entry.exp:7.2f}\n"
                    )
                else:
                    output += f"            Region expected:    {'NA':>7}\n"

                if hasattr(min_oe_entry, "oe") and min_oe_entry.oe is not None:
                    output += (
                        f"            Region O/E:         {min_oe_entry.oe:7.2f}\n"
                    )
                else:
                    output += f"            Region O/E:         {'NA':>7}\n"

                if (
                    hasattr(min_oe_entry, "oe_upper")
                    and min_oe_entry.oe_upper is not None
                ):
                    output += f"            Region OE Upper:    {min_oe_entry.oe_upper:7.2f}\n"
                else:
                    output += f"            Region OE Upper:    {'NA':>7}\n"

                # Add region information if available
                if (
                    hasattr(min_oe_entry, "region")
                    and min_oe_entry.region is not None
                    and len(min_oe_entry.region) > 0
                ):
                    region_str = ", ".join([f"{r}" for r in min_oe_entry.region])
                    output += f"            Region residues ({len(min_oe_entry.region)}): {region_str}\n"
                else:
                    # If region is not available, build it from dist_mat up to
                    # min_oe_upper_idx
                    region_residues = [
                        entry.residue_index
                        for entry in row.dist_mat[: min_oe_upper_idx + 1]
                    ]
                    region_str = ", ".join([f"{r}" for r in region_residues])
                    output += f"            Region residues ({len(region_residues)}): {region_str}\n"

        # Store output for this center residue
        # Include header and legend in the output string
        full_output = "\n\n" + header_string + legend_string + output
        output_dict[center_idx] = full_output

    return output_dict


def get_3d_residue(
    center_residue_index_expr: hl.expr.Int32Expression,
    dist_mat_expr: hl.expr.ArrayExpression,
    oe_expr: hl.expr.ArrayExpression,
    alpha: float = 0.05,
    oe_upper_method: str = "gamma",
    min_exp_mis: int = MIN_EXP_MIS,
    pae_cutoff_method: str = "truncate_on_pairwise_pae_with_center",
    plddt_cutoff_method: Optional[str] = None,
    debug: bool = False,
    max_pae: Optional[float] = None,
    min_plddt: Optional[float] = None,
) -> hl.expr.StructExpression:
    """
    Get the 3D residue with the lowest upper bound of the OE confidence interval.

    :param dist_mat_expr: Array[Struct] with at least `residue_index` and `dist`.
    :param oe_expr: Array expression with observed and expected values.
    :param pae_expr: Optional array expression with PAE (Predicted Aligned Error) values.
        - For center-based options: 1D array of PAE(center_residue, neighbor_residue).
        - For region-based option: 2D array of PAE(region_i, region_j) where
          pae_expr[i][j] = PAE(residue_i, residue_j) for this region.
        Default is None.
    :param plddt_expr: Optional array expression with pLDDT (per-residue confidence) values.
        1D array of pLDDT scores for each residue, parallel to dist_mat_expr.
        Default is None.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :param oe_upper_method: Method to use for calculating OE upper bound. Default is "gamma".
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is MIN_EXP_MIS.
    :param max_pae: Maximum allowed PAE for including residues in the region. Default is 15.
    :param pae_cutoff_method: Strategy for handling residues with high PAE values. Options:
        - "truncate_on_pairwise_pae_with_center": Remove all residues after the first one
          with PAE(center, neighbor) > max_pae (sequential cutoff).
        - "filter_on_pairwise_pae_with_center": Remove only residues with
          PAE(center, neighbor) > max_pae, keep others (filter individual residues).
        - "filter_on_pairwise_pae_in_region": Filter out residues whose maximum pairwise
          PAE to any residue in the region exceeds max_pae. Requires 2D pae_expr.
        - "remove_all_residues_after_first_over_cutoff": Alias for
          "truncate_on_pairwise_pae_with_center" (for backward compatibility).
        - "remove_residues_after_first_high_pairwise_pae_with_current_residue": Alias for
          "truncate_on_pairwise_pae_with_center" (for backward compatibility).
        - "remove_residues_with_high_pairwise_pae_with_current_residue": Alias for
          "filter_on_pairwise_pae_with_center" (for backward compatibility).
        - "remove_residues_with_high_pairwise_pae_in_region": Alias for
          "filter_on_pairwise_pae_in_region" (for backward compatibility).
        Default is "truncate_on_pairwise_pae_with_center".
    :param min_plddt: Minimum allowed pLDDT for including residues in the region.
        Default is None.
    :param plddt_cutoff_method: Strategy for handling residues with low pLDDT values. Options:
        - "truncate_at_first_low_plddt": Remove all residues after the first one with
          pLDDT < min_plddt (sequential cutoff).
        - "remove_low_plddt_residues": Remove only residues with pLDDT < min_plddt,
          keep others (filter individual residues).
        - "mask_low_confidence_plddt": Mask (set to missing) residues with pLDDT < min_plddt
          instead of removing them.
        - "exclude_low_plddt_from_stats": Keep low pLDDT residues in regions (for assignment
          based on distance) but exclude them from statistical calculations (Obs, Exp, OE,
          OE upper, nLL, etc.).
        Default is None (no pLDDT filtering).
    :return: Struct expression with the 3D residue with the lowest upper bound of the OE
        confidence interval.
    """
    # Dictionary to collect all debug outputs, keyed by center residue index
    # Structure: {center_res_idx: {step_name: output_string}}
    debug_outputs_by_residue = {}

    # Extract UniProt ID and Transcript ID once for debug output
    debug_uniprot_id = None
    debug_transcript_id = None
    if debug:
        _ht_debug = dist_mat_expr._indices.source
        first_row = _ht_debug.head(1)
        if first_row.count() > 0:
            debug_uniprot_id = first_row.uniprot_id.collect()[0]
            debug_transcript_id = first_row.enst.collect()[0]
        # Print UniProt ID and Transcript ID once at the beginning
        logger.info(
            f"\nUniProt ID: {debug_uniprot_id}, Transcript ID: {debug_transcript_id}\n"
        )
        debug_print_oe_table(
            oe_expr, min_plddt=min_plddt, max_pae=max_pae, dist_mat_expr=dist_mat_expr
        )

    # Annotate neighbor order per residue.
    dist_mat_expr = add_idx_to_array(
        dist_mat_expr, "residue_index", element_name="dist"
    )

    # Calculate global min/max for consistent coloring across debug statements
    # We need to collect these values once and reuse them
    debug_min_max = None
    if debug:
        # Get min/max values from the original dist_mat_expr for consistent coloring
        # We'll calculate these from the first row's dist_mat
        _ht_debug = dist_mat_expr._indices.source
        _ht_debug = _ht_debug.annotate(dist_mat=dist_mat_expr).head(1)
        if _ht_debug.count() > 0:
            debug_data = _ht_debug.select("dist_mat").collect()
            if debug_data and len(debug_data[0].dist_mat) > 0:
                all_residue_indices = [
                    entry.residue_index for entry in debug_data[0].dist_mat
                ]
                all_distances = [entry.dist for entry in debug_data[0].dist_mat]
                if all_residue_indices and all_distances:
                    debug_min_max = {
                        "min_res_idx": min(all_residue_indices),
                        "max_res_idx": max(all_residue_indices),
                        "min_dist": min(all_distances),
                        "max_dist": max(all_distances),
                    }

    # Debug: Show dist_mat_expr before sorting by nearest neighbor.
    if debug:
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: Before sorting by nearest neighbor",
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output (key -1 is for full matrix, otherwise it's center_res_idx)
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx]["Before sorting"] = output_str

    dist_mat_expr = hl.sorted(dist_mat_expr, key=lambda x: x.dist)

    # Debug: Show sorting by nearest neighbor.
    if debug:
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: After sorting by nearest neighbor",
            max_pae=max_pae,
            min_plddt=min_plddt,
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx]["After sorting"] = output_str

    drop_fields = []
    # Apply pLDDT filtering FIRST (before PAE filtering)
    # This makes sense because pLDDT is a per-residue confidence score - if a residue
    # has low confidence, we should filter it out before doing pairwise PAE checks
    if min_plddt is not None:
        # Check if pLDDT data is available in dist_mat
        has_plddt_field = "plddt" in dist_mat_expr.dtype.element_type.fields

        if not has_plddt_field:
            raise ValueError("pLDDT data is not available in dist_mat_expr!")

        if plddt_cutoff_method == "truncate_at_first_low_plddt":
            # Remove all residues after the first one with pLDDT < min_plddt
            plddt_cutoff_idx = hl.enumerate(dist_mat_expr).find(
                lambda x: hl.is_defined(x[1].plddt) & (x[1].plddt < min_plddt)
            )
            dist_mat_expr = hl.if_else(
                hl.is_defined(plddt_cutoff_idx),
                dist_mat_expr[: plddt_cutoff_idx[0]],
                dist_mat_expr,
            )

        elif plddt_cutoff_method == "remove_low_plddt_residues":
            # Remove only residues with pLDDT < min_plddt, keep others
            dist_mat_expr = dist_mat_expr.filter(
                lambda x: ~hl.is_defined(x.plddt) | (x.plddt >= min_plddt)
            )

        elif plddt_cutoff_method == "exclude_low_plddt_from_stats":
            # Keep all residues but mark low confidence ones to exclude from stats
            # Mark residues with exclude_from_stats flag but keep them in the array
            dist_mat_expr = dist_mat_expr.map(
                lambda x: x.annotate(
                    exclude_from_stats=hl.is_defined(x.plddt) & (x.plddt < min_plddt)
                )
            )

        else:
            raise ValueError(f"Unknown plddt_cutoff_method: {plddt_cutoff_method}")

    # Apply PAE filtering AFTER pLDDT filtering
    if max_pae is not None:
        # Check if PAE data is available in dist_mat (as arrays)
        # Each residue may have a 'pae_array' field (full PAE array) or 'pae'
        # field (single value)
        has_pae_arrays = (
            "pae" in dist_mat_expr.dtype.element_type.fields
            or "pae_array" in dist_mat_expr.dtype.element_type.fields
        )

        if not has_pae_arrays:
            raise ValueError("PAE data is not available in dist_mat_expr!")

        if pae_cutoff_method == "truncate_on_pairwise_pae_with_center":
            # Find first neighbor where PAE(center, neighbor) > max_pae
            # Each residue has a 'pae' field which is PAE from center to that residue
            pae_cutoff_idx = hl.enumerate(dist_mat_expr).find(
                lambda x: hl.is_defined(x[1].pae) & (x[1].pae > max_pae)
            )
            if pae_cutoff_idx is not None:
                dist_mat_expr = hl.if_else(
                    hl.is_defined(pae_cutoff_idx),
                    dist_mat_expr[: pae_cutoff_idx[0]],
                    dist_mat_expr,
                )
        elif pae_cutoff_method == "filter_on_pairwise_pae_with_center":
            # Keep only neighbors with PAE(center, neighbor) <= max_pae
            # Each residue has a 'pae' field which is PAE from center to that residue
            dist_mat_expr = dist_mat_expr.filter(
                lambda x: ~hl.is_defined(x.pae) | (x.pae <= max_pae)
            )
        elif pae_cutoff_method == "filter_on_pairwise_pae_in_region":
            # Debug: Show PAE matrix before filtering (pLDDT filtering has already
            # been applied)
            if debug:
                debug_dict = debug_print_pae_matrix_for_region(
                    dist_mat_expr,
                    center_residue_index_expr,
                    max_pae=max_pae,
                    title="get_3d_residue: PAE matrix before filter_on_pairwise_pae_in_region (after pLDDT filtering)",
                    min_plddt=min_plddt,
                    plddt_cutoff_method=plddt_cutoff_method,
                )
                # Store debug output
                for center_idx, output_str in debug_dict.items():
                    if center_idx not in debug_outputs_by_residue:
                        debug_outputs_by_residue[center_idx] = {}
                    debug_outputs_by_residue[center_idx][
                        "PAE matrix before filter"
                    ] = output_str

            # Build region by scanning: for each residue, check if its PAE to any residue
            # already in the region exceeds max_pae. If not, add it to the region.
            # Note: pLDDT filtering has already been applied, so we're only considering
            # high-confidence residues (or all residues if exclude_low_plddt_from_stats was used)
            # Check if pae_array field exists (for pairwise PAE)
            has_pae_array_field = "pae_array" in dist_mat_expr.dtype.element_type.fields
            if not has_pae_array_field:
                raise ValueError("PAE array data is not available in dist_mat_expr!")
            # Check PAE from new residue to each residue in region
            # residue.pae_array[x.residue_index] gives PAE from residue to x
            # Logic: add residue if (pae_array is None) OR (NOT any(PAE > max_pae to residues in region))
            # This means: add if pae_array is None OR if all PAE values to residues in
            # region are <= max_pae
            dist_mat_expr = dist_mat_expr.scan(
                lambda region, residue: hl.if_else(
                    (residue.pae_array is None)
                    | ~hl.any(
                        region.map(
                            lambda x: hl.is_defined(residue.pae_array)
                            & (residue.pae_array[x.residue_index] > max_pae)
                        )
                    ),
                    region.append(residue),
                    region,
                ),
                hl.empty_array(dist_mat_expr.dtype.element_type),
            )

            # Get the final region (last element of the scan result)
            dist_mat_expr = dist_mat_expr[-1]

        else:
            raise ValueError(f"Unknown pae_cutoff_method: {pae_cutoff_method}")

    # Save excluded_residues dict before dropping fields (if exclude_from_stats exists)
    excluded_residues_dict_expr = None
    if "exclude_from_stats" in dist_mat_expr.dtype.element_type.fields:
        # Create a dict mapping residue_index to exclude_from_stats flag
        # Only include residues that are actually excluded (exclude_from_stats == True)
        excluded_residues_dict_expr = hl.dict(
            dist_mat_expr.filter(lambda x: x.exclude_from_stats).map(
                lambda x: (x.residue_index, hl.bool(True))
            )
        )

    dist_mat_expr = dist_mat_expr.map(lambda x: x.drop(*drop_fields))

    # Annotate neighbor observed and expected, cumulative observed and expected, and
    # upper bound of OE confidence interval.
    if "exclude_from_stats" in dist_mat_expr.dtype.element_type.fields:
        oe_expr = dist_mat_expr.map(
            lambda x: x.annotate(
                **hl.or_missing(~x.exclude_from_stats, oe_expr[x.residue_index])
            )
        )
    else:
        oe_expr = dist_mat_expr.map(lambda x: x.annotate(**oe_expr[x.residue_index]))

    oe_expr = get_cumulative_oe(oe_expr)
    oe_expr = calculate_oe_upper(oe_expr, alpha=alpha, oe_upper_method=oe_upper_method)

    # Debug: Show OE calculations after calculate_oe_upper
    if debug:
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: After calculate_oe_upper",
            max_pae=max_pae,
            min_plddt=min_plddt,
            oe_expr=oe_expr,
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx][
                "After calculate_oe_upper"
            ] = output_str

    # Print all debug outputs grouped by residue
    if debug and debug_outputs_by_residue:
        debug_string = ""
        # First, handle the full matrix case (key -1) if it exists
        if -1 in debug_outputs_by_residue:
            debug_string += debug_outputs_by_residue[-1]["Before sorting"]
            del debug_outputs_by_residue[-1]

        # Then print per-residue outputs
        for center_idx in sorted(debug_outputs_by_residue.keys()):
            residue_outputs = debug_outputs_by_residue[center_idx]

            # Print center residue index once before all steps.
            combined_output = f"\n{BOLD}Center residue index: {center_idx}{RESET}\n"

            # Combine all steps for this residue.
            for step_name in [
                "Before sorting",
                "After sorting",
                "PAE matrix before filter",
                "After calculate_oe_upper",
            ]:
                if step_name in residue_outputs:
                    combined_output += residue_outputs[step_name] + "\n"
            if combined_output:
                debug_string += combined_output

        logger.info(debug_string)

    # Get the 3D region with the lowest upper bound of the OE confidence interval for
    # each residue.
    min_moeuf_expr = get_min_oe_upper(oe_expr, min_exp_mis=min_exp_mis)

    # If excluded_residues dict was created, add it to the return value
    if excluded_residues_dict_expr is not None:
        min_moeuf_expr = min_moeuf_expr.annotate(
            excluded_residues=excluded_residues_dict_expr
        )

    return min_moeuf_expr


def determine_regions_with_min_oe_upper(
    af2_ht: hl.Table,
    oe_codon_ht: hl.Table,
    pae_ht: Optional[hl.Table] = None,
    plddt_ht: Optional[hl.Table] = None,
    alpha: float = 0.05,
    oe_upper_method: str = "gamma",
    min_exp_mis: int = MIN_EXP_MIS,
    max_pae: Optional[float] = 15,
    pae_cutoff_method: str = "truncate_on_pairwise_pae_with_center",
    min_plddt: Optional[float] = None,
    plddt_cutoff_method: Optional[str] = None,
    debug: bool = False,
) -> hl.Table:
    """
    Determine the most intolerant region for each UniProt ID and residue index.

    :param af2_ht: Hail Table with AlphaFold2 data.
    :param oe_codon_ht: Hail Table with observed and expected values for codons.
    :param pae_ht: Optional Hail Table with PAE (Predicted Aligned Error) data.
        Default is None.
    :param plddt_ht: Optional Hail Table with pLDDT (per-residue confidence) data.
        Default is None.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :param oe_upper_method: Method to use for calculating OE upper bound. Default is "gamma".
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is MIN_EXP_MIS.
    :param max_pae: Maximum allowed PAE for including residues in the region. Default is 15.
    :param pae_cutoff_method: Strategy for handling residues with high PAE values. Options:
        - "truncate_on_pairwise_pae_with_center": Remove all residues after the first one
          with PAE(center, neighbor) > max_pae (sequential cutoff).
        - "filter_on_pairwise_pae_with_center": Remove only residues with
          PAE(center, neighbor) > max_pae, keep others (filter individual residues).
        - "filter_on_pairwise_pae_in_region": Filter out residues whose maximum pairwise
          PAE to any residue in the region exceeds max_pae. Requires 2D pae_expr.
        Default is "truncate_on_pairwise_pae_with_center".
    :param min_plddt: Minimum allowed pLDDT for including residues in the region.
        Default is None.
    :param plddt_cutoff_method: Strategy for handling residues with low pLDDT values. Options:
        - "truncate_at_first_low_plddt": Remove all residues after the first one with
          pLDDT < min_plddt (sequential cutoff).
        - "remove_low_plddt_residues": Remove only residues with pLDDT < min_plddt,
          keep others (filter individual residues).
        - "mask_low_confidence_plddt": Mask (set to missing) residues with pLDDT < min_plddt
          instead of removing them.
        - "exclude_low_plddt_from_stats": Keep low pLDDT residues in regions (for assignment
          based on distance) but exclude them from statistical calculations (Obs, Exp, OE,
          OE upper, nLL, etc.).
        Default is None (no pLDDT filtering).
    :return: Hail Table with the most intolerant region for each UniProt ID and residue
        index
    """
    ann_expr = {"oe": oe_codon_ht[af2_ht.uniprot_id].oe_by_transcript}
    if pae_ht is not None:
        # For pairwise PAE in region, we need each residue's PAE array, not just the center's
        # Create an array of structs (residue_index, pae_array) for lookup
        _pae_lookup_ht = pae_ht.group_by("uniprot_id").aggregate(
            pae_lookup=hl.agg.collect(
                hl.struct(residue_index=pae_ht.aa_index, pae_array=pae_ht.pae)
            )
        )
        pae_lookup_array = _pae_lookup_ht[af2_ht.uniprot_id].pae_lookup
        # Annotate pae_lookup onto the HT for debugging
        ann_expr["pae_lookup"] = pae_lookup_array
        # Get center's PAE array for single PAE value
        center_pae_array = pae_ht[af2_ht.uniprot_id, af2_ht.aa_index].pae
        # Annotate each element of dist_mat with its corresponding PAE array
        ann_expr["dist_mat"] = af2_ht.dist_mat.map(
            lambda x: x.annotate(
                pae_array=hl.or_missing(
                    hl.is_defined(pae_lookup_array),
                    pae_lookup_array.find(
                        lambda y: y.residue_index == x.residue_index
                    ).pae_array,
                ),
                pae=hl.or_missing(
                    hl.is_defined(center_pae_array), center_pae_array[x.residue_index]
                ),  # Keep single PAE value (from center to this residue) for backward compatibility
            )
        )

    if plddt_ht is not None:
        # Create lookup for pLDDT values (for debugging)
        _plddt_lookup_ht = plddt_ht.group_by("uniprot_id").aggregate(
            plddt_lookup=hl.agg.collect(
                hl.struct(residue_index=plddt_ht.aa_index, plddt=plddt_ht.plddt)
            ),
        )
        ann_expr["plddt_lookup"] = _plddt_lookup_ht[af2_ht.uniprot_id].plddt_lookup

        # Add pLDDT to dist_mat (similar to how PAE is added)
        plddt_lookup_array = _plddt_lookup_ht[af2_ht.uniprot_id].plddt_lookup
        # Annotate each element of dist_mat with its corresponding pLDDT value
        # dist_mat may have already been created with PAE, or we need to create it
        if "dist_mat" in ann_expr:
            # dist_mat already exists (from PAE annotation), add pLDDT to it
            ann_expr["dist_mat"] = ann_expr["dist_mat"].map(
                lambda x: x.annotate(
                    plddt=hl.or_missing(
                        hl.is_defined(plddt_lookup_array),
                        plddt_lookup_array.find(
                            lambda y: y.residue_index == x.residue_index
                        ).plddt,
                    )
                )
            )
        else:
            # dist_mat hasn't been created yet (no PAE), so create it with pLDDT
            ann_expr["dist_mat"] = af2_ht.dist_mat.map(
                lambda x: x.annotate(
                    plddt=hl.or_missing(
                        hl.is_defined(plddt_lookup_array),
                        plddt_lookup_array.find(
                            lambda y: y.residue_index == x.residue_index
                        ).plddt,
                    )
                )
            )

    ht = af2_ht.annotate(**ann_expr)
    ht = ht.explode(ht.oe)
    # After explode, ht.oe is a struct with 'enst' and 'oe' (array)
    # Unpack the struct fields - this adds 'enst' and 'oe' (array) as top-level fields
    # But ht.oe still refers to the struct, so we need to drop it and use the
    # unpacked 'oe' field
    ht = ht.annotate(oe_array=ht.oe.oe, enst=ht.oe.enst).drop("oe")
    ht = ht.rename({"oe_array": "oe"})

    ht = ht.annotate(
        transcript_id=ht.enst,
        min_oe_upper=get_3d_residue(
            ht.aa_index,
            ht.dist_mat,
            ht.oe,
            alpha=alpha,
            oe_upper_method=oe_upper_method,
            min_exp_mis=min_exp_mis,
            pae_cutoff_method=pae_cutoff_method,
            plddt_cutoff_method=plddt_cutoff_method,
            debug=debug,
            max_pae=max_pae if pae_ht is not None else None,
            min_plddt=min_plddt if plddt_ht is not None else None,
        ),
    ).drop("enst")

    ht = ht.group_by("uniprot_id", "transcript_id").aggregate(
        oe=hl.agg.take(ht.oe, 1)[0],
        min_oe_upper=hl.agg.filter(
            hl.is_defined(ht.min_oe_upper),
            hl.agg.collect(
                ht.min_oe_upper.annotate(residue_index=ht.aa_index).drop("dist")
            ),
        ),
    )
    # TODO: minimize on oe or oe_upper?
    ht = ht.annotate(
        oe=hl.enumerate(ht.oe).map(lambda x: x[1].annotate(residue_index=x[0])),
        min_oe_upper=hl.sorted(ht.min_oe_upper, key=lambda x: x.oe),
    )

    # Debug: Show OE and regions for first uniprot/transcript
    if debug:
        debug_print_oe_and_regions(
            ht,
            title="determine_regions_with_min_oe_upper: Final output",
        )

    return ht


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


def prep_region_struct(region_expr, oe_expr, excluded_residues=None):
    """
    Prepare a region struct.

    :param region_expr: Region expression.
    :param oe_expr: OE expression.
    :param excluded_residues: Optional dict mapping residue_index to exclude_from_stats flag.
        If provided, excluded residues will be filtered out from statistical calculations.
    :return: Region struct expression.
    """
    # Filter out excluded residues from stats if excluded_residues dict is provided
    # We need to filter region_expr first, then annotate with OE
    if excluded_residues is not None:
        region_expr_for_stats = region_expr.filter(
            lambda res_idx: hl.if_else(
                hl.is_defined(excluded_residues.get(res_idx)),
                ~excluded_residues.get(res_idx),  # If in dict, exclude if True
                True,  # If not in dict, include it
            )
        )
    else:
        region_expr_for_stats = region_expr

    oe_expr_for_stats = annotate_region_with_oe(region_expr_for_stats, oe_expr)

    oe_agg_expr = get_agg_oe_for_region(oe_expr_for_stats)
    nll_expr = calculate_neg_log_likelihood(oe_expr_for_stats, oe_agg_expr.oe)
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


def debug_print_forward_round(
    ht: hl.Table,
    uniprot_id: str,
    transcript_id: str,
    model_comparison_method: str,
    title: str = "Forward Algorithm Round",
) -> str:
    """
    Print debug output for a round of the forward algorithm.

    :param ht: Hail Table with forward algorithm state
    :param uniprot_id: UniProt ID
    :param transcript_id: Transcript ID
    :param model_comparison_method: Model comparison method being used
    :param title: Title for the debug output
    """
    output_string = "\n"

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter(
        (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    output_string += f"{BOLD}    === {title} ==={RESET}\n"

    # Collect data
    # Note: key fields (idx, uniprot_id, transcript_id) are automatically included in collect()
    # Check if excluded_residues is available
    has_excluded_res_debug = "excluded_residues" in _ht_debug.row.dtype.fields
    debug_data = _ht_debug.select(
        "region",
        "selected",
        "selected_nll",
        "null_model",
        "best_aic",
        "found_best",
        "_region",
        "_updated_null",
        "oe",
        "center_residue_index",
        *["excluded_residues"] if has_excluded_res_debug else [],
    ).collect()

    if not debug_data:
        return ""

    # Print current state (use first row for state info, should be same for all)
    row = debug_data[0] if debug_data else None
    if row is None:
        return ""

    output_string += f"        {BOLD}Current State:{RESET}\n"
    output_string += (
        f"            Selected regions: {len(row.selected) if row.selected else 0}\n"
    )
    if row.selected_nll is None:
        selected_nll_str = "NA"
    elif row.selected_nll == 0:
        selected_nll_str = "NA (set to 0.0 as initial state)"
    else:
        selected_nll_str = f"{row.selected_nll:.4f}"
    output_string += f"            Selected NLL: {selected_nll_str}\n"
    best_aic_str = f"{row.best_aic:.4f}" if row.best_aic is not None else "NA"
    output_string += f"            Best AIC: {best_aic_str}\n"
    output_string += f"            Found best: {row.found_best}\n"

    if row.null_model:
        output_string += f"\n        {BOLD}Null Model (catch-all):{RESET}\n"
        obs_str = (
            f"{row.null_model.obs:7d}" if row.null_model.obs is not None else "     NA"
        )
        output_string += f"            Observed: {obs_str}\n"
        exp_str = (
            f"{row.null_model.exp:7.2f}"
            if row.null_model.exp is not None
            else "     NA"
        )
        output_string += f"            Expected: {exp_str}\n"
        # Color aggregate O/E values
        if row.null_model.oe is not None:
            oe_str = f"{row.null_model.oe:7.2f}"
        else:
            oe_str = "     NA"
        output_string += f"            O/E:      {oe_str}\n"
        # Color NLL values (more negative = better fit)
        if row.null_model.nll is not None:
            nll_str = f"{BOLD}{UNDERLINE}{row.null_model.nll:.4f}{RESET}"
        else:
            nll_str = "NA"
        output_string += f"            NLL:      {nll_str}\n\n"
        if hasattr(row.null_model, "region") and row.null_model.region:
            region_str = ", ".join([f"{r:7d}" for r in row.null_model.region])
            output_string += f"            Region residues: {region_str}\n"

            # Print per-residue observed and expected values for null model
            # Get oe from first row (should be same for all)
            if debug_data and len(debug_data) > 0:
                first_row = debug_data[0]
                if hasattr(first_row, "oe") and first_row.oe:
                    obs_list = []
                    exp_list = []
                    # Get excluded_residues if available
                    excluded_residues_dict = None
                    if (
                        hasattr(first_row, "excluded_residues")
                        and first_row.excluded_residues is not None
                    ):
                        excluded_residues_dict = first_row.excluded_residues
                    for res_idx in row.null_model.region:
                        # Check if this residue is excluded from stats
                        is_excluded = False
                        if excluded_residues_dict is not None:
                            if hasattr(excluded_residues_dict, "get"):
                                exclude_flag = excluded_residues_dict.get(res_idx)
                                is_excluded = (
                                    exclude_flag is True
                                    if exclude_flag is not None
                                    else False
                                )
                            elif isinstance(excluded_residues_dict, dict):
                                is_excluded = (
                                    excluded_residues_dict.get(res_idx, False) is True
                                )

                        if is_excluded:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                        elif res_idx < len(first_row.oe):
                            oe_entry = first_row.oe[res_idx]
                            obs_val = (
                                oe_entry.obs
                                if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                                else None
                            )
                            exp_val = (
                                oe_entry.exp
                                if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                                else None
                            )

                            obs_list.append(
                                f"{obs_val:7d}" if obs_val is not None else "    NA"
                            )
                            exp_list.append(
                                f"{exp_val:7.2f}" if exp_val is not None else "    NA"
                            )
                        else:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                    output_string += f"            Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                    output_string += f"            Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

    # Print all candidate regions being evaluated
    # First, filter to only valid candidates (those with valid _region data)
    valid_candidates = []
    for row in debug_data:
        if (
            hasattr(row, "_region")
            and row._region
            and hasattr(row._region, "obs")
            and row._region.obs is not None
        ):
            valid_candidates.append(row)

    # Find the minimum total NLL to highlight the best candidate
    min_total_nll = None
    min_total_nll_idx = None
    for row_idx, row in enumerate(valid_candidates):
        if (
            hasattr(row, "_region")
            and row._region
            and hasattr(row._region, "nll")
            and row._region.nll is not None
        ):
            if min_total_nll is None or row._region.nll < min_total_nll:
                min_total_nll = row._region.nll
                min_total_nll_idx = row_idx

    output_string += (
        f"\n        {BOLD}Candidate Regions ({len(valid_candidates)} total):{RESET}\n"
    )
    for row_idx, row in enumerate(valid_candidates):
        if hasattr(row, "_region") and row._region:
            idx_str = str(row.idx) if row.idx is not None else "NA"
            center_res_str = (
                f" (center residue: {row.center_residue_index})"
                if hasattr(row, "center_residue_index")
                and row.center_residue_index is not None
                else ""
            )
            output_string += (
                f"\n            {BOLD}Candidate {row_idx + 1}{center_res_str}:{RESET}\n"
            )
            obs_str = (
                f"{row._region.obs:7d}"
                if (hasattr(row._region, "obs") and row._region.obs is not None)
                else "     NA"
            )
            output_string += f"                Observed: {obs_str}\n"
            exp_str = (
                f"{row._region.exp:7.2f}"
                if (hasattr(row._region, "exp") and row._region.exp is not None)
                else "     NA"
            )
            output_string += f"                Expected: {exp_str}\n"
            # Keep O/E uncolored for candidates (selection is based on NLL, not O/E)
            oe_str = (
                f"{row._region.oe:7.2f}"
                if (hasattr(row._region, "oe") and row._region.oe is not None)
                else "     NA"
            )
            output_string += f"                O/E:      {oe_str}\n"
            # Color NLL values (more negative = better fit)
            if (
                hasattr(row._region, "region_nll")
                and row._region.region_nll is not None
            ):
                region_nll_str = f"{row._region.region_nll:.4f}"
            else:
                region_nll_str = "NA"
            output_string += f"                Region NLL: {region_nll_str}\n"
            # Highlight the best (lowest) total NLL
            if hasattr(row._region, "nll") and row._region.nll is not None:
                if min_total_nll_idx is not None and row_idx == min_total_nll_idx:
                    total_nll_str = f"{HIGHLIGHT}{row._region.nll:.4f} (BEST){RESET}"
                else:
                    total_nll_str = f"{BOLD}{UNDERLINE}{row._region.nll:.4f}{RESET}"
            else:
                total_nll_str = "NA"
            output_string += f"                {BOLD}{UNDERLINE}Total NLL: {total_nll_str}{RESET}\n\n"
            if hasattr(row._region, "region") and row._region.region:
                region_str = ", ".join([f"{r:7d}" for r in row._region.region])
                output_string += f"                Region residues: {region_str}\n"

                # Print per-residue observed and expected values
                if hasattr(row, "oe") and row.oe:
                    obs_list = []
                    exp_list = []
                    # Get excluded_residues if available
                    excluded_residues_dict = None
                    if (
                        hasattr(row, "excluded_residues")
                        and row.excluded_residues is not None
                    ):
                        excluded_residues_dict = row.excluded_residues
                    for res_idx in row._region.region:
                        # Check if this residue is excluded from stats
                        is_excluded = False
                        if excluded_residues_dict is not None:
                            # excluded_residues_dict is a dict-like structure
                            # Check if res_idx is in the dict and if its value is True
                            if hasattr(excluded_residues_dict, "get"):
                                exclude_flag = excluded_residues_dict.get(res_idx)
                                is_excluded = (
                                    exclude_flag is True
                                    if exclude_flag is not None
                                    else False
                                )
                            elif isinstance(excluded_residues_dict, dict):
                                is_excluded = (
                                    excluded_residues_dict.get(res_idx, False) is True
                                )

                        if is_excluded:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                        elif res_idx < len(row.oe):
                            oe_entry = row.oe[res_idx]
                            obs_val = (
                                oe_entry.obs
                                if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                                else None
                            )
                            exp_val = (
                                oe_entry.exp
                                if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                                else None
                            )

                            obs_list.append(
                                f"{obs_val:7d}" if obs_val is not None else "    NA"
                            )
                            exp_list.append(
                                f"{exp_val:7.2f}" if exp_val is not None else "    NA"
                            )
                        else:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                    output_string += f"                Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                    output_string += f"                Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

            # Print updated null model after removing candidate region
            if hasattr(row, "_updated_null") and row._updated_null:
                output_string += f"\n                {BOLD}Updated Null Model:{RESET}\n"
                obs_str = (
                    f"{row._updated_null.obs:7d}"
                    if (
                        hasattr(row._updated_null, "obs")
                        and row._updated_null.obs is not None
                    )
                    else "     NA"
                )
                output_string += f"                    Observed: {obs_str}\n"
                exp_str = (
                    f"{row._updated_null.exp:7.2f}"
                    if (
                        hasattr(row._updated_null, "exp")
                        and row._updated_null.exp is not None
                    )
                    else "     NA"
                )
                output_string += f"                    Expected: {exp_str}\n"
                # Color aggregate O/E
                if (
                    hasattr(row._updated_null, "oe")
                    and row._updated_null.oe is not None
                ):
                    oe_str = f"{row._updated_null.oe:7.2f}"
                else:
                    oe_str = "     NA"
                output_string += f"                    O/E:      {oe_str}\n"
                # Color NLL values
                if (
                    hasattr(row._updated_null, "nll")
                    and row._updated_null.nll is not None
                ):
                    nll_str = f"{row._updated_null.nll:.4f}"
                else:
                    nll_str = "NA"
                output_string += f"                    NLL:      {nll_str}\n\n"
                if hasattr(row._updated_null, "region") and row._updated_null.region:
                    region_str = ", ".join(
                        [f"{r:7d}" for r in row._updated_null.region]
                    )
                    output_string += (
                        f"                    Region residues: {region_str}\n"
                    )

                    # Print per-residue observed and expected values for updated null
                    # model
                    if hasattr(row, "oe") and row.oe:
                        obs_list = []
                        exp_list = []
                        # Get excluded_residues if available
                        excluded_residues_dict = None
                        if (
                            hasattr(row, "excluded_residues")
                            and row.excluded_residues is not None
                        ):
                            excluded_residues_dict = row.excluded_residues
                        for res_idx in row._updated_null.region:
                            # Check if this residue is excluded from stats
                            is_excluded = False
                            if excluded_residues_dict is not None:
                                if hasattr(excluded_residues_dict, "get"):
                                    exclude_flag = excluded_residues_dict.get(res_idx)
                                    is_excluded = (
                                        exclude_flag is True
                                        if exclude_flag is not None
                                        else False
                                    )
                                elif isinstance(excluded_residues_dict, dict):
                                    is_excluded = (
                                        excluded_residues_dict.get(res_idx, False)
                                        is True
                                    )

                            if is_excluded:
                                obs_list.append("     NA")
                                exp_list.append("     NA")
                            elif res_idx < len(row.oe):
                                oe_entry = row.oe[res_idx]
                                obs_val = (
                                    oe_entry.obs
                                    if hasattr(oe_entry, "obs")
                                    and oe_entry.obs is not None
                                    else None
                                )
                                exp_val = (
                                    oe_entry.exp
                                    if hasattr(oe_entry, "exp")
                                    and oe_entry.exp is not None
                                    else None
                                )

                                obs_list.append(
                                    f"{obs_val:7d}" if obs_val is not None else "    NA"
                                )
                                exp_list.append(
                                    f"{exp_val:7.2f}"
                                    if exp_val is not None
                                    else "    NA"
                                )
                            else:
                                obs_list.append("     NA")
                                exp_list.append("     NA")
                        output_string += f"                    Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                        output_string += f"                    Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

    # Print model comparison info if available (from first row)
    if hasattr(row, "aic_cand") and row.aic_cand is not None:
        output_string += (
            f"\n    {BOLD}Model Comparison ({model_comparison_method}):{RESET}\n"
        )
        output_string += f"      Candidate AIC: {row.aic_cand:.4f}\n"
        if hasattr(row, "p_lrt") and row.p_lrt is not None:
            output_string += f"      LRT p-value: {row.p_lrt:.6f}\n"
        if hasattr(row, "w_cand") and row.w_cand is not None:
            output_string += f"      AIC weight: {row.w_cand:.4f}\n"

    return output_string


def run_forward(
    ht,
    min_exp_mis=MIN_EXP_MIS,
    oe_upper_method: str = "gamma",
    model_comparison_method: str = "aic",
    lrt_alpha: float = 0.001,
    lrt_df_added: int = 1,
    bonferroni_per_round: bool = True,
    aic_weight_thresh: float = 0.80,
    debug: bool = False,
):
    """
    Run the forward algorithm to find the most intolerant region.


    model_comparison_method:
      - "aic"         : choose region if AIC(candidate) < AIC(current)
      - "aic_weight"  : choose region if Akaike weight(candidate) >= aic_weight_thresh
      - "lrt"         : choose region if LRT p-value <= (alpha [Bonferroni-adjusted per round if enabled])

    :param ht: Hail Table with the most intolerant region for each UniProt ID and residue
        index.
    :param min_exp_mis: Minimum number of expected missense variants in a region to be
        considered for constraint calculation. Default is MIN_EXP_MIS.
    :param oe_upper_method: Method to use for calculating the upper bound of the OE
        confidence interval. Default is "gamma".
    :param model_comparison_method: Method to use for comparing the candidate region to
        the current region. Options are "aic", "aic_weight", and "lrt". Default is
        "aic".
    :param lrt_alpha: Significance level for the LRT. Default is 0.001.
    :param lrt_df_added: Degrees of freedom added per region for the LRT. Default is 1.
    :param bonferroni_per_round: Whether to adjust the alpha level by the number of
        candidates scanned this round. Default is True.
    :param aic_weight_thresh: Threshold for the Akaike weight. Default is 0.80.
    :return: Hail Table annotated with the observed and expected values for each residue.
    """
    if oe_upper_method not in ["gamma", "chisq"]:
        raise ValueError(f"Invalid OE upper method: {oe_upper_method}")
    if model_comparison_method not in {"aic", "aic_weight", "lrt"}:
        raise ValueError(f"Invalid model comparison method: {model_comparison_method}")

    num_residues = ht.oe.length()
    null_region = hl.range(num_residues)

    # Get excluded_residues dict if present (for pLDDT exclusion method)
    has_excluded_res = "excluded_residues" in ht.min_oe_upper.dtype.element_type.fields
    excluded_residues_global = None
    if has_excluded_res:
        # Get excluded_residues from first region (should be same for all)
        excluded_residues_global = ht.min_oe_upper[0].excluded_residues

    null_model = prep_region_struct(
        null_region, ht.oe, excluded_residues=excluded_residues_global
    )
    ht = ht.select(
        "oe",
        num_residues=num_residues,
        null_model=null_model,
        regions=hl.enumerate(
            ht.min_oe_upper.map(
                lambda x: x.select(
                    "region",
                    "residue_index",
                    *["excluded_residues"] if has_excluded_res else [],
                )
            ).filter(lambda x: x.region.length() < num_residues)
        ),
        selected=hl.empty_array(null_model.dtype),
        selected_nll=0.0,
        best_aic=getAIC(null_model, null_model.nll),
        found_best=False,
    )

    ht = ht.explode("regions")
    ht = ht.transmute(
        idx=ht.regions[0],
        region=ht.regions[1].region,
        center_residue_index=ht.regions[1].residue_index,
        **(
            {"excluded_residues": ht.regions[1].get("excluded_residues")}
            if has_excluded_res
            else {}
        ),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "idx")
    # ht = ht.repartition(200, shuffle=True)  # ht = ht.repartition(50, shuffle=True)
    ht = ht.repartition(1, shuffle=True)
    ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_explode", "ht"))

    # Collect debug output to print at the end
    debug_outputs = []

    round_num = 1
    while ht.aggregate(hl.agg.any(hl.is_defined(ht.region))):
        # For each region in regions, update the list of selected by
        # removing the residues in the region from the "catch all remaining"
        # region, and adding the new region to the selected list.
        region_expr = prep_region_struct(
            ht.region,
            ht.oe,
            excluded_residues=(ht.excluded_residues if has_excluded_res else None),
        )
        ht = ht.annotate(_region=region_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.prep1", "ht")
        )
        region_expr = ht._region
        updated_null_expr = remove_residues_from_region(ht.null_model, region_expr)
        ht = ht.annotate(_updated_null=updated_null_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.prep2", "ht")
        )
        updated_null_expr = prep_region_struct(
            ht._updated_null.region,
            ht.oe,
            excluded_residues=(ht.excluded_residues if has_excluded_res else None),
        )
        ht = ht.annotate(_updated_null=updated_null_expr).checkpoint(
            hl.utils.new_temp_file(f"forward_round_{round_num}.prep3", "ht")
        )
        updated_null_expr = ht._updated_null
        region_expr = ht._region
        region_expr = region_expr.annotate(
            null_model=updated_null_expr,
            region_nll=region_expr.nll,
            nll=updated_null_expr.nll + ht.selected_nll + region_expr.nll,
        )

        if debug:
            debug_outputs.append(
                f"\n\n{BOLD}=== run_forward: Round {round_num} ==={RESET}\n"
            )
            # Get first key for debugging
            uniprot_id = None
            transcript_id = None
            first_row = ht.head(1)
            if first_row.count() > 0:
                uniprot_id = first_row.uniprot_id.collect()[0]
                transcript_id = first_row.transcript_id.collect()[0]
            debug_outputs.append(
                f"\nUniProt ID: {uniprot_id}, Transcript ID: {transcript_id}\n"
            )

            # Debug after preparing candidate region
            debug_output = debug_print_forward_round(
                ht.annotate(_region=region_expr),
                uniprot_id,
                transcript_id,
                model_comparison_method,
                title=f"run_forward: Round {round_num} - After preparing candidate region",
            )
            if debug_output:
                debug_outputs.append(debug_output)

        # Scan candidates this round; keep those with enough expected missense.
        ht2 = ht.select(exp=region_expr.exp, nll=region_expr.nll)
        ht2 = ht2.filter(hl.is_defined(ht2.nll) & (ht2.exp >= min_exp_mis))

        # For each (uniprot, transcript), find argmin by NLL + count candidates.
        ht2 = (
            ht2.group_by("uniprot_id", "transcript_id")
            .aggregate(
                m_candidates=hl.agg.count(),
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
                ),
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

        best_region = ht2[ht.uniprot_id, ht.transcript_id].best_region
        m_candidates = hl.or_else(ht2[ht.uniprot_id, ht.transcript_id].m_candidates, 1)

        if debug:
            # Debug after finding best candidate
            _ht_debug = ht.head(1)
            if _ht_debug.count() > 0:
                uniprot_id = _ht_debug.uniprot_id.collect()[0]
                transcript_id = _ht_debug.transcript_id.collect()[0]
                _ht_debug = ht.filter(
                    (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
                )
                # Only print if this uniprot/transcript has rows with defined regions
                has_regions = _ht_debug.aggregate(
                    hl.agg.any(hl.is_defined(_ht_debug.region))
                )
                if has_regions:
                    best_region_debug = (
                        ht2.filter(
                            (ht2.uniprot_id == uniprot_id)
                            & (ht2.transcript_id == transcript_id)
                        )
                        .select("best_region", "m_candidates", "min_idx")
                        .collect()
                    )
                    _ht_debug_oe = _ht_debug.select("oe").collect()
                    # Get center_residue_index for the best candidate
                    center_residue_idx = None
                    if best_region_debug and best_region_debug[0].min_idx is not None:
                        _ht_center = (
                            _ht_debug.filter(
                                (_ht_debug.uniprot_id == uniprot_id)
                                & (_ht_debug.transcript_id == transcript_id)
                                & (_ht_debug.idx == best_region_debug[0].min_idx)
                            )
                            .select("center_residue_index")
                            .collect()
                        )
                        if _ht_center and len(_ht_center) > 0:
                            center_residue_idx = _ht_center[0].center_residue_index
                    if best_region_debug:
                        br_row = best_region_debug[0]
                        output_string = f"\n    {BOLD}=== run_forward: Round {round_num} - Best Candidate ==={RESET}\n\n"
                        output_string += f"        Number of candidates evaluated: {br_row.m_candidates}\n"
                        if br_row.best_region:
                            center_res_str = (
                                f" (center residue: {center_residue_idx})"
                                if center_residue_idx is not None
                                else ""
                            )
                            output_string += f"\n        {BOLD}Best Candidate Region{center_res_str}:{RESET}\n"
                            output_string += (
                                f"            Observed: {br_row.best_region.obs:7d}\n"
                            )
                            output_string += (
                                f"            Expected: {br_row.best_region.exp:7.2f}\n"
                            )
                            # Color aggregate O/E
                            if br_row.best_region.oe is not None:
                                oe_str = f"{br_row.best_region.oe:7.2f}"
                            else:
                                oe_str = "     NA"
                            output_string += f"            O/E:      {oe_str}\n"
                            # Color NLL values
                            if br_row.best_region.region_nll is not None:
                                region_nll_str = f"{br_row.best_region.region_nll:.4f}"
                            else:
                                region_nll_str = "NA"
                            output_string += (
                                f"            Region NLL: {region_nll_str}\n"
                            )
                            if br_row.best_region.nll is not None:
                                total_nll_str = f"{br_row.best_region.nll:.4f}"
                            else:
                                total_nll_str = "NA"
                            output_string += f"            {BOLD}{UNDERLINE}Total NLL: {total_nll_str}{RESET}\n\n"
                            if (
                                hasattr(br_row.best_region, "region")
                                and br_row.best_region.region
                            ):
                                region_str = ", ".join(
                                    [f"{r:7d}" for r in br_row.best_region.region]
                                )
                                output_string += (
                                    f"            Region residues: {region_str}\n"
                                )

                                # Print per-residue observed and expected values
                                if (
                                    _ht_debug_oe
                                    and len(_ht_debug_oe) > 0
                                    and _ht_debug_oe[0].oe
                                ):
                                    obs_list = []
                                    exp_list = []
                                    for res_idx in br_row.best_region.region:
                                        if res_idx < len(_ht_debug_oe[0].oe):
                                            oe_entry = _ht_debug_oe[0].oe[res_idx]
                                            obs_val = (
                                                oe_entry.obs
                                                if hasattr(oe_entry, "obs")
                                                and oe_entry.obs is not None
                                                else None
                                            )
                                            exp_val = (
                                                oe_entry.exp
                                                if hasattr(oe_entry, "exp")
                                                and oe_entry.exp is not None
                                                else None
                                            )

                                            obs_list.append(
                                                f"{obs_val:7d}"
                                                if obs_val is not None
                                                else "    NA"
                                            )
                                            exp_list.append(
                                                f"{exp_val:7.2f}"
                                                if exp_val is not None
                                                else "    NA"
                                            )
                                    output_string += f"            Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                                    output_string += f"            Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"
                        debug_outputs.append(output_string)

        # Totals for LRT (candidate only if best_region defined).
        current_nll = ht.selected_nll + ht.null_model.nll
        candidate_nll = hl.or_missing(hl.is_defined(best_region), best_region.nll)

        # Candidate AIC (only if defined).
        aic_cand = hl.or_missing(
            hl.is_defined(best_region),
            getAIC(best_region.null_model, 0)
            + getAIC(
                ht.selected.append(best_region.drop("null_model")), best_region.nll
            ),
        )

        # Build found_best as a fresh per-row expression.
        found_best_expr = ht.found_best | hl.is_missing(best_region)

        if model_comparison_method == "lrt":
            lrt_stat = hl.or_missing(
                hl.is_defined(best_region), 2 * (current_nll - candidate_nll)
            )
            p_lrt = hl.or_missing(
                hl.is_defined(best_region), hl.pchisqtail(lrt_stat, lrt_df_added)
            )
            adj_alpha = hl.if_else(
                bonferroni_per_round, lrt_alpha / hl.max(1, m_candidates), lrt_alpha
            )
            found_best_expr |= p_lrt > adj_alpha

        elif model_comparison_method == "aic":
            found_best_expr |= aic_cand >= ht.best_aic

        elif model_comparison_method == "aic_weight":
            w_cand = 1 / (1 + hl.exp(0.5 * (aic_cand - ht.best_aic)))
            found_best_expr |= w_cand < aic_weight_thresh

        if debug:
            # Debug model comparison after expressions are computed
            _ht_debug = ht.head(1)
            if _ht_debug.count() > 0:
                uniprot_id = _ht_debug.uniprot_id.collect()[0]
                transcript_id = _ht_debug.transcript_id.collect()[0]
                _ht_debug = ht.filter(
                    (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
                )
                # Only print if this uniprot/transcript has rows with defined regions
                has_regions = _ht_debug.aggregate(
                    hl.agg.any(hl.is_defined(_ht_debug.region))
                )
                if has_regions:
                    _ht_debug2 = ht2.filter(
                        (ht2.uniprot_id == uniprot_id)
                        & (ht2.transcript_id == transcript_id)
                    )
                    # Re-derive expressions from the filtered table
                    _best_region_debug = _ht_debug2[
                        (_ht_debug.uniprot_id, _ht_debug.transcript_id)
                    ].best_region
                    _m_candidates_debug = hl.or_else(
                        _ht_debug2[
                            (_ht_debug.uniprot_id, _ht_debug.transcript_id)
                        ].m_candidates,
                        1,
                    )
                    _current_nll_debug = (
                        _ht_debug.selected_nll + _ht_debug.null_model.nll
                    )
                    _candidate_nll_debug = hl.or_missing(
                        hl.is_defined(_best_region_debug), _best_region_debug.nll
                    )
                    _aic_cand_debug = hl.or_missing(
                        hl.is_defined(_best_region_debug),
                        getAIC(_best_region_debug.null_model, 0)
                        + getAIC(
                            _ht_debug.selected.append(
                                _best_region_debug.drop("null_model")
                            ),
                            _best_region_debug.nll,
                        ),
                    )
                    _found_best_debug = _ht_debug.found_best | hl.is_missing(
                        _best_region_debug
                    )
                    if model_comparison_method == "lrt":
                        _lrt_stat_debug = hl.or_missing(
                            hl.is_defined(_best_region_debug),
                            2 * (_current_nll_debug - _candidate_nll_debug),
                        )
                        _p_lrt_debug = hl.or_missing(
                            hl.is_defined(_best_region_debug),
                            hl.pchisqtail(_lrt_stat_debug, lrt_df_added),
                        )
                        _adj_alpha_debug = hl.if_else(
                            bonferroni_per_round,
                            lrt_alpha / hl.max(1, _m_candidates_debug),
                            lrt_alpha,
                        )
                        _found_best_debug |= _p_lrt_debug > _adj_alpha_debug
                    elif model_comparison_method == "aic":
                        _found_best_debug |= _aic_cand_debug >= _ht_debug.best_aic
                    elif model_comparison_method == "aic_weight":
                        _w_cand_debug = 1 / (
                            1 + hl.exp(0.5 * (_aic_cand_debug - _ht_debug.best_aic))
                        )
                        _found_best_debug |= _w_cand_debug < aic_weight_thresh

                    debug_data = (
                        _ht_debug.annotate(
                            current_nll=_current_nll_debug,
                            candidate_nll=_candidate_nll_debug,
                            aic_cand=_aic_cand_debug,
                            found_best=_found_best_debug,
                        )
                        .select(
                            "selected_nll",
                            "null_model",
                            "best_aic",
                            "current_nll",
                            "candidate_nll",
                            "aic_cand",
                            "found_best",
                        )
                        .collect()
                    )
                    debug_data2 = _ht_debug2.select(
                        "best_region",
                        "m_candidates",
                    ).collect()
                    if debug_data and debug_data2:
                        row = debug_data[0]
                        row2 = debug_data2[0]
                    output_string = f"\n{BOLD}    === run_forward: Round {round_num} - Model Comparison ==={RESET}\n"
                    output_string += f"        {BOLD}Current Model:{RESET}\n"
                    # Color NLL values
                    if row.current_nll is not None:
                        current_nll_str = (
                            f"{BOLD}{UNDERLINE}{row.current_nll:.4f}{RESET}"
                        )
                    else:
                        current_nll_str = "NA"
                    output_string += f"            Current NLL: {current_nll_str}\n"
                    # Color AIC values
                    if row.best_aic is not None:
                        best_aic_str = f"{row.best_aic:.4f}"
                    else:
                        best_aic_str = "NA"
                    output_string += f"            Best AIC: {best_aic_str}\n"
                    if row2.best_region and row.candidate_nll is not None:
                        output_string += f"\n          {BOLD}Candidate Model:{RESET}\n"
                        # Color candidate NLL
                        if row.candidate_nll is not None:
                            candidate_nll_str = (
                                f"{BOLD}{UNDERLINE}{row.candidate_nll:.4f}{RESET}"
                            )
                        else:
                            candidate_nll_str = "NA"
                        output_string += (
                            f"            Candidate NLL: {candidate_nll_str}\n"
                        )
                        if row.aic_cand is not None:
                            # Color candidate AIC - highlight if better (lower) than
                            # current
                            if row.best_aic is not None and row.aic_cand < row.best_aic:
                                aic_cand_str = f"{GREEN}{row.aic_cand:.4f}{RESET}"
                            else:
                                aic_cand_str = f"{row.aic_cand:.4f}"
                            output_string += (
                                f"            Candidate AIC: {aic_cand_str}\n"
                            )
                        if model_comparison_method == "lrt":
                            lrt_stat = 2 * (row.current_nll - row.candidate_nll)
                            # Use Hail's pchisqtail for consistency
                            p_lrt_val = hl.eval(hl.pchisqtail(lrt_stat, lrt_df_added))
                            adj_alpha = (
                                lrt_alpha / max(1, row2.m_candidates)
                                if bonferroni_per_round
                                else lrt_alpha
                            )
                            output_string += (
                                f"            LRT statistic: {lrt_stat:.4f}\n"
                            )
                            output_string += (
                                f"            LRT p-value: {p_lrt_val:.6f}\n"
                            )
                            output_string += (
                                f"            Adjusted alpha: {adj_alpha:.6f}\n"
                            )
                            # Color acceptance status
                            accept = p_lrt_val <= adj_alpha
                            if accept:
                                accept_str = f"{GREEN}{accept}{RESET}"
                            else:
                                accept_str = f"{RED}{accept}{RESET}"
                            output_string += (
                                f"            Accept candidate: {accept_str}\n"
                            )
                        elif model_comparison_method == "aic":
                            accept = (
                                row.aic_cand < row.best_aic
                                if row.aic_cand is not None
                                else False
                            )
                            if accept:
                                accept_str = f"{GREEN}{accept}{RESET}"
                            else:
                                accept_str = f"{RED}{accept}{RESET}"
                            output_string += (
                                f"            Accept candidate: {accept_str}\n"
                            )
                        elif model_comparison_method == "aic_weight":
                            if row.aic_cand is not None:
                                w_cand_val = 1 / (
                                    1 + np.exp(0.5 * (row.aic_cand - row.best_aic))
                                )
                                output_string += (
                                    f"            AIC weight: {w_cand_val:.4f}\n"
                                )
                                accept = w_cand_val >= aic_weight_thresh
                                if accept:
                                    accept_str = f"{GREEN}{accept}{RESET}"
                                else:
                                    accept_str = f"{RED}{accept}{RESET}"
                                output_string += (
                                    f"            Accept candidate: {accept_str}\n"
                                )
                    output_string += f"\n        Found best (stop): {row.found_best}\n"
                    debug_outputs.append(output_string)

        # Current (no-change) state.
        curr_vals = hl.struct(
            null_model=ht.null_model,
            selected=ht.selected,
            selected_nll=ht.selected_nll,
            best_aic=ht.best_aic,
            region=hl.missing(ht.region.dtype),
        )

        # Build updated state LAZILY so we don't touch best_region when stopping.
        def _updated_vals():
            next_region = remove_residues_from_region(
                hl.struct(region=ht.region), best_region
            ).region
            next_region = hl.or_missing(
                hl.is_defined(next_region) & (next_region.length() > 0), next_region
            )
            return hl.struct(
                null_model=best_region.null_model,
                selected=ht.selected.append(best_region.drop("null_model")),
                selected_nll=ht.selected_nll + best_region.region_nll,
                best_aic=aic_cand,
                region=next_region,
            )

        if debug:
            # Collect candidates BEFORE updating region (for showing removed candidates)
            _ht_debug_before_update = ht.head(1)
            if _ht_debug_before_update.count() > 0:
                uniprot_id_before = _ht_debug_before_update.uniprot_id.collect()[0]
                transcript_id_before = _ht_debug_before_update.transcript_id.collect()[
                    0
                ]
                _ht_candidates_before = ht.filter(
                    (ht.uniprot_id == uniprot_id_before)
                    & (ht.transcript_id == transcript_id_before)
                )
                # Collect region (candidate regions) before they're removed, along with
                # oe data
                candidates_before_filter = _ht_candidates_before.select(
                    "center_residue_index", "region", "oe"
                ).collect()
            else:
                candidates_before_filter = []
                uniprot_id_before = None
                transcript_id_before = None

        ht = ht.annotate(
            **hl.if_else(found_best_expr, curr_vals, _updated_vals()),
            found_best=found_best_expr,
        )

        # Keep only rows that still have remaining candidate regions or the idx==0.
        ht = ht.filter(hl.is_defined(ht.region) | (ht.idx == 0))

        if debug:
            # Debug after round completion
            _ht_debug = ht.head(1)
            if _ht_debug.count() > 0:
                uniprot_id = _ht_debug.uniprot_id.collect()[0]
                transcript_id = _ht_debug.transcript_id.collect()[0]
                _ht_debug = ht.filter(
                    (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
                )
                # Only print if this uniprot/transcript has rows with defined regions
                has_regions = _ht_debug.aggregate(
                    hl.agg.any(hl.is_defined(_ht_debug.region))
                )
                if has_regions:
                    debug_data = _ht_debug.select(
                        "selected", "selected_nll", "best_aic", "found_best", "region"
                    ).collect()
                    if debug_data:
                        row = debug_data[0]
                        output_string = f"\n\n\n{BOLD}    === run_forward: Round {round_num} - After Update ==={RESET}\n"
                        output_string += f"        Selected regions: {len(row.selected) if row.selected else 0}\n"
                        # Show the newly selected region (last one in the array)
                        if row.selected and len(row.selected) > 0:
                            latest_region = row.selected[-1]
                            output_string += (
                                f"\n        {BOLD}Newly Selected Region:{RESET}\n"
                            )
                            if (
                                hasattr(latest_region, "obs")
                                and latest_region.obs is not None
                            ):
                                output_string += (
                                    f"            Observed: {latest_region.obs:7d}\n"
                                )
                            if (
                                hasattr(latest_region, "exp")
                                and latest_region.exp is not None
                            ):
                                output_string += (
                                    f"            Expected: {latest_region.exp:7.2f}\n"
                                )
                            if (
                                hasattr(latest_region, "oe")
                                and latest_region.oe is not None
                            ):
                                output_string += (
                                    f"            O/E:      {latest_region.oe:7.2f}\n"
                                )
                            if (
                                hasattr(latest_region, "region")
                                and latest_region.region
                            ):
                                region_str = ", ".join(
                                    [f"{r:7d}" for r in latest_region.region]
                                )
                                output_string += f"            {BOLD}Residues ({len(latest_region.region)}): {region_str}{RESET}\n"
                        # Color selected NLL
                        if row.selected_nll is not None:
                            selected_nll_str = (
                                f"{BOLD}{UNDERLINE}{row.selected_nll:.4f}{RESET}"
                            )
                        else:
                            selected_nll_str = "NA"
                        output_string += (
                            f"\n            Selected NLL: {selected_nll_str}\n"
                        )
                        # Color best AIC
                        if row.best_aic is not None:
                            best_aic_str = f"{row.best_aic:.4f}"
                        else:
                            best_aic_str = "NA"
                        output_string += f"            Best AIC: {best_aic_str}\n"
                        output_string += f"            Found best: {row.found_best}\n"

                        # Always show all candidates from this round with chosen region
                        # residues colored red
                        if (
                            row.selected
                            and len(row.selected) > 0
                            and candidates_before_filter
                        ):
                            latest_region = row.selected[-1]
                            chosen_residues = (
                                set(latest_region.region)
                                if hasattr(latest_region, "region")
                                and latest_region.region
                                else set()
                            )
                            output_string += f"\n        {BOLD}All Candidates from This Round (with chosen region residues colored red):{RESET}\n"
                            for cand_idx, cand_row in enumerate(
                                candidates_before_filter
                            ):
                                if hasattr(cand_row, "region") and cand_row.region:
                                    center_res_str = (
                                        f" (center residue: {cand_row.center_residue_index})"
                                        if hasattr(cand_row, "center_residue_index")
                                        and cand_row.center_residue_index is not None
                                        else ""
                                    )
                                    output_string += f"\n            {BOLD}Candidate {cand_idx + 1}{center_res_str}:{RESET}\n"
                                    # Calculate obs/exp/oe from oe array for this
                                    # candidate's residues
                                    region_residues = cand_row.region
                                    if hasattr(cand_row, "oe") and cand_row.oe:
                                        total_obs = 0
                                        total_exp = 0.0
                                        for res_idx in region_residues:
                                            if res_idx < len(cand_row.oe):
                                                oe_entry = cand_row.oe[res_idx]
                                                if (
                                                    hasattr(oe_entry, "obs")
                                                    and oe_entry.obs is not None
                                                ):
                                                    total_obs += oe_entry.obs
                                                if (
                                                    hasattr(oe_entry, "exp")
                                                    and oe_entry.exp is not None
                                                ):
                                                    total_exp += oe_entry.exp
                                        if total_exp > 0:
                                            total_oe = total_obs / total_exp
                                        else:
                                            total_oe = None
                                        output_string += f"                Observed: {total_obs:7d}\n"
                                        output_string += f"                Expected: {total_exp:7.2f}\n"
                                        if total_oe is not None:
                                            output_string += f"                O/E:      {total_oe:7.2f}\n"
                                    # Show region residues with chosen region residues
                                    # colored red
                                    residue_strs = []
                                    for res in region_residues:
                                        if res in chosen_residues:
                                            residue_strs.append(f"{RED}{res}{RESET}")
                                        else:
                                            residue_strs.append(f"{res}")
                                    region_str = ", ".join(residue_strs)
                                    output_string += f"                {BOLD}Residues ({len(region_residues)}): {region_str}{RESET}\n"
                                    # Calculate remaining residues after removing chosen
                                    # region
                                    remaining_residues = [
                                        r
                                        for r in region_residues
                                        if r not in chosen_residues
                                    ]
                                    if remaining_residues:
                                        remaining_str = ", ".join(
                                            [f"{r}" for r in remaining_residues]
                                        )
                                        output_string += f"                Remaining after removal: {remaining_str} ({len(remaining_residues)} residues)\n"
                                    else:
                                        output_string += f"                {RED}(All residues removed - candidate eliminated){RESET}\n"
                        debug_outputs.append(output_string)

        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}", "ht"))
        round_num += 1

    if debug:
        # Debug final state before flattening
        _ht_debug = ht.head(1)
        if _ht_debug.count() > 0:
            uniprot_id = _ht_debug.uniprot_id.collect()[0]
            transcript_id = _ht_debug.transcript_id.collect()[0]
            _ht_debug = ht.filter(
                (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
            )
            debug_data = _ht_debug.select(
                "selected", "selected_nll", "best_aic", "null_model", "oe"
            ).collect()
            if debug_data:
                row = debug_data[0]
                output_string = f"\n\n\n{BOLD}=== run_forward: Final State ==={RESET}\n"
                output_string += f"    {BOLD}Total selected regions: {len(row.selected) if row.selected else 0}{RESET}\n"
                # Color final selected NLL
                if row.selected_nll is not None:
                    selected_nll_str = f"{BOLD}{UNDERLINE}{row.selected_nll:.4f}{RESET}"
                else:
                    selected_nll_str = "NA"
                output_string += (
                    f"    {BOLD}Final selected NLL:{RESET} {selected_nll_str}\n"
                )
                # Color final best AIC
                if row.best_aic is not None:
                    best_aic_str = f"{row.best_aic:.4f}"
                else:
                    best_aic_str = "NA"
                output_string += f"    {BOLD}Final best AIC:{RESET} {best_aic_str}\n"
                if row.selected:
                    output_string += f"\n    {BOLD}Selected Regions:{RESET}\n"
                    for idx, region in enumerate(row.selected):
                        output_string += f"        {BOLD}Region {idx}:{RESET}\n"
                        output_string += (
                            f"            {BOLD}Observed:{RESET} {region.obs:7d}\n"
                        )
                        output_string += (
                            f"            {BOLD}Expected:{RESET} {region.exp:7.2f}\n"
                        )
                        # Color aggregate O/E
                        if region.oe is not None:
                            oe_str = f"{region.oe:7.2f}"
                        else:
                            oe_str = "     NA"
                        output_string += f"            {BOLD}O/E: {oe_str}{RESET}\n"
                        if hasattr(region, "region") and region.region:
                            region_str = ", ".join([f"{r:7d}" for r in region.region])
                            output_string += (
                                f"            {BOLD}Residues: {region_str}{RESET}\n"
                            )

                            # Print per-residue observed and expected values
                            if row.oe:
                                obs_list = []
                                exp_list = []
                                for res_idx in region.region:
                                    if res_idx < len(row.oe):
                                        oe_entry = row.oe[res_idx]
                                        obs_val = (
                                            oe_entry.obs
                                            if hasattr(oe_entry, "obs")
                                            and oe_entry.obs is not None
                                            else None
                                        )
                                        exp_val = (
                                            oe_entry.exp
                                            if hasattr(oe_entry, "exp")
                                            and oe_entry.exp is not None
                                            else None
                                        )

                                        obs_list.append(
                                            f"{obs_val:7d}"
                                            if obs_val is not None
                                            else "    NA"
                                        )
                                        exp_list.append(
                                            f"{exp_val:7.2f}"
                                            if exp_val is not None
                                            else "    NA"
                                        )
                                output_string += f"            {BOLD}Observed ({len(obs_list):3d}):  {', '.join(obs_list)}{RESET}\n"
                                output_string += f"            {BOLD}Expected ({len(exp_list):3d}):  {', '.join(exp_list)}{RESET}\n"
                output_string += f"\n\n    {BOLD}Null Model (catch-all):{RESET}\n"
                output_string += (
                    f"        {BOLD}Observed: {row.null_model.obs:7d}{RESET}\n"
                )
                output_string += (
                    f"        {BOLD}Expected: {row.null_model.exp:7.2f}{RESET}\n"
                )
                # Color aggregate O/E
                if row.null_model.oe is not None:
                    oe_str = f"{row.null_model.oe:7.2f}"
                else:
                    oe_str = "     NA"
                output_string += f"        {BOLD}O/E: {oe_str}{RESET}\n"
                if (
                    hasattr(row.null_model, "region")
                    and row.null_model.region
                    and row.oe
                ):
                    obs_list = []
                    exp_list = []
                    for res_idx in row.null_model.region:
                        if res_idx < len(row.oe):
                            oe_entry = row.oe[res_idx]
                            obs_val = (
                                oe_entry.obs
                                if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                                else None
                            )
                            exp_val = (
                                oe_entry.exp
                                if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                                else None
                            )

                            obs_list.append(
                                f"{obs_val:7d}" if obs_val is not None else "    NA"
                            )
                            exp_list.append(
                                f"{exp_val:7.2f}" if exp_val is not None else "    NA"
                            )
                    output_string += f"        {BOLD}Observed ({len(obs_list):3d}): {', '.join(obs_list)}{RESET}\n"
                    output_string += f"        {BOLD}Expected ({len(exp_list):3d}): {', '.join(exp_list)}{RESET}\n"
                debug_outputs.append(output_string)

        # Print all collected debug output
        logger.info("\n".join(debug_outputs))

    # Flatten final selected + null, explode to residues, compute per-residue stats.
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
    oe_upper_func = gamma_upper_ci if oe_upper_method == "gamma" else chisq_upper_ci
    ht = ht.annotate(
        residue_index=ht.region,
        oe_upper=oe_upper_func(ht.obs, ht.exp, 0.05),
        chisq=chisq_expr,
        p_value=hl.pchisqtail(chisq_expr, 1),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index").select(
        "region_index", "obs", "exp", "oe", "oe_upper", "chisq", "p_value", "is_null"
    )

    return ht


def run_forward_no_catch_all(ht, min_exp_mis=MIN_EXP_MIS):
    num_residues = ht.oe.length()

    # keep a full set for set ops; no growing catch-all
    ht = ht.select(
        "oe",
        num_residues=num_residues,
        full_region=hl.struct(region=hl.range(num_residues)),
        regions=hl.enumerate(
            ht.min_oe_upper.map(lambda x: x.select("region")).filter(
                lambda x: x.region.length() < num_residues
            )
        ),
        selected=hl.empty_array(
            prep_region_struct(hl.range(num_residues), ht.oe).dtype
        ),
        selected_nll=0.0,
        found_best=False,
    )

    ht = ht.explode("regions")
    ht = ht.transmute(idx=ht.regions[0], region=ht.regions[1].region)
    ht = ht.annotate(
        _region=prep_region_struct(ht.region, ht.oe),
        _bg_after_selected=ht.full_region,
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "idx")

    ###ht = ht.repartition(200, shuffle=True)
    ht = ht.repartition(1, shuffle=True)

    ht = ht.checkpoint(hl.utils.new_temp_file("forward_explode", "ht"))
    round_num = 1

    while ht.aggregate(hl.agg.any(hl.is_defined(ht.region))):
        # Baseline this round: background after selected = full − union(selected)
        # Candidate background = (full − selected) − candidate
        ht = ht.annotate(
            _bg_after_both=remove_residues_from_region(
                ht._bg_after_selected, ht._region
            ),
            _current_model_nll=ht.selected_nll
            + prep_region_struct(ht._bg_after_selected.region, ht.oe).nll,
        )
        ht = ht.annotate(
            _cand_bg_model=prep_region_struct(ht._bg_after_both.region, ht.oe)
        )

        # candidate total nll vs baseline for this round
        ht = ht.annotate(
            _region=ht._region.annotate(
                region_nll=ht._region.nll,
                nll=ht.selected_nll + ht._region.nll + ht._cand_bg_model.nll,
            ),
        ).drop("_bg_after_selected", "_bg_after_both", "_cand_bg_model")
        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}.1", "ht"))

        # pick best (min nll) candidate per (uniprot, transcript)
        ht2 = ht.select("_region", "_current_model_nll")
        ht2 = ht2.filter(
            hl.is_defined(ht2._region.nll) & (ht2._region.exp >= min_exp_mis)
        )
        ht2 = (
            ht2.group_by("uniprot_id", "transcript_id")
            .aggregate(
                **hl.agg.fold(
                    hl.missing(
                        hl.tstruct(
                            min_idx=hl.tint,
                            min_nll=hl.tfloat,
                            best_region=ht2._region.dtype,
                            current_model_nll=ht2._current_model_nll.dtype,
                        )
                    ),
                    lambda accum: (
                        hl.case()
                        .when(
                            hl.is_missing(accum) | (accum.min_nll > ht2._region.nll),
                            hl.struct(
                                min_idx=ht2.idx,
                                min_nll=ht2._region.nll,
                                best_region=ht2._region,
                                current_model_nll=ht2._current_model_nll,
                            ),
                        )
                        .default(accum)
                    ),
                    lambda a, b: (
                        hl.case()
                        .when(hl.is_missing(a), b)
                        .when(hl.is_missing(b), a)
                        .when(a.min_nll <= b.min_nll, a)
                        .default(b)
                    ),
                )
            )
            .checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}.2", "ht"))
        )

        # remove chosen best_region's residues from *every* remaining candidate region
        ht2_keyed = ht2[ht.uniprot_id, ht.transcript_id]
        region_expr = remove_residues_from_region(
            hl.struct(region=ht.region), ht2_keyed.best_region
        ).region
        ht = ht.annotate(
            region=hl.or_missing(region_expr.length() > 0, region_expr),
            _best_region=ht2_keyed.best_region,
            _found_best=(
                hl.is_missing(ht2_keyed.best_region)
                | (ht2_keyed.best_region.nll >= ht2_keyed.current_model_nll)
            ),
        )

        # accept if candidate improves over current baseline
        ht = ht.annotate(
            **hl.if_else(
                ht._found_best,
                hl.struct(
                    selected=ht.selected,
                    selected_nll=ht.selected_nll,
                    region=hl.missing(ht.region.dtype),
                ),
                hl.struct(
                    selected=ht.selected.append(ht._best_region),
                    selected_nll=ht._best_region.nll,
                    region=ht.region,
                ),
            ),
        ).drop("_found_best", "_best_region")

        ht = ht.filter(hl.is_defined(ht.region) | (ht.idx == 0))
        ht = ht.annotate(_region=prep_region_struct(ht.region, ht.oe))
        ht = ht.filter(
            (hl.is_defined(ht._region.nll) & (ht._region.exp >= min_exp_mis))
            | (ht.idx == 0)
        )
        ht = ht.annotate(
            _bg_after_selected=hl.fold(
                lambda acc, r: remove_residues_from_region(acc, r),
                ht.full_region,
                ht.selected,
            )
        )

        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}", "ht"))
        round_num += 1

    # finalize
    # scoring params
    alpha = 1.0  # distance weight (Å)
    beta = 0.5  # absolute OE difference weight

    result_ht = ht.annotate(selected=add_idx_to_array(ht.selected, "region_index"))
    ht = result_ht.select(
        "oe",
        "num_residues",
        "selected",
        "selected_nll",
        residual=result_ht._bg_after_selected.region,
    )
    ht = ht.explode("residual")

    af2_dist = hl.read_table(
        "gs://gnomad/v4.1/constraint/proemis3d/test_gene_set_run/af2_dist.ht"
    )
    af2_dist = af2_dist.key_by("uniprot_id", "aa_index")
    ht = ht.annotate(
        dist_mat=af2_dist[ht.uniprot_id, ht.residual].dist_mat,
        residual=ht.oe[ht.residual],
    ).drop("oe")
    ht = ht.annotate(
        selected=ht.selected.map(
            lambda x: x.annotate(
                min_dist=hl.min(x.region.map(lambda i: ht.dist_mat[i])),
                delta=hl.abs((ht.residual.obs / ht.residual.exp) - x.oe),
            )
        ).map(lambda x: x.annotate(score=alpha * x.min_dist + beta * x.delta))
    ).drop("dist_mat")

    ht = ht.select(
        **ht.residual,
        region_index=hl.sorted(ht.selected, key=lambda x: x.score)[0].region_index,
    )
    ht = ht.group_by("uniprot_id", "transcript_id").aggregate(
        residuals=hl.agg.group_by(
            ht.region_index,
            hl.struct(
                obs=hl.agg.sum(ht.obs),
                exp=hl.agg.sum(ht.exp),
                residues=hl.agg.collect(
                    hl.struct(residue_index=ht.residue_index, assigned=True)
                ),
            ),
        )
    )
    ht = result_ht.annotate(
        residuals=ht[result_ht.uniprot_id, result_ht.transcript_id].residuals,
    )
    ht = ht.select(
        selected=ht.selected.map(
            lambda x: x.annotate(
                region=x.region.map(
                    lambda r: hl.struct(residue_index=r, assigned=False)
                ),
                add_res=ht.residuals.get(x.region_index),
            )
        )
        .map(
            lambda x: hl.if_else(
                hl.is_missing(x.add_res),
                x.annotate(source="selected"),
                x.annotate(
                    region=x.region.extend(x.add_res.residues),
                    obs=x.obs + x.add_res.obs,
                    exp=x.exp + x.add_res.exp,
                    source="selected_plus_assigned",
                ),
            ).drop("add_res")
        )
        .map(
            lambda x: prep_region_struct(
                x.region.map(lambda r: r.residue_index), ht.oe
            ).annotate(
                region=x.region,
                region_index=x.region_index,
                source=hl.or_else(x.source, "selected"),
            )
        )
    )
    ht = ht.explode("selected")
    ht = ht.select(**ht.selected)
    ht = ht.annotate(region_length=hl.len(ht.region))
    ht = ht.explode("region")
    ht = ht.annotate(**ht.region).drop("region")
    chisq_expr = calculate_oe_neq_1_chisq(ht.obs, ht.exp)
    ht = ht.annotate(
        oe_upper=gamma_upper_ci(ht.obs, ht.exp, 0.05),
        chisq=chisq_expr,
        p_value=hl.pchisqtail(chisq_expr, 1),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index").select(
        "region_index",
        "obs",
        "exp",
        "oe",
        "region_nll",
        "nll",
        "oe_upper",
        "chisq",
        "p_value",
        "source",
        "assigned",
        "region_length",
    )

    return ht


def run_forward_no_catch_all_standardized(
    ht,
    min_exp_mis: int = MIN_EXP_MIS,
    # post-hoc assignment scoring params
    alpha: float = 1.0,  # weight for standardized min distance
    beta: float = 1.0,  # weight for standardized ΔOE
):
    """
    num_residues = ht.oe.length()

    # keep a full set for set ops; no growing catch-all
    ht = ht.select(
        "oe",
        num_residues=num_residues,
        full_region=hl.struct(region=hl.range(num_residues)),
        regions=hl.enumerate(
            ht.min_oe_upper.map(lambda x: x.select("region")).filter(
                lambda x: x.region.length() < num_residues
            )
        ),
        selected=hl.empty_array(
            prep_region_struct(hl.range(num_residues), ht.oe).dtype
        ),
        selected_nll=0.0,
        found_best=False,
    )

    ht = ht.explode("regions")
    ht = ht.transmute(idx=ht.regions[0], region=ht.regions[1].region)
    ht = ht.annotate(
        _region=prep_region_struct(ht.region, ht.oe),
        _bg_after_selected=ht.full_region,
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "idx")
    ht = ht.repartition(200, shuffle=True)
    ht = ht.checkpoint(hl.utils.new_temp_file("forward_explode", "ht"))
    round_num = 1

    while ht.aggregate(hl.agg.any(hl.is_defined(ht.region))):
        # Baseline this round: background after selected = full − union(selected)
        # Candidate background = (full − selected) − candidate
        ht = ht.annotate(
            _bg_after_both=remove_residues_from_region(
                ht._bg_after_selected, ht._region
            ),
            _current_model_nll=ht.selected_nll
            + prep_region_struct(ht._bg_after_selected.region, ht.oe).nll,
        )
        ht = ht.annotate(
            _cand_bg_model=prep_region_struct(ht._bg_after_both.region, ht.oe)
        )

        # candidate total nll vs baseline for this round
        ht = ht.annotate(
            _region=ht._region.annotate(
                region_nll=ht._region.nll,
                nll=ht.selected_nll + ht._region.nll + ht._cand_bg_model.nll,
            ),
        ).drop("_bg_after_selected", "_bg_after_both", "_cand_bg_model")
        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}.1", "ht"))

        # pick best (min nll) candidate per (uniprot, transcript)
        ht2 = ht.select("_region", "_current_model_nll")
        ht2 = ht2.filter(
            hl.is_defined(ht2._region.nll) & (ht2._region.exp >= min_exp_mis)
        )
        ht2 = (
            ht2.group_by("uniprot_id", "transcript_id")
            .aggregate(
                **hl.agg.fold(
                    hl.missing(
                        hl.tstruct(
                            min_idx=hl.tint,
                            min_nll=hl.tfloat,
                            best_region=ht2._region.dtype,
                            current_model_nll=ht2._current_model_nll.dtype,
                        )
                    ),
                    lambda accum: (
                        hl.case()
                        .when(
                            hl.is_missing(accum) | (accum.min_nll > ht2._region.nll),
                            hl.struct(
                                min_idx=ht2.idx,
                                min_nll=ht2._region.nll,
                                best_region=ht2._region,
                                current_model_nll=ht2._current_model_nll,
                            ),
                        )
                        .default(accum)
                    ),
                    lambda a, b: (
                        hl.case()
                        .when(hl.is_missing(a), b)
                        .when(hl.is_missing(b), a)
                        .when(a.min_nll <= b.min_nll, a)
                        .default(b)
                    ),
                )
            )
            .checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}.2", "ht"))
        )

        # remove chosen best_region's residues from *every* remaining candidate region
        ht2_keyed = ht2[ht.uniprot_id, ht.transcript_id]
        region_expr = remove_residues_from_region(
            hl.struct(region=ht.region), ht2_keyed.best_region
        ).region
        ht = ht.annotate(
            region=hl.or_missing(region_expr.length() > 0, region_expr),
            _best_region=ht2_keyed.best_region,
            _found_best=(
                hl.is_missing(ht2_keyed.best_region)
                | (ht2_keyed.best_region.nll >= ht2_keyed.current_model_nll)
            ),
        )

        # accept if candidate improves over current baseline
        ht = ht.annotate(
            **hl.if_else(
                ht._found_best,
                hl.struct(
                    selected=ht.selected,
                    selected_nll=ht.selected_nll,
                    region=hl.missing(ht.region.dtype),
                ),
                hl.struct(
                    selected=ht.selected.append(ht._best_region),
                    selected_nll=ht._best_region.nll,
                    region=ht.region,
                ),
            ),
        ).drop("_found_best", "_best_region")

        ht = ht.filter(hl.is_defined(ht.region) | (ht.idx == 0))
        ht = ht.annotate(_region=prep_region_struct(ht.region, ht.oe))
        ht = ht.filter(
            (hl.is_defined(ht._region.nll) & (ht._region.exp >= min_exp_mis))
            | (ht.idx == 0)
        )
        ht = ht.annotate(
            _bg_after_selected=hl.fold(
                lambda acc, r: remove_residues_from_region(acc, r),
                ht.full_region,
                ht.selected,
            )
        )

        ht = ht.checkpoint(hl.utils.new_temp_file(f"forward_round_{round_num}", "ht"))
        round_num += 1

    # -------- post-hoc assignment with standardized score --------
    # Tag selected with region_index for later bookkeeping
    result_ht = ht.annotate(
        selected=add_idx_to_array(ht.selected, "region_index")
    ).checkpoint(
        "gs://gnomad/v4.1/constraint/proemis3d/test_gene_set_run/forward_no_catch_all_standardized.before_posthoc_assignment.ht",
        overwrite=True,
    )
    """
    result_ht = hl.read_table(
        "gs://gnomad/v4.1/constraint/proemis3d/test_gene_set_run/forward_no_catch_all_standardized.before_posthoc_assignment.ht"
    )

    result_ht = result_ht.annotate(
        oe=result_ht.oe.map(
            lambda x: x.annotate(
                oe=hl.if_else(x.obs == 0, 0, divide_null(x.obs, x.exp))
            )
        )
    )
    result_ht = result_ht.annotate(
        oe=hl.enumerate(result_ht.oe).map(
            lambda x: x[1].annotate(
                d_oe_med=hl.mean(  # hl.median(
                    hl.enumerate(result_ht.oe).map(
                        lambda y: hl.or_missing(x[0] != y[0], hl.abs(x[1].oe - y[1].oe))
                    )
                )
            )
        )
    )
    result_ht = result_ht.annotate(
        oe=hl.enumerate(result_ht.oe).map(
            lambda x: x[1].annotate(
                d_oe_mad=hl.mean(  # hl.median(
                    hl.enumerate(result_ht.oe).map(
                        lambda y: hl.or_missing(
                            x[0] != y[0],
                            hl.abs(hl.abs(x[1].oe - y[1].oe) - x[1].d_oe_med),
                        )
                    )
                )
            )
        )
    )

    """
    # Residual residues (not in any selected region)
    ht = result_ht.select(
        "oe",
        "num_residues",
        "selected",
        "selected_nll",
        residual=result_ht._bg_after_selected.region,
    )
    ht = ht.explode("residual")

    # distance vectors: one row per (protein, residual_residue) → dist to all residues
    af2_dist_ht = hl.read_table(
        "gs://gnomad/v4.1/constraint/proemis3d/test_gene_set_run/af2_dist.ht"
    )
    dist_mat = (
        hl.enumerate(af2_dist_ht.dist_mat)
        .filter(lambda x: x[0] != af2_dist_ht.aa_index)
        .map(lambda x: x[1])
    )

    d_med = hl.median(dist_mat)
    af2_dist_ht = (
        af2_dist_ht.annotate(
            d_med=d_med,
            d_min=hl.min(dist_mat),
            d_max=hl.max(dist_mat),
            d_mad=hl.median(dist_mat.map(lambda x: hl.abs(x - d_med))),
        )
        .key_by("uniprot_id", "aa_index")
        .checkpoint(
            "gs://gnomad-tmp-4day/persist_TableOIn8fLYabs.all_genes.ht",
            _read_if_exists=True,
            # overwrite=True,
        )
    )

    # pass 1: per-protein medians
    af2_dist_agg_ht = af2_dist_ht.group_by("uniprot_id").aggregate(
        d_min_all=hl.agg.collect(af2_dist_ht.d_min)
    )
    af2_dist_agg_ht = af2_dist_agg_ht.annotate(
        d_min_med=hl.median(af2_dist_agg_ht.d_min_all),
    )

    # pass 2: per-protein MADs (median absolute deviation)
    af2_dist_agg_ht = (
        af2_dist_agg_ht.annotate(
            d_min_mad=hl.median(
                af2_dist_agg_ht.d_min_all.map(
                    lambda x: hl.abs(x - af2_dist_agg_ht.d_min_med)
                )
            ),
        )
        .drop("d_min_all")
        .checkpoint(
            "gs://gnomad-tmp-4day/persist_TableOBsbLYS6NH.all_genes.ht",
            _read_if_exists=True,
            # overwrite=True,
        )
    )

    ht = ht.annotate(
        **af2_dist_ht[ht.uniprot_id, ht.residual],
        **af2_dist_agg_ht[ht.uniprot_id],
        residual=ht.oe[ht.residual],
    ).drop("oe")

    # candidate features per selected region for THIS residual residue
    ht = ht.annotate(
        selected=ht.selected.map(
            lambda x: x.annotate(
                min_dist=hl.min(x.region.map(lambda i: ht.dist_mat[i])),
                delta=hl.abs((ht.residual.obs / ht.residual.exp) - x.oe),
            )
        )
    ).drop("dist_mat")

    # nearest per-residual (for robust scalers)
    ht = ht.annotate(
        nearest_d=hl.min(ht.selected.map(lambda c: c.min_dist)),
        nearest_dOE=hl.min(ht.selected.map(lambda c: c.delta)),
    ).checkpoint(
        "gs://gnomad-tmp-4day/persist_TableOBsbLYS6NH.finalize2.all_genes.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    def zscore(val, med, mad):
        return (val - med) / mad

    # add standardized terms and score to each candidate
    ht = ht.annotate(
        selected=ht.selected.map(
            lambda c: c.annotate(
                d_z=zscore(c.min_dist, ht.d_med, ht.d_mad),
                deo_z=zscore(c.delta, ht.residual.d_oe_med, ht.residual.d_oe_mad),
            )
        ).map(lambda c: c.annotate(score=alpha * c.d_z + beta * hl.abs(c.deo_z)))
    ).checkpoint(
        "gs://gnomad-tmp-4day/persist_TableOBsbLYS6NH.finalize3.all_genes.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    # filter by z-thresholds, with fallback to nearest-distance only
    ht = ht.annotate(
        _choice=hl.or_missing(
            hl.len(ht.selected) > 0, hl.sorted(ht.selected, key=lambda c: c.score)[0]
        ),
        _fallback=hl.sorted(ht.selected, key=lambda c: c.min_dist)[0],
        has_selected=hl.len(ht.selected) > 0,
    )

    # final region assignment for this residual residue
    ht = ht.annotate(
        region_index=hl.if_else(
            ~ht.has_selected,
            hl.missing(hl.tint32),  # no regions exist for this protein
            hl.if_else(
                hl.is_defined(ht._choice),
                ht._choice.region_index,
                ht._fallback.region_index,
            ),
        )
    ).drop("_choice", "_fallback")

    # Group residuals by assigned region_index (drop missings)
    ht_assigned = ht.filter(hl.is_defined(ht.region_index))
    ht_assigned = (
        ht_assigned.group_by("uniprot_id", "transcript_id").aggregate(
            residuals=hl.agg.group_by(
                ht_assigned.region_index,
                hl.struct(
                    obs=hl.agg.sum(ht_assigned.residual.obs),
                    exp=hl.agg.sum(ht_assigned.residual.exp),
                    residues=hl.agg.collect(
                        hl.struct(
                            residue_index=ht_assigned.residual.residue_index,
                            assigned=True,
                        )
                    ),
                ),
            ),
        )
    ).checkpoint(
        "gs://gnomad-tmp-4day/persist_TableOBsbLYS6NH.finalize4.all_genes.ht",
        _read_if_exists=True,
        # overwrite=True,
    )
    """
    ht_assigned = hl.read_table(
        "gs://gnomad-tmp-4day/persist_TableOBsbLYS6NH.finalize4.all_genes.ht"
    )

    # stitch assignments back to selection result
    assigned_keyed = ht_assigned[result_ht.uniprot_id, result_ht.transcript_id]
    ht = result_ht.annotate(
        residuals=assigned_keyed.residuals,
    )

    # extend each selected region with its assigned residuals and recompute stats
    ht = ht.select(
        selected=ht.selected.map(
            lambda x: x.annotate(
                region=x.region.map(
                    lambda r: hl.struct(
                        residue_index=r,
                        assigned=False,
                    )
                ),
                add_res=ht.residuals.get(x.region_index),
            )
        )
        .map(
            lambda x: hl.if_else(
                hl.is_missing(x.add_res),
                x.annotate(source="selected"),
                x.annotate(
                    region=x.region.extend(x.add_res.residues),
                    obs=x.obs + x.add_res.obs,
                    exp=x.exp + x.add_res.exp,
                    source="selected_plus_assigned",
                ),
            ).drop("add_res")
        )
        .map(
            lambda x: prep_region_struct(
                x.region.map(lambda r: r.residue_index), ht.oe
            ).annotate(
                region=x.region,
                region_index=x.region_index,
                source=hl.or_else(x.source, "selected"),
            )
        )
    )
    ht = ht.explode("selected")
    ht = ht.select(**ht.selected)
    ht = ht.annotate(region_length=hl.len(ht.region))
    ht = ht.explode("region")
    ht = ht.annotate(**ht.region).drop("region")

    chisq_expr = calculate_oe_neq_1_chisq(ht.obs, ht.exp)
    ht = ht.annotate(
        oe_upper=gamma_upper_ci(ht.obs, ht.exp, 0.05),
        chisq=chisq_expr,
        p_value=hl.pchisqtail(chisq_expr, 1),
    )

    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index").select(
        "region_index",
        "obs",
        "exp",
        "oe",
        "region_nll",
        "nll",
        "oe_upper",
        "chisq",
        "p_value",
        "source",
        "assigned",
        "region_length",
    )

    return ht


def create_missense_viewer_input_ht(
    pos_ht: hl.Table,
    proemis3d_hts: Dict[str, hl.Table],
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

    proemis3d_annotated_hts = {}
    for model, proemis3d_ht in proemis3d_hts.items():
        chisq_expr = calculate_oe_neq_1_chisq(proemis3d_ht.obs, proemis3d_ht.exp)
        proemis3d_ht = proemis3d_ht.annotate(
            chisq=chisq_expr, p_value=hl.pchisqtail(chisq_expr, 1)
        )

        # Key by all fields except 'pos' and collect by key into a field named 'pos'.
        proemis3d_ht = proemis3d_ht.key_by(
            "uniprot_id", "transcript_id", "region_index", "is_null"
        ).collect_by_key("pos")

        # Sort the 'pos' field in ascending order.
        proemis3d_ht = proemis3d_ht.annotate(
            pos=hl.sorted(proemis3d_ht.pos, key=lambda x: x.residue_index)
        )

        # Annotate with 'start' and 'stop' positions for regions by merging adjacent
        # positions.
        proemis3d_ht = proemis3d_ht.annotate(
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
                        proemis3d_ht.pos[0].residue_index,
                        proemis3d_ht.pos[0].residue_index,
                        proemis3d_ht.pos[0].obs,
                        proemis3d_ht.pos[0].exp,
                        proemis3d_ht.pos[0].oe,
                        proemis3d_ht.pos[0].oe_upper,
                        proemis3d_ht.pos[0].chisq,
                        proemis3d_ht.pos[0].p_value,
                    )
                ],
                proemis3d_ht.pos[1:],
            )
        )
        proemis3d_ht = proemis3d_ht.explode("pos")

        # Key by 'gene_id' and transform 'pos' into 'start' and 'stop' fields.
        proemis3d_ht = proemis3d_ht.key_by(
            "uniprot_id", "transcript_id", "region_index", "is_null"
        )
        proemis3d_ht = proemis3d_ht.transmute(
            start=proemis3d_ht.pos[0],
            stop=proemis3d_ht.pos[1],
            obs_mis=proemis3d_ht.pos[2],
            exp_mis=proemis3d_ht.pos[3],
            obs_exp=proemis3d_ht.pos[4],
            oe_upper=proemis3d_ht.pos[5],
            chisq=proemis3d_ht.pos[6],
            p_value=proemis3d_ht.pos[7],
        )

        # Select fields in preferred order and collect by key into a field named
        # 'regions'.
        proemis3d_ht = proemis3d_ht.collect_by_key("regions")
        proemis3d_ht = proemis3d_ht.annotate(
            **pos_ht[proemis3d_ht.uniprot_id, proemis3d_ht.transcript_id]
        )
        proemis3d_ht = proemis3d_ht.annotate(
            regions=proemis3d_ht.regions.map(
                lambda x: x.select(
                    chrom=proemis3d_ht.locus_by_aapos[x.start].contig,
                    start=hl.if_else(
                        proemis3d_ht.locus_by_aapos[x.start].position
                        <= proemis3d_ht.locus_by_aapos[x.stop].position,
                        proemis3d_ht.locus_by_aapos[x.start].position,
                        proemis3d_ht.locus_by_aapos[x.stop].position,
                    ),
                    stop=hl.if_else(
                        proemis3d_ht.locus_by_aapos[x.start].position
                        <= proemis3d_ht.locus_by_aapos[x.stop].position,
                        proemis3d_ht.locus_by_aapos[x.stop].position + 2,
                        proemis3d_ht.locus_by_aapos[x.start].position + 2,
                    ),
                    aa_start=x.start,
                    aa_stop=x.stop,
                    obs_mis=x.obs_mis,
                    exp_mis=x.exp_mis,
                    obs_exp=x.obs_exp,
                    oe_upper=x.oe_upper,
                    region_index=proemis3d_ht.region_index,
                    is_null=proemis3d_ht.is_null,
                    chisq_diff_null=x.chisq,
                    p_value=x.p_value,
                )
            )
        )
        # proemis3d_ht = proemis3d_ht.group_by("transcript_id", "uniprot_id").aggregate(
        # **{f"gnomad_3d_constraint.{model}": hl.struct(
        # has_no_rmc_evidence=False,
        # passed_qc=True,
        # regions=hl.flatten(hl.agg.collect(proemis3d_ht.regions)),
        # )
        # regions=hl.flatten(hl.agg.collect(ht.regions))
        # )
        proemis3d_ht = proemis3d_ht.key_by("transcript_id")

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
    proemis3d_ht = proemis3d_ht.annotate(
        gnomad_regional_missense_constraint=hl.struct(
            has_no_rmc_evidence=False,
            passed_qc=True,
            regions=rmc_ht[proemis3d_ht.transcript_id].regions,
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
    7. Annotates flags indicating for each row if it's the "one uniprot per transcript" and/or "one transcript per gene".
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


def annotate_proemis3d_with_af2_metrics(
    proemis3D_ht: hl.Table,
    af2_plddt_ht: hl.Table,
    af2_pae_ht: hl.Table,
    af2_dist_ht: hl.Table,
) -> hl.Table:
    """
    Annotate a PROEMIS3D Hail Table with per-residue AlphaFold2 metrics (pLDDT, pAE, dist).

    This function:
    1. Explodes the pLDDT array into (residue_index, score) format.
    2. Keys the pAE and distance matrices by (uniprot_id, aa_index).
    3. Filters PROEMIS3D rows to coding transcripts (ENST).
    4. Annotates each region residue with pLDDT, pAE, and dist.
    5. Aggregates residue-level and region-level metrics.
    6. Outputs both residue-level and region-level structured annotations.

    :param proemis3D_ht: PROEMIS3D region Hail Table keyed by (uniprot_id, transcript_id,
        residue_index).
    :param af2_plddt_ht: AlphaFold2 pLDDT Hail Table with array of scores.
    :param af2_pae_ht: AlphaFold2 predicted aligned error matrix Hail Table.
    :param af2_dist_ht: AlphaFold2 distance matrix Hail Table.
    :return: Annotated PROEMIS3D Hail Table keyed by (uniprot_id, transcript_id,
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
    proemis3D_ht = proemis3D_ht.filter(proemis3D_ht.transcript_id.startswith("ENST"))

    # Annotate with per-residue metrics.
    proemis3D_ht = proemis3D_ht.annotate(
        plddt=af2_plddt_ht[proemis3D_ht.uniprot_id, proemis3D_ht.residue_index].plddt,
        pae=af2_pae_ht[proemis3D_ht.uniprot_id, proemis3D_ht.residue_index].pae,
        dist=af2_dist_ht[proemis3D_ht.uniprot_id, proemis3D_ht.residue_index].dist_mat,
    ).checkpoint(hl.utils.new_temp_file("proemis3D.annotated", "ht"))

    # Group by region.
    proemis3D_ht = proemis3D_ht.key_by("uniprot_id", "transcript_id", "region_index")
    proemis3D_ht = proemis3D_ht.collect_by_key("by_residue")

    # Compute per-residue references
    residues_expr = proemis3D_ht.by_residue.map(lambda x: x.residue_index)
    by_residue_expr = proemis3D_ht.by_residue.map(
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

    proemis3D_ht = proemis3D_ht.annotate(
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
    proemis3D_ht = proemis3D_ht.checkpoint(
        hl.utils.new_temp_file("proemis3D.region_annotated", "ht")
    )

    # Flatten by residue
    proemis3D_ht = proemis3D_ht.explode("by_residue")
    proemis3D_ht = (
        proemis3D_ht.transmute(**proemis3D_ht.by_residue)
        .key_by("uniprot_id", "transcript_id", "residue_index")
        .checkpoint(hl.utils.new_temp_file("proemis3D.by_residue.exploded", "ht"))
    )

    # Final structured output
    proemis3D_ht = proemis3D_ht.select(
        residue_level_annotations=hl.struct(
            residue_to_region_aa_dist_stats=proemis3D_ht.res.aggregate(
                lambda x: hl.agg.stats(x)
            ).annotate(median=hl.median(proemis3D_ht.res)),
            alphafold2_info=hl.struct(
                residue_plddt=proemis3D_ht.plddt,
                residue_to_region_pae_stats=proemis3D_ht.pae.aggregate(
                    lambda x: hl.agg.stats(x)
                ).annotate(median=hl.median(proemis3D_ht.pae)),
                residue_to_region_dist_stats=proemis3D_ht.dist.aggregate(
                    lambda x: hl.agg.stats(x)
                ).annotate(median=hl.median(proemis3D_ht.dist)),
            ),
        ),
        region_level_annotations=hl.struct(
            region_index=proemis3D_ht.region_index,
            region_residues=proemis3D_ht.region_residues,
            region_length=proemis3D_ht.region_length,
            obs=proemis3D_ht.obs,
            exp=proemis3D_ht.exp,
            oe=proemis3D_ht.oe,
            oe_upper=proemis3D_ht.oe_upper,
            oe_ci=oe_confidence_interval(proemis3D_ht.obs, proemis3D_ht.exp),
            chisq=proemis3D_ht.chisq,
            p_value=proemis3D_ht.p_value,
            is_null=proemis3D_ht.is_null,
            region_aa_dist_stats=proemis3D_ht.region_aa_dist_stats,
            alphafold2_info=proemis3D_ht.alphafold2_info,
        ),
    )

    return proemis3D_ht


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
    ht = ht.distinct()

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
            annotation_name=c.get("annotation_name"),
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
        dd_denovo_no_transcript_match=var_expr.dd_denovo_no_transcript_match.annotate(
            **get_kaplanis_sig_gene_annotations(ht.gene_symbol)
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
    proemis3d_ht: hl.Table,
) -> hl.Table:
    """
    Combine residue-level annotations by joining PROEMIS3D regions with COSMIS (multiple sources) and InterPro data.

    Adds the following annotations:

        - interpro
        - mtr3d
        - cosmis_alphafold
        - cosmis_pdb
        - cosmis_swiss_model
        - proemis3d

    :param ht: Input Hail Table.
    :param proemis3d_ht: PROEMIS3D Hail Table.
    :return: Annotated Hail Table.
    """
    annotation_hts = {
        **RESIDUE_LEVEL_ANNOTATION_CONFIG,
        "proemis3d": {
            "ht": proemis3d_ht,
            "keys": ["transcript_id", "uniprot_id", "residue_index"],
            "annotation_name": "proemis3d",
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
    proemis3d_ht: hl.Table,
    af2_plddt_ht: hl.Table,
    af2_pae_ht: hl.Table,
    af2_dist_ht: hl.Table,
) -> hl.Table:
    """
    Create a fully annotated per-SNV Hail Table with structured variant-, residue-, and gene-level annotations.

    :param ht: All possible SNVs Hail Table.
    :param proemis3d_ht: PROEMIS3D Hail Table.
    :param af2_plddt_ht: AlphaFold2 pLDDT Hail Table.
    :param af2_pae_ht: AlphaFold2 PAE Hail Table.
    :param af2_dist_ht: AlphaFold2 distance matrix Hail Table.
    :param partition_intervals: Partition intervals to read annotation Hail Tables with.
    :return: Annotated and checkpointed Hail Table.
    """
    hl._set_flags(use_new_shuffle="1")
    ht = annotate_snvs_with_variant_level_data(ht).naive_coalesce(5000).cache()
    hl._set_flags(use_new_shuffle=None)

    ht = ht.annotate(residue_ref=ht.aminoacid_ref)

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
        annotate_proemis3d_with_af2_metrics(
            proemis3d_ht, af2_plddt_ht, af2_pae_ht, af2_dist_ht
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
        proemis3d=ht.residue_level_annotations.proemis3d,
        gene_level_annotations=ht.gene_level_annotations,
    )

    return ht.naive_coalesce(2000)


def create_per_proemis3d_region_ht_from_residue_ht(ht: hl.Table) -> hl.Table:
    """
    Create a PROEMIS3D region-level Hail Table from a per-residue annotated Hail Table.

    This function:
    1. Extracts region_index from PROEMIS3D residue annotations.
    2. Groups rows by (transcript_id, uniprot_id, region_index).
    3. Aggregates exomes coverage, collects RMC regions, computes per-region pLDDT stats.
    4. Annotates PROEMIS3D region-level alphaFold2 info with region-level pLDDT stats.

    :param ht: Hail Table with residue-level annotations including PROEMIS3D and coverage.
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

    ht = ht.annotate(region_index=ht.proemis3d.region_level_annotations.region_index)

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
        proemis3d=hl.agg.take(ht.proemis3d.region_level_annotations, 1)[0].annotate(
            region_plddt_stats=hl.agg.stats(
                ht.proemis3d.residue_level_annotations.alphafold2_info.residue_plddt
            ).annotate(
                median=hl.median(
                    hl.agg.collect(
                        ht.proemis3d.residue_level_annotations.alphafold2_info.residue_plddt
                    )
                )
            )
        ),
        gene_level_annotations=hl.agg.take(ht.gene_level_annotations, 1)[0],
    )

    # Move pLDDT stats into alphaFold2 info.
    ht = ht.annotate(
        proemis3d=ht.proemis3d.annotate(
            alphafold2_info=ht.proemis3d.alphafold2_info.annotate(
                region_plddt_stats=ht.proemis3d.region_plddt_stats
            )
        ).drop("region_plddt_stats")
    )

    return ht.naive_coalesce(100)
