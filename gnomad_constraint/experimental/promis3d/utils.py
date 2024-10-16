"""Script with utility functions for the Promis3D pipeline."""

import os
from io import StringIO
from typing import List, Optional, Union

import hail as hl
import numpy as np
import pandas as pd
from Bio.PDB import PPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from hail.utils.misc import divide_null
from pyspark.sql.functions import col, explode, pandas_udf, rtrim, split
from pyspark.sql.types import StringType, StructField, StructType

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

COLNAMES_TRANSLATIONS = [
    "enst",
    "ensg",
    "havana_g",
    "havana_t",
    "transcript",
    "gene",
    "aalength",
]
"""
Column names for the GENCODE translations Hail Table.
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
def get_structure_peptide(structure) -> str:
    """
    Get the sequence from a structure.

    :param structure:
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
) -> Union[str, np.ndarray]:
    """
    Process AlphaFold2 MMCIF content.

    :param uniprot_id: UniProt ID.
    :param mmcif_content: MMCIF content.
    :param distance_matrix: Whether to return the distance matrix. Default is False,
        which returns the peptide sequence.
    :return: Peptide sequence as a string or distance matrix as a NumPy array.
    """
    # Create an MMCIFParser object with quiet mode to suppress warnings.
    parser = MMCIFParser(QUIET=True)

    # Create a StringIO object from the file content.
    mmcif_io = StringIO(mmcif_content)

    # Parse the MMCIF content.
    structure = parser.get_structure(uniprot_id, mmcif_io)

    if distance_matrix:
        return get_structure_dist_matrix(structure)
    else:
        return get_structure_peptide(structure)


def process_af2_structures(
    gcs_bucket_glob: str, distance_matrix: bool = False
) -> hl.Table:
    """
    Process AlphaFold2 structures from a GCS bucket.

    .. note::

        All files in the bucket must be in CIF format with a '.cif.gz' extension.

    :param gcs_bucket_glob: GCS bucket glob pattern.
    :param distance_matrix: Whether to return the distance matrix. Default is False,
        which returns the peptide sequence.
    :return: Hail Table with UniProt IDs and sequences or distance matrices.
    """
    # Get Spark session for file distribution and processing.
    spark = hl.utils.java.Env.spark_session()

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
    if distance_matrix:
        col_name = "dist_mat"
        rtype = "array<array<float>>"
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
            af2_data = process_af2_mmcif(uniprot_id, file_content, distance_matrix)
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
    tmp_path = hl.utils.new_temp_file("process_af2_structures", "parquet")
    result_df.write.mode("overwrite").save(tmp_path)
    result_df = spark.read.parquet(tmp_path)

    # Convert the Spark DataFrame to a Hail Table.
    ht = hl.Table.from_spark(result_df, key=["af2_file", "uniprot_id"])

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
) -> hl.Table:
    """
    Get GENCODE positions for the given transcripts and translations.

    :param transcripts_ht: Hail Table with GENCODE transcripts.
    :param translations_ht: Hail Table with GENCODE translations.
    :param gencode_gtf_ht: Hail Table with GENCODE GTF data.
    :return: Hail Table with GENCODE positions.
    """
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
    cds_len_ht = cds_len_ht.filter(
        (cds_len_ht.gtf_cds_len == cds_len_ht.cds_len) & (cds_len_ht.cds_len % 3 == 0)
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
            lambda gp, r, ap: hl.struct(locus=hl.locus(ht.chrom, gp), ref=r, aapos=ap),
            ht.gtf_cds_pos,
            hl.range(hl.len(ht.ref)).map(lambda i: ht.ref[i]),
            ht.aapos,
        )
    ).explode("positions")

    return ht.transmute(**ht.positions)


########################################################################################
# Functions to perform tasks from run_greedy.R
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

    return oe_codon_ht


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


def get_3d_residue(
    dist_mat_expr: hl.expr.ArrayExpression,
    oe_expr: hl.expr.ArrayExpression,
    alpha: float = 0.05,
) -> hl.expr.StructExpression:
    """
    Get the 3D residue with the lowest upper bound of the OE confidence interval.

    :param dist_mat_expr: Array expression with distance matrix.
    :param oe_expr: Array expression with observed and expected values.
    :param alpha: Significance level for the OE confidence interval. Default is 0.05.
    :return: Struct expression with the 3D residue with the lowest upper bound of the OE
        confidence interval.
    """
    dist_mat_expr = add_idx_to_array(dist_mat_expr, "aa_index", element_name="dist")
    dist_mat_expr = hl.sorted(dist_mat_expr, key=lambda x: x.dist)
    dist_mat_expr = add_idx_to_array(dist_mat_expr, "dist_index")

    oe_expr = dist_mat_expr.map(
        lambda x: x.annotate(**oe_expr[x.aa_index]).drop("dist")
    )
    oe_expr = hl.array_scan(
        lambda i, j: j.annotate(obs=i.obs + j.obs, exp=i.exp + j.exp),
        oe_expr[0],
        oe_expr[1:],
    )

    # Calculate upper bound of oe confidence interval.
    min_loeuf_expr = hl.sorted(
        oe_expr.map(
            lambda x: x.annotate(
                oe=divide_null(x.obs, x.exp),
                oe_upper=(
                    hl.qchisqtail(1 - alpha / 2, 2 * (x.obs + 1), lower_tail=True)
                    / (2 * x.exp)
                ),
            )
        ),
        key=lambda x: x.oe_upper,
    )[0]
    oe_expr = oe_expr[: min_loeuf_expr.dist_index]
    min_loeuf_expr = min_loeuf_expr.drop("dist_index")

    return min_loeuf_expr.annotate(region=oe_expr.map(lambda x: x.aa_index))


def run_greedy(af2_ht: hl.Table, oe_codon_ht: hl.Table) -> hl.Table:
    """
    Run the greedy algorithm to find the most intolerant region.

    :param af2_ht: Hail Table with AlphaFold2 data.
    :param oe_codon_ht: Hail Table with observed and expected values for codons.
    :return: Hail Table with the most intolerant region for each UniProt ID and residue
        index
    """
    af2_ht = af2_ht.annotate(oe=oe_codon_ht[af2_ht.uniprot_id].oe)

    min_loeuf_expr = af2_ht.dist_mat.map(lambda x: get_3d_residue(x, af2_ht.oe))
    min_loeuf_expr = hl.sorted(min_loeuf_expr, key=lambda x: x.oe)
    min_loeuf_expr = add_idx_to_array(min_loeuf_expr, "region_index")

    initial_score_expr = min_loeuf_expr.map(
        lambda x: hl.missing(min_loeuf_expr.dtype.element_type)
    )
    score_expr = hl.fold(
        lambda i, j: hl.enumerate(i).map(
            lambda x: hl.coalesce(x[1], hl.or_missing(j.region.contains(x[0]), j))
        ),
        initial_score_expr,
        min_loeuf_expr,
    )
    score_expr = add_idx_to_array(score_expr, "residue_index")
    ann_keep = ["residue_index", "region_index", "obs", "exp", "oe", "oe_upper"]
    af2_ht = af2_ht.select(score=score_expr.map(lambda x: x.select(*ann_keep)))
    af2_ht = af2_ht.explode("score")

    return af2_ht.select(**af2_ht.score).key_by("uniprot_id", "residue_index")
