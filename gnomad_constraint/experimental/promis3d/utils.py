"""Script with utility functions for the Promis3D pipeline."""

import logging
import os
from io import StringIO
from typing import Iterator, List, Optional, Union

import hail as hl
import numpy as np
import pandas as pd
from Bio.PDB import PPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from hail.utils.misc import divide_null
from pyspark.sql.functions import col, explode, pandas_udf, rtrim, split
from pyspark.sql.types import StringType, StructField, StructType

from gnomad_constraint.experimental.promis3d.constants import MIN_EXP_MIS

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
    if distance_matrix:
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
    if distance_matrix:
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
    oe_expr = hl.array_scan(
        lambda i, j: j.annotate(obs=i.obs + j.obs, exp=i.exp + j.exp),
        oe_expr[0],
        oe_expr[1:],
    )

    return oe_expr


def calculate_oe_upper(oe_expr, alpha=0.05):
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
    return region_expr.map(lambda x: oe_expr[x])


def get_agg_oe_for_region(region_expr):
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
    # Calculate negative log-likelihood a region.
    return hl.sum(
        region_expr.map(lambda x: -hl.dpois(x.obs, gamma_expr * x.exp, log_p=True))
    )


def getAIC(region_expr, nll):
    if isinstance(region_expr, hl.expr.ArrayExpression):
        region_count = region_expr.length()
    else:
        region_count = hl.int(region_expr.region.length() > 0)

    return 2 * region_count + 2 * nll


def remove_residues_from_region(region_expr, remove_region_expr):
    remove_region_residues = hl.set(remove_region_expr.region)
    updated_region_expr = hl.set(region_expr.region).difference(remove_region_residues)
    return hl.or_missing(
        hl.is_defined(region_expr.region) & hl.is_defined(remove_region_expr.region),
        region_expr.annotate(region=hl.array(updated_region_expr)),
    )


def get_min_region(regions_expr, min_field="oe_upper", min_exp_mis=None):
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
    oe_expr = annotate_region_with_oe(region_expr, oe_expr)
    oe_agg_expr = get_agg_oe_for_region(oe_expr)
    nll_expr = calculate_neg_log_likelihood(oe_expr, oe_agg_expr.oe)
    return hl.struct(
        region=region_expr,
        **oe_agg_expr,
        region_nll=nll_expr,
        nll=nll_expr,
    )


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

    ht = ht.select(
        selected=add_idx_to_array(ht.selected.append(ht.null_model), "region_index")
    )
    ht = ht.explode("selected")
    ht = ht.select(**ht.selected).explode("region")
    ht = ht.annotate(
        residue_index=ht.region,
        oe_upper=(
            hl.qchisqtail(1 - 0.05 / 2, 2 * (ht.obs + 1), lower_tail=True)
            / (2 * ht.exp)
        ),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index").select(
        "region_index", "obs", "exp", "oe", "oe_upper"
    )

    return ht
