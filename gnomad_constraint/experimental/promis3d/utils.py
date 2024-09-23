from typing import List
import gzip
import os
import hail as hl
from pyspark.sql import SparkSession
from pyspark.sql.functions import pandas_udf, col, explode, split, rtrim
from pyspark.sql.types import StringType, StructType, StructField
import pandas as pd
import numpy as np
from pyspark.sql.types import StringType
from io import StringIO
from Bio import SeqIO, SeqFeature
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PPBuilder


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
    df = df.select(explode(split(df['value'], '>')).alias("sequence"))

    # Convert the Spark DataFrame to a Hail Table.
    ht = hl.Table.from_spark(df)
    ht = ht.filter(ht.sequence != "")

    split_expr = ht.sequence.split('\n')
    split_info_expr = split_expr[0].split('\|')
    ht = ht.select(
        **{
            n: hl.or_missing(
                (split_info_expr.length() > i) & ~hl.array(['', '-']).contains(
                    split_info_expr[i]),
                split_info_expr[i]
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
    names = ['utr5', 'cds', 'utr3']
    ht = ht.rename({f"index{i + 1}": n for i, n in enumerate(names)})

    # Split the 'utr5', 'cds', and 'utr3' annotations into start and end positions.
    ht = ht.annotate(
        **{
            n: hl.bind(lambda x: x.map(hl.int), ht[n].split(":")[1].split("-"))
            for n in names
        }
    )

    # Trim the sequence to the CDS range.
    ht = ht.annotate(cds_sequence=ht.sequence[ht.cds[0] - 1:ht.cds[1]])

    return ht


########################################################################################
# Note that the functionality in split_context_obs_exp.R is not implemented here
# because it just splits the tsv file by transcript, we will be directly using the
# Hail Table instead.
########################################################################################


########################################################################################
# Functions to perform tasks from read_af2_sequences.R
########################################################################################
def read_structure_peptide(uniprot_id: str, mmcif_content: str) -> str:
    """
    Read the sequence from a structure file in MMCIF format.

    :param uniprot_id: UniProt ID.
    :param mmcif_content: MMCIF content.
    :return: Sequence as a string.
    """
    # Create an MMCIFParser object with quiet mode to suppress warnings.
    parser = MMCIFParser(QUIET=True)

    # Create a StringIO object from the file content.
    mmcif_io = StringIO(mmcif_content)

    # Parse the MMCIF content.
    structure = parser.get_structure(uniprot_id, mmcif_io)

    ppb = PPBuilder()

    # Return the sequence as a string.
    return ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(structure)])


@pandas_udf("uniprot_id string, sequence string")
def process_file(
    file_path_series: pd.Series,
    file_content_series: pd.Series
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
        uniprot_id = os.path.basename(file_path).split('-')[1]

        # Process the file content.
        sequence = read_structure_peptide(uniprot_id, file_content)
        result.append((uniprot_id, sequence))

    return pd.DataFrame(result, columns=["uniprot_id", "sequence"])


def process_af2_structures(gcs_bucket_path: str) -> hl.Table:
    """
    Process AlphaFold2 structures from a GCS bucket.

    .. note::

        All files in the bucket must be in CIF format with a '.cif.gz' extension.

    :param gcs_bucket_path: GCS bucket path with AlphaFold2 structures.
    :return: Hail Table with UniProt IDs and sequences.
    """
    # Get Spark session for file distribution and processing.
    spark = hl.utils.java.Env.spark_session()

    # Define schema for loading the files.
    schema = StructType([
        StructField("file_content", StringType(), True),
        StructField("af2_file", StringType(), True)
    ])

    # Use Spark to read files in parallel.
    # This reads the entire content of each file as a (filename, content) pair.
    af2_files_df = spark.read.format("text").load(f"{gcs_bucket_path}/*.cif.gz", schema=schema, wholetext=True).withColumn("af2_file", col("_metadata.file_path"))

    # Apply the Pandas UDF to process the files.
    result_df = af2_files_df.withColumn("uniprot_id_sequence", process_file(col("af2_file"), col("file_content")))

    # Split the resulting column into separate columns.
    result_df = result_df.select("af2_file", col("uniprot_id_sequence.uniprot_id"), col("uniprot_id_sequence.sequence"))

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
    uniprots_with_multifrags = ht.filter(ht.af2_file.contains("-F2-")).select().distinct()

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
    ht1 = ht1.key_by('sequence')
    ht2 = ht2.key_by('sequence')

    return ht1.join(ht2, how='inner')
