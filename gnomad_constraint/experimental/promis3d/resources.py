"""Resource definitions for the promis3d pipeline."""

import logging

import hail as hl
from gnomad.resources.grch37.reference_data import gencode as grch37_gencode
from gnomad.resources.grch38.reference_data import gencode as grch38_gencode
from gnomad.resources.resource_utils import (
    BaseResource,
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("promis3d_pipeline")
logger.setLevel(logging.INFO)

VERSIONS = ["2.1.1", "4.1"]
"""Possible gnomAD versions for the promis3d pipeline."""

CURRENT_VERSION = "2.1.1"
"""Current gnomAD version for the promis3d pipeline."""

GENCODE_VERSION_MAP = {
    "2.1.1": "19",
    "4.1": "39",
}
"""GENCODE version map for each gnomAD version."""

def get_promis3d_root(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Get root path to promis3d resources.

    :param version: Version of promis3d resources to use.
    :param test: Whether to use a tmp path for testing.
    :return: Root path to promis3d resources.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/constraint/promis3d"
        if test
        else f"gs://gnomad/v{version}/constraint/promis3d"
    )


def get_gencode_fasta(
    version: str = CURRENT_VERSION,
    name: str = "pc_transcripts",
) -> str:
    """
    Get GENCODE FASTA file path.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param name: Name of the type of GENCODE FASTA file to get. One of 'pc_transcripts',
        'pc_translations'. Default is 'pc_transcripts'.
    :return: GENCODE FASTA file path.
    """
    # TODO: Change this path when moved to a more permanent location.
    return (
        f'gs://gnomad-julia/promis3d/resources/'
        f'gencode.v{GENCODE_VERSION_MAP[version]}.{name}.fa.gz'
    )


def get_alpha_fold2_dir(version: str = CURRENT_VERSION) -> str:
    """
    Get AlphaFold2 directory path.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :return: AlphaFold2 directory path.
    """
    # TODO: Change this path when moved to a more permanent location and add any
    #  needed versioning.
    return 'gs://gnomad-julia/alphafold2'


def get_gencode_seq_ht(
    version: str = CURRENT_VERSION,
    name: str = "pc_transcripts",
    test: bool = False,
) -> TableResource:
    """
    Get GENCODE sequences Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param name: Name of the type of GENCODE sequences to get. One of 'pc_transcripts',
        'pc_translations'. Default is 'pc_transcripts'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE sequences Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/"
        f"gencode_sequences.{name}.ht"
    )


def get_af2_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get alphafold2 Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE sequences Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/af2.ht"
    )


def get_gencode_translations_matched_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get GENCODE translations matched Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE translations matched Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/"
        f"gencode_translations_matched.ht"
    )


def get_gencode_cds_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get GENCODE CDS Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :return: GENCODE CDS Hail Table resource.
    """
    if version == "2.1.1":
        gencode = grch37_gencode
    elif version == "4.1":
        gencode = grch38_gencode
    else:
        raise ValueError(f"Invalid version: {version}")

    ht = gencode.versions[f"v{GENCODE_VERSION_MAP[version]}"].ht()
    ht = ht.filter(
        (ht.feature == "CDS") & (ht.transcript_type == "protein_coding")
    )

    return ht
