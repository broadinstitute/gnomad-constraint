# noqa: D100

from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
    DataException,
)
from gnomad.utils.file_utils import file_exists

VERSIONS = ["2.0", "4.0"]
CURRENT_VERSION = "2.0"
constraint_tmp_prefix = "gs://gnomad-tmp/constraint"

# A context table annotated with VEP, coverage, and methylation information.
annotated_context_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.0": TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht",
        ),
        "4.0": TableResource(
            path="",
        ),
    },
)


def get_preprocessed_ht(
    data_type: str, sex_chr: str = None, version: str = CURRENT_VERSION
) -> TableResource:
    """
    Return TableResource of preprocessed genome, exomes, and context Table.

    The exome and genome Table will have annotations added by `prepare_ht_for_constraint_calculations` and VEP annotation from context Table.

    The context Table will have annotations added by `prepare_ht_for_constraint_calculations`.

    :param data_type: One of "exomes" or "genomes".
    :param sex_chr: Which sex chromosome the dataset has. It will only be used in the file name. Defaults to None.
    :return: TableResource of processed genomes or exomes Table.
    """
    if sex_chr and sex_chr not in ("chrx", "chry"):
        raise ValueError("sex_chr must be one of: 'chrx' or 'chry'!")
    if version not in VERSIONS:
        raise ValueError("The requested version doesn't exist!")
    preprocessed_ht_path = f"{constraint_tmp_prefix}/{version}/model/{data_type}_processed{'' if sex_chr is None else f'_{sex_chr}'}.ht"
    if file_exists(preprocessed_ht_path):
        return TableResource(preprocessed_ht_path)
    else:
        raise DataException(
            f"No file or directory found at {preprocessed_ht_path}. Please add --preprocess-data to the script and rerun the pipeline."
        )


def get_logging_path(name: str) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file
    :return: Output log path
    """
    return f"{constraint_tmp_prefix}/{name}.log"
