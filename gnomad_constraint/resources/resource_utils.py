# noqa: D100

from gnomad.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
    DataException,
)
from gnomad.utils.file_utils import file_exists


constraint_tmp_prefix = "gs://gnomad-tmp/constraint"
# A context table annotated with VEP, coverage, and methylation information.
annotated_context_ht = VersionedTableResource(
    default_version="85",
    versions={
        "85": TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht",
        ),
    },
)


def preprocessed_ht(data_type: str) -> TableResource:
    """
    Return TableResource of preprocessed genome, exomes, and context Table.

    The exome and genome Table will have annotations added by `prepare_ht_for_constraint_calculations` and VEP annotation from context Table.

    The context Table will have annotations added by `prepare_ht_for_constraint_calculations`.

    :param data_type: One of "exomes" or "genomes".
    :return: TableResource of processed genomes or exomes Table.
    """
    preprocessed_ht_path = f"{constraint_tmp_prefix}/model/{data_type}_processed.ht"
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
