# noqa: D100
constraint_tmp_prefix = "gs://gnomad-tmp/constraint"
# This is an annotated v2 context table with coverage and methylation information.
context_ht_path = "gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht"


def get_processed_ht_path(data_type: str) -> str:
    """
    Return path to processed genomes or exomes Table.

    :param data_type: One of "exomes" or "genomes".
    :return: Path to processed genomes or exomes Table.
    """
    return f"{constraint_tmp_prefix}/model/{data_type}_processed.ht"


def get_logging_path(name: str) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file
    :return: Output log path
    """
    return f"{constraint_tmp_prefix}/{name}.log"
