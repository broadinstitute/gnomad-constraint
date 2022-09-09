# noqa: D100
constraint_tmp_prefix = "gs://gnomad-tmp/constraint"
context_ht_path = ""


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
    :param version: Version of annotation path to return
    :return: Output log path
    """
    return f"{constraint_tmp_prefix}/{name}.log"
