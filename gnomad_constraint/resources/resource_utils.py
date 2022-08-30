# noqa: D100
constraint_tmp_prefix = "gs://gnomad-tmp/constraint"


def get_processed_ht_path(data_type: str) -> str:
    """
    Return path to processed genomes or exomes Table.

    :param data_type: One of "exomes" or "genomes".
    :return: Path to processed genomes or exomes Table.
    """
    return f"{constraint_tmp_prefix}/model/genomes_processed.ht"
