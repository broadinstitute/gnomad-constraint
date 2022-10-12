# noqa: D100
import logging

from gnomad.resources.grch37.gnomad import public_release as public_release_grch37
from gnomad.resources.grch38.gnomad import public_release as public_release_grch38
from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from gnomad.utils.file_utils import file_exists

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)

VERSIONS = ["2.1.1"]
CURRENT_VERSION = "2.1.1"
DATA_TYPES = ["context", "exomes", "genomes"]
GENOMIC_REGIONS = ["autosome_par", "chrx_nonpar", "chry_nonpar"]

constraint_tmp_prefix = "gs://gnomad-tmp/constraint"

# Context table annotated with VEP, coverage, and methylation information.
annotated_context_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.1.1": TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht",
        ),
    },
)


def get_sites_resource(data_type: str, version: str = CURRENT_VERSION):
    """
    Return genomes or exomes sites Table.

    :param data_type: One of "exomes" or "genomes".
    :param version: The version of the Table. Defaults to CURRENT_VERSION.
    :return: Genome or exomes sites Table.
    """
    if version == "2.1.1":
        return public_release_grch37(data_type).versions[version]
    else:
        return public_release_grch38(data_type).versions[version]


def get_preprocessed_ht(
    data_type: str,
    version: str = CURRENT_VERSION,
    genomic_region: str = "autosome_par",
    test: bool = False,
) -> TableResource:
    """
    Return TableResource of preprocessed genome, exomes, and context Table.

    The exome and genome Table will have annotations added by `prepare_ht_for_constraint_calculations` and VEP annotation from context Table.
    The context Table will have annotations added by `prepare_ht_for_constraint_calculations`.

    :param data_type: One of "exomes", "genomes" or "context.
    :param version: One of the release versions (`VERSIONS`). Defaults to `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par", "chrx_non_par", or "chry_non_par". Defaults to "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in chr20, chrX, and chrY. Defaults to False.
    :return: TableResource of processed genomes, exomes, or context Table.
    """
    if genomic_region not in GENOMIC_REGIONS:
        raise ValueError(f"genomic_region must be one of: {GENOMIC_REGIONS}!")
    if version not in VERSIONS:
        raise ValueError("The requested version doesn't exist!")
    preprocessed_ht_path = f"{constraint_tmp_prefix}/{version}/model/{data_type}_processed.{genomic_region}{'.test' if test else ''}.ht"
    if not file_exists(preprocessed_ht_path):
        logger.info(
            "No file or directory found at %s. Please ensure --preprocess-data is included on command line",
            preprocessed_ht_path,
        )
    return TableResource(preprocessed_ht_path)


def get_logging_path(name: str) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file
    :return: Output log path
    """
    return f"{constraint_tmp_prefix}/{name}.log"
