"""Script containing resource utility constants, reference resources, and resources of intermediate files generated by the constraint pipeline."""

import logging
from typing import Optional, Tuple

from gnomad.resources.grch37.gnomad import public_release as public_release_grch37
from gnomad.resources.grch38.gnomad import public_release as public_release_grch38
from gnomad.resources.resource_utils import (
    GnomadPublicTableResource,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import check_file_exists_raise_error

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
POPS = ("global", "afr", "amr", "eas", "nfe", "sas")
"""
Stands for global, African-American/African, Latino, East Asian, Non-Finnish European,
and South Asian population. Labels are from gnomAD.
"""

# The temporary folder to store intermediate file generated by the constraint pipeline.
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
# Mutation rate Table that include the baseline mutation rate for each substitution
# and context.
mutation_rate_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.1.1": TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht",
        ),
    },
)


def get_sites_resource(
    data_type: str, version: str = CURRENT_VERSION
) -> GnomadPublicTableResource:
    """
    Return genomes or exomes sites Table.

    :param data_type: One of "exomes" or "genomes".
    :param version: The version of the Table. Default is CURRENT_VERSION.
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

    The exome and genome Table will have annotations added by
    `prepare_ht_for_constraint_calculations()` and VEP annotation from context Table.

    The context Table will have annotations added by
    `prepare_ht_for_constraint_calculations()`.

    :param data_type: One of "exomes", "genomes" or "context.
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is "autosome_par".
    :param test: Whether the Table is for testing purposes and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: TableResource of processed genomes, exomes, or context Table.
    """
    check_param_scope(version, genomic_region, data_type)
    return TableResource(
        f"{constraint_tmp_prefix}/{version}/model/{data_type}_processed.{genomic_region}{'.test' if test else ''}.ht"
    )


def get_training_dataset(
    version: str = CURRENT_VERSION,
    genomic_region: str = "autosome_par",
    test: bool = False,
) -> TableResource:
    """
    Return TableResource of training dataset with observed and possible variant count.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: TableResource of training dataset.
    """
    check_param_scope(version, genomic_region)
    return TableResource(
        f"{constraint_tmp_prefix}/{version}/model/training/constraint_training.{genomic_region}{'.test' if test else ''}.ht"
    )


def check_resource_existence(
    input_pipeline_step: Optional[str] = None,
    output_pipeline_step: Optional[str] = None,
    input_resources: Optional[Tuple] = None,
    output_resources: Optional[Tuple] = None,
    overwrite: bool = False,
) -> None:
    """
    Check whether all the input files exist and the overwrite parameter is set to True when writing the output files.

    If no parameters are passed to the function, nothing is done.

    :param input_pipeline_step: The pipeline step that generates input files. Default
        is None.
    :param output_pipeline_step: The pipeline step that generates output files. Default
        is None.
    :param input_resources: Paths of the input files to check the existence of. Default
        is None.
    :param output_resources: Paths of the output files to check the existence of.
        Default is None.
    :param overwrite: The overwrite parameter used when writing the output files.
        Default is None.
    :return: None.
    """
    # Check if the input resources exist.
    if input_pipeline_step and input_resources:
        check_file_exists_raise_error(
            [r.path for r in input_resources],
            error_if_not_exists=True,
            error_if_not_exists_msg=(
                f"Not all input resources exist. Please add {input_pipeline_step} to "
                "the command line. The following files are missing: "
            ),
        )

    # Check if the output resources exist when `overwrite` is False.
    if not overwrite and output_pipeline_step and output_resources:
        check_file_exists_raise_error(
            [r.path for r in output_resources],
            error_if_exists=True,
            error_if_exists_msg=(
                "Some of the output resources that will be created by "
                f"{output_pipeline_step} already exist and the --overwrite argument "
                f"was not set. Please rerun {output_pipeline_step} with --overwrite. "
                "The following files already exist: "
            ),
        )


def check_param_scope(
    version: Optional[str] = None,
    genomic_region: Optional[str] = None,
    data_type: Optional[str] = None,
) -> None:
    """
    Check if the specified version, genomic region, and data type are in the scope of the constraint pipeline.

    :param version: One of the release versions (`VERSIONS`). Default is None.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is None.
    :param data_type: One of "exomes", "genomes" or "context". Default is None.
    """
    if version and version not in VERSIONS:
        raise ValueError("The requested version doesn't exist!")
    if genomic_region and genomic_region not in GENOMIC_REGIONS:
        raise ValueError(f"genomic_region must be one of: {GENOMIC_REGIONS}!")
    if data_type and data_type not in DATA_TYPES:
        raise ValueError(f"data_type must be one of: {DATA_TYPES}!")


def get_logging_path(name: str) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :return: Output log path.
    """
    return f"{constraint_tmp_prefix}/{name}.log"
