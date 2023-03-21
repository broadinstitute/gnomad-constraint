"""Script containing resource utility constants, reference resources, and resources of intermediate files generated by the constraint pipeline."""

import logging
from typing import Dict, Optional, Tuple

from gnomad.resources.grch37.gnomad import public_release as public_release_grch37
from gnomad.resources.grch38.gnomad import public_release as public_release_grch38
from gnomad.resources.resource_utils import (
    ExpressionResource,
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
MODEL_TYPES = ["plateau", "coverage"]
GENOMIC_REGIONS = ["full", "autosome_par", "chrx_nonpar", "chry_nonpar"]

CUSTOM_VEP_ANNOTATIONS = ["transcript_consequences", "worst_csq_by_gene"]
"""
VEP annotations used when applying models.

"transcript_consequences" option will annotate the Table with 'annotation', 'gene',
'coverage', 'transcript', and 'canonical' annotations using 'transcript_consequences'
VEP annotation.

"worst_csq_by_gene" option will annotate the Table with 'annotation', 'gene', and
'coverage' annotations using 'worst_csq_by_gene' VEP annotation.
"""

POPS = ("global", "afr", "amr", "eas", "nfe", "sas")
"""
Population labels from gnomAD.

Abbreviations stand for: global (all populations), African-American/African, Latino, East Asian, Non-Finnish European, and South Asian.
"""

COVERAGE_CUTOFF = 40
"""
Minimum median exome coverage differentiating high coverage sites from low coverage sites.

Low coverage sites require an extra calibration when computing the proportion of expected variation.
"""

# The temporary folder to store intermediate file generated by the constraint pipeline.
constraint_tmp_prefix = "gs://gnomad-tmp/constraint"

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
# Methylation Table
# TODO: decide path to methylation Table
methylation_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.1.1": TableResource(
            path="",
        ),
    },
)
# Gerp Table
# TODO: decide path to Gerp Table
gerp_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.1.1": TableResource(
            path="",
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


def get_coverage_ht(
    data_type: str,
    version: str = CURRENT_VERSION,
) -> TableResource:
    """
    Return TableResource of coverage Table.

    :param data_type: One of "exomes", "genomes" or "context.
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of coverage Table.
    """
    check_param_scope(version, data_type)
    # TODO: decide path to coverage Table
    return TableResource("")


def get_mutation_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
    use_v2_release_mutation_ht: bool = False,
) -> TableResource:
    """
    Return mutation Table that includes the baseline mutation rate for each substitution and context.

    :param version: The version of the Table. Default is CURRENT_VERSION.
    :param test: Whether the Table is for testing purposes and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :param use_v2_release_mutation_ht: Whether to use the precomputed gnomAD v2.1.1 released mutation rate table.
    :return: Mutation rate Table.
    """
    if use_v2_release_mutation_ht:
        return TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht",
        )
    else:
        check_param_scope(version)
        return TableResource(
            f"{constraint_tmp_prefix if test else get_constraint_release_prefix(version)}/gnomad.v{version}.constraint.mutation_rate.{'.test' if test else ''}.ht"
        )


def get_annotated_context_ht(
    version: str = CURRENT_VERSION,
    use_old_data: bool = False,
) -> TableResource:
    """
    Return TableResource of annotated context Table.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param use_old_data: Whether to use old annotated context Table. Default is False.
    :return: TableResource of annotated context Table.
    """
    if use_old_data:
        return TableResource(
            "gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht"
        )

    check_param_scope(version)
    return TableResource(
        f"{constraint_tmp_prefix}/{version}/preprocessed_data/annotated_context.ht"
    )


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
        f"{constraint_tmp_prefix}/{version}/preprocessed_data/{data_type}_processed.{genomic_region}{'.test' if test else ''}.ht"
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
        f"{constraint_tmp_prefix}/{version}/training_data/constraint_training.{genomic_region}{'.test' if test else ''}.ht"
    )


def get_models(
    model_type: str,
    version: str = CURRENT_VERSION,
    genomic_region: str = "autosome_par",
    test: bool = False,
) -> ExpressionResource:
    """
    Return path to a HailExpression that contains desired model type.

    :param model_type: The type of model. One of "plateau", "coverage". Default is None.
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: Path to the specified model.
    """
    check_param_scope(
        version=version, genomic_region=genomic_region, model_type=model_type
    )
    return ExpressionResource(
        f"{constraint_tmp_prefix}/{version}/models/{model_type}.{genomic_region}{'.test' if test else ''}.he"
    )


def get_predicted_proportion_observed_dataset(
    custom_vep_annotation: str = "transcript_consequences",
    version: str = CURRENT_VERSION,
    genomic_region: str = "autosome_par",
    test: bool = False,
) -> str:
    """
    Return TableResource containing the expected variant counts and observed:expected ratio.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint model (one of "transcript_consequences" or "worst_csq_by_gene").
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: Path of the model.
    """
    check_param_scope(
        version=version,
        genomic_region=genomic_region,
        custom_vep_annotation=custom_vep_annotation,
    )
    return TableResource(
        f"{constraint_tmp_prefix}/{version}/predicted_proportion_observed/{custom_vep_annotation}/predicted_proportion_observed.{genomic_region}{'.test' if test else ''}.ht"
    )


def get_constraint_metrics_dataset(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> str:
    """
    Return TableResource of testing dataset with expected variant count and observed:expected ratio.

    :param custom_model: The customized model (one of "standard" or "worst_csq").
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: Path of the model.
    """
    check_param_scope(
        version=version,
    )
    return TableResource(
        f"gs://gnomad/{version}/constraint/metrics/constraint_metrics.{'.test' if test else ''}.ht"
    )


def check_resource_existence(
    input_step_resources: Optional[Dict[str, Tuple]] = None,
    output_step_resources: Optional[Dict[str, Tuple]] = None,
    overwrite: bool = False,
) -> None:
    """
    Check whether all the input files exist and the overwrite parameter is set to True when writing the output files.

    If no parameters are passed to the function, nothing is done.

    :param input_step_resources: A dictionary with keys as pipeline steps that generate
        input files and the value as the input files to check the existence of. Default
        is None.
    :param output_step_resources: A dictionary with keys as pipeline step that generates
        output files and the value as the output files to check the existence of.
        Default is None.
    :param overwrite: The overwrite parameter used when writing the output files.
        Default is None.
    :return: None.
    """
    # Check if the input resources exist.
    if input_step_resources:
        for step, input_resources in input_step_resources.items():
            check_file_exists_raise_error(
                input_resources,
                error_if_not_exists=True,
                error_if_not_exists_msg=(
                    f"Not all input resources exist. Please add {step} to "
                    "the command line. The following files are missing: "
                ),
            )

    # Check if the output resources exist when `overwrite` is False.
    if not overwrite and output_step_resources:
        for step, output_resources in output_step_resources.items():
            check_file_exists_raise_error(
                output_resources,
                error_if_exists=True,
                error_if_exists_msg=(
                    "Some of the output resources that will be created by "
                    f"{step} already exist and the --overwrite argument "
                    f"was not set. Please rerun {step} with --overwrite. "
                    "The following files already exist: "
                ),
            )


def check_param_scope(
    version: Optional[str] = None,
    genomic_region: Optional[str] = None,
    data_type: Optional[str] = None,
    model_type: Optional[str] = None,
    custom_vep_annotation: Optional[str] = None,
) -> None:
    """
    Check if the specified version, genomic region, and data type are in the scope of the constraint pipeline.

    :param version: One of the release versions (`VERSIONS`). Default is None.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is None.
    :param data_type: One of "exomes", "genomes" or "context". Default is None.
    :param model_type: One of "plateau", "coverage". Default is None.
    :param custom_vep_annotation: The VEP annotation used to customize the constraint model (one of "transcript_consequences" or "worst_csq_by_gene").
        Default is None.
    """
    if version and version not in VERSIONS:
        raise ValueError("The requested version doesn't exist!")
    if genomic_region and genomic_region not in GENOMIC_REGIONS:
        raise ValueError(f"genomic_region must be one of: {GENOMIC_REGIONS}!")
    if data_type and data_type not in DATA_TYPES:
        raise ValueError(f"data_type must be one of: {DATA_TYPES}!")
    if model_type and model_type not in MODEL_TYPES:
        raise ValueError(f"model_type must be one of: {MODEL_TYPES}!")
    if custom_vep_annotation and custom_vep_annotation not in CUSTOM_VEP_ANNOTATIONS:
        raise ValueError(
            f"custom_vep_annotation must be one of: {CUSTOM_VEP_ANNOTATIONS}!"
        )


def get_logging_path(name: str) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :return: Output log path.
    """
    return f"{constraint_tmp_prefix}/{name}.log"
