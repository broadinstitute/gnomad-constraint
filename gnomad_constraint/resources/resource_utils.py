"""Script containing resource utility constants, reference resources, and resources of intermediate files generated by the constraint pipeline."""

import logging
from typing import Dict, Optional, Tuple, Union

import gnomad.resources.grch37.gnomad as gnomad_grch37
import gnomad.resources.grch37.reference_data as ref_grch37
import gnomad.resources.grch38.gnomad as gnomad_grch38
import gnomad.resources.grch38.reference_data as ref_grch38
import hail as hl
from gnomad.resources.resource_utils import (
    BaseResource,
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)
from gnomad_qc.v4.resources.release import release_coverage, release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)

VERSIONS = ["2.1.1", "4.0"]
CURRENT_VERSION = "2.1.1"
DATA_TYPES = ["context", "exomes", "genomes"]
MODEL_TYPES = ["plateau", "coverage"]
GENOMIC_REGIONS = ["autosome_par", "chrx_nonpar", "chry_nonpar"]

CUSTOM_VEP_ANNOTATIONS = ["transcript_consequences", "worst_csq_by_gene"]
"""
VEP annotations used when applying models.

"transcript_consequences" option will annotate the Table with 'annotation', 'gene',
'coverage', 'transcript', and either 'canonical' or 'mane_select' annotations using 'transcript_consequences'
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

# VEP context Table.
vep_context_ht = VersionedTableResource(
    CURRENT_VERSION,
    versions={
        "2.1.1": ref_grch37.vep_context.versions["85"],
        "4.0": ref_grch38.vep_context.versions["105"],
    },
)


def get_constraint_root(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Return path to constraint root folder.

    :param version: Version of constraint path to return.
    :param test: Whether to use a tmp path.
    :return: Root path to constraint resources.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/constraint"
        if test
        else f"gs://gnomad/v{version}/constraint"
    )


def get_sites_resource(data_type: str, version: str = CURRENT_VERSION) -> BaseResource:
    """
    Return genomes or exomes sites Table.

    :param data_type: One of "exomes" or "genomes".
    :param version: The version of the Table. Default is CURRENT_VERSION.
    :return: Genome or exomes sites Table.
    """
    build = check_param_scope(version=version, data_type=data_type)
    if build == "GRCh37":
        return gnomad_grch37.public_release(data_type).versions[version]
    elif version == "4.0":
        # TODO: This can change when gnomAD v4.0 is publicly released.
        if data_type == "genomes":
            return gnomad_grch38.public_release(data_type).versions["3.1.2"]
        else:
            return release_sites().versions[version]
    else:
        raise ValueError(
            "The sites resource has not been defined for the specified version!"
        )


def get_gerp_ht(build: str) -> hl.Table:
    """
    Retrieve publicly released GERP scores.

    :param build: Build of the reference genomes to use. One of "GRCh37" or "GRCh38".
    :return: Table of GERP scores stored in 'S' annotation for the specified build.
    """
    return hl.experimental.load_dataset(
        name="gerp_scores", version="hg19", reference_genome=build
    )


def get_methylation_ht(build: str) -> TableResource:
    """
    Retrieve methylation data TableResource.

    :param build: Build of the reference genomes to use. One of "GRCh37" or "GRCh38".
    :return: TableResource of methylation data for the specified build.
    """
    if build == "GRCh37":
        return ref_grch37.methylation_sites
    elif build == "GRCh38":
        methylation_chrx = methylation_sites_chrx.ht()
        methylation_autosomes = methylation_sites.ht()
        methylation_ht = methylation_autosomes.union(methylation_chrx)
        methylation_ht = methylation_ht.checkpoint("gs://gnomad-tmp/methylation_ht.ht")
        return TableResource(
            path="gs://gnomad-tmp/methylation_ht.ht",
        )
    else:
        raise ValueError("Build must be one of 'GRCh37' or 'GRCh38'.")


def get_coverage_ht(
    data_type: str,
    version: str = CURRENT_VERSION,
) -> VersionedTableResource:
    """
    Return TableResource of coverage Table.

    :param data_type: One of "exomes", "genomes".
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of coverage Table.
    """
    build = check_param_scope(version=version, data_type=data_type)
    if build == "GRCh37":
        return gnomad_grch37.coverage(data_type)
    else:
        # TODO: This can change when gnomAD v4.0 is publicly released.
        if data_type == "genomes":
            return gnomad_grch38.coverage(data_type)
        elif version == "4.0":
            return release_coverage().versions[version]
        else:
            raise ValueError(
                "The coverage Table resource has not been defined for exomes of "
                f"version: {version}!"
            )


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
    :param use_v2_release_mutation_ht: Whether to use the precomputed gnomAD v2.1.1
        released mutation rate table.
    :return: Mutation rate Table.
    """
    if use_v2_release_mutation_ht:
        return TableResource(
            path="gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht",
        )
    else:
        check_param_scope(version)
        return TableResource(
            f"{get_constraint_root(version, test)}/mutation_rate/gnomad.v{version}.mutation_rate.ht"
        )


def get_annotated_context_ht(
    version: str = CURRENT_VERSION,
    use_v2_context_ht: bool = False,
) -> TableResource:
    """
    Return TableResource of annotated context Table.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param use_v2_context_ht: Whether to use annotated context Table that was produced
        for gnomAD v2. Default is False.
    :return: TableResource of annotated context Table.
    """
    if use_v2_context_ht:
        return TableResource(
            "gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht"
        )

    check_param_scope(version)
    return TableResource(
        f"{get_constraint_root(version)}/preprocessed_data/annotated_context.ht"
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
        "chrx_nonpar", "chry_nonpar". Default is "autosome_par".
    :param test: Whether the Table is for testing purposes and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: TableResource of processed genomes, exomes, or context Table.
    """
    check_param_scope(version, genomic_region, data_type)
    return TableResource(
        f"{get_constraint_root(version, test)}/preprocessed_data/gnomad.v{version}.{data_type}.preprocessed.{genomic_region}.ht"
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
        "chrx_nonpar", or "chry_nonpar". Default is "autosome_par".
    :param test: Whether the Table is for testing purpose and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: TableResource of training dataset.
    """
    check_param_scope(version, genomic_region)
    return TableResource(
        f"{get_constraint_root(version, test)}/training_data/gnomad.v{version}.constraint_training.{genomic_region}.ht"
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
        f"{get_constraint_root(version, test)}/models/gnomad.v{version}.{model_type}.{genomic_region}.he"
    )


def get_predicted_proportion_observed_dataset(
    custom_vep_annotation: str = "transcript_consequences",
    version: str = CURRENT_VERSION,
    genomic_region: str = "autosome_par",
    test: bool = False,
) -> TableResource:
    """
    Return TableResource containing the expected variant counts and observed:expected ratio.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene").
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
        f"{get_constraint_root(version, test)}/predicted_proportion_observed/{custom_vep_annotation}/gnomad.v{version}.predicted_proportion_observed.{genomic_region}.ht"
    )


def get_constraint_metrics_dataset(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Return TableResource of pLI scores, observed:expected ratio, 90% confidence interval around the observed:expected ratio, and z scores.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param test: Whether the Table is for testing purposes and only contains sites in
        chr20, chrX, and chrY. Default is False.
    :return: TableResource of constraint metrics.
    """
    check_param_scope(version=version)

    return TableResource(
        f"{get_constraint_root(version, test)}/metrics/gnomad.v{version}.constraint_metrics.ht"
    )


def check_param_scope(
    version: Optional[str] = None,
    genomic_region: Optional[str] = None,
    data_type: Optional[str] = None,
    model_type: Optional[str] = None,
    custom_vep_annotation: Optional[str] = None,
) -> Union[str, None]:
    """
    Check if the specified version, genomic region, and data type are in the scope of the constraint pipeline.

    If version is specified, return the genome build of the version as a string.

    :param version: One of the release versions (`VERSIONS`). Default is None.
    :param genomic_region: The genomic region of the resource. One of "autosome_par",
        "chrx_non_par", or "chry_non_par". Default is None.
    :param data_type: One of "exomes", "genomes" or "context". Default is None.
    :param model_type: One of "plateau", "coverage". Default is None.
    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene"). Default is None.
    :return: Genome build of version as a string or None.
    """
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
    if version and version not in VERSIONS:
        raise ValueError("The requested version doesn't exist!")
    else:
        if version.startswith("2"):
            return "GRCh37"
        else:
            return "GRCh38"


def get_logging_path(name: str, version: str = CURRENT_VERSION) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: Output log path.
    """
    return f"{get_constraint_root(version, test=True)}/logging/{name}.log"


def get_checkpoint_path(
    name: str, version: str = CURRENT_VERSION, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable.
    :param version: Version of path to return.
    :param bool mt: Whether path is for a MatrixTable. Default is False.
    :return: Output checkpoint path.
    """
    return (
        f'{get_constraint_root(version, test=True)}/checkpoint_files/{name}.{"mt" if mt else "ht"}'
    )
