"""Resource utility functions and resource definitions for the constraint pipeline."""

import logging
from typing import Optional, Union

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
    import_gencode,
)

from gnomad_constraint.resources.constants import (
    CURRENT_VERSION,
    CUSTOM_VEP_ANNOTATIONS,
    DATA_TYPES,
    EXTENSIONS,
    MODEL_TYPES,
    VERSIONS,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_pipeline")
logger.setLevel(logging.INFO)


def check_param_scope(
    version: Optional[str] = None,
    model_type: Optional[str] = None,
    custom_vep_annotation: Optional[str] = None,
    extension: Optional[str] = None,
    data_type: Optional[str] = None,
) -> Union[str, None]:
    """
    Check if the specified version, genomic region, and other parameters are in the scope of the constraint pipeline.

    If version is specified, return the genome build of the version as a string.

    :param version: One of the release versions (`VERSIONS`). Default is None.
    :param model_type: One of "plateau", "coverage". Default is None.
    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene"). Default is None.
    :param extension: File extension. Default is None.
    :param data_type: One of "exomes", "genomes". Default is None.
    :return: Genome build of version as a string or None.
    """
    if data_type and data_type not in DATA_TYPES:
        raise ValueError(f"data_type must be one of: {DATA_TYPES}!")
    if model_type and model_type not in MODEL_TYPES:
        raise ValueError(f"model_type must be one of: {MODEL_TYPES}!")
    if custom_vep_annotation and custom_vep_annotation not in CUSTOM_VEP_ANNOTATIONS:
        raise ValueError(
            f"custom_vep_annotation must be one of: {CUSTOM_VEP_ANNOTATIONS}!"
        )
    if extension and extension not in EXTENSIONS:
        raise ValueError(f"extension must be one of: {EXTENSIONS}!")
    if version:
        if version not in VERSIONS:
            raise ValueError("The requested version doesn't exist!")
        else:
            if version.startswith("2"):
                return "GRCh37"
            else:
                return "GRCh38"


def get_vep_context_ht(version: str) -> TableResource:
    """
    Return VEP context Table corresponding to specified gnomAD version.

    gnomAD version 2.1.1 uses build GRCh37 and VEP version 85. gnomAD v4 versions use build GRCh38 and VEP version 105.

    :param version: Version of gnomAD.
    :return: VEP context Table Resource.
    """
    if version == "2.1.1":
        return ref_grch37.vep_context.versions["85"]
    elif int(version[0]) == 4:
        return ref_grch38.vep_context.versions["105"]
    else:
        raise ValueError("Not a valid gnomAD version -- must be either 2.1.1 or 4.x!")


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
    elif int(version[0]) == 4:
        # Continue to use v3.1.2 for genomes as downsamplings are dropped in v4
        # versions.
        if data_type == "genomes":
            return gnomad_grch38.public_release(data_type).versions["3.1.2"]
        else:
            return gnomad_grch38.public_release(data_type).versions[version]
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
        return ref_grch38.methylation_sites
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
        return gnomad_grch38.coverage(data_type)


def get_gencode_ht(version: str) -> hl.Table:
    """
    Retrieve GENCODE Table with transcript version annotations.

    Re-imports the GENCODE GTF with ``include_version=True`` so that both
    ``transcript_id_version`` and ``gene_id_version`` fields are present,
    then checkpoints the result so subsequent calls read from the checkpoint.

    :param version: gnomAD version. If version 2, GENCODE v19 will be
        loaded. If version 4, GENCODE v39 will be re-imported with version
        fields and checkpointed.
    :return: Table of GENCODE data with version annotations.
    """
    if int(version[0]) == 2:
        return ref_grch37.gencode.ht()
    elif int(version[0]) == 4:
        gencode_resource = ref_grch38.gencode
        import_args = gencode_resource.versions[
            gencode_resource.default_version
        ].import_args
        ht = import_gencode(**import_args, include_version=True)
        checkpoint_path = "gs://gnomad-tmp/gencode_v39_with_versions.ht"

        return ht.checkpoint(checkpoint_path, _read_if_exists=True)
    else:
        raise ValueError("Version must be within gnomAD v2 or v4.")


def get_gencode_cds_ht(
    version: str = CURRENT_VERSION,
) -> TableResource:
    """Build and checkpoint a per-locus GENCODE CDS transcript ID table.

    Calls :func:`get_gencode_ht` to retrieve the GENCODE table, filters to
    CDS features, explodes each CDS interval into individual locus positions,
    and groups by locus to produce an array of transcript IDs per position.
    The result is checkpointed (read if it already exists) and returned as a
    :class:`TableResource`.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of GENCODE CDS positions, keyed by locus with
        ``transcript_id`` (array of transcript IDs whose CDS covers that
        position).
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version, temp=True)
    path = f"{root}/gencode_cds_positions.ht"

    gencode_ht = get_gencode_ht(version)
    gencode_ht = gencode_ht.filter(gencode_ht.feature == "CDS").select("transcript_id")
    gencode_ht = gencode_ht.annotate(
        positions=hl.range(
            gencode_ht.interval.start.position,
            gencode_ht.interval.end.position + 1,
        )
    ).explode("positions")
    gencode_ht = gencode_ht.key_by(
        locus=hl.locus(
            gencode_ht.interval.start.contig,
            gencode_ht.positions,
            reference_genome="GRCh38",
        )
    ).select("transcript_id")
    gencode_ht = gencode_ht.group_by("locus").aggregate(
        transcript_id=hl.agg.collect(gencode_ht.transcript_id)
    )
    gencode_ht.checkpoint(path, _read_if_exists=True)

    return TableResource(path)


def get_constraint_root(
    version: str = CURRENT_VERSION,
    test: bool = False,
    post_fix: Optional[str] = None,
    temp: bool = False,
    sub_dir: Optional[str] = None,
) -> str:
    """
    Return path to constraint root folder.

    :param version: Version of constraint path to return. Default is CURRENT_VERSION.
    :param test: Whether to use a tmp path. Default is False.
    :param post_fix: Postfix to append to the path. Default is None.
    :param temp: Whether to use a temp path. Default is False.
    :param sub_dir: Subdirectory to append to the path. Default is None.
    :return: Root path to constraint resources folder.
    """
    post_fix = post_fix or ""
    if post_fix:
        post_fix = f"_{post_fix}"

    sub_dir = sub_dir or ""
    if sub_dir:
        sub_dir = f"/{sub_dir}"

    constraint_dir = f"constraint{post_fix}{sub_dir}"

    if test:
        return f"gs://gnomad-tmp/gnomad_v{version}_testing/{constraint_dir}"
    if temp:
        return f"gs://gnomad-tmp/gnomad_v{version}/{constraint_dir}"

    return f"gs://gnomad/v{version}/{constraint_dir}"


def get_constraint_data(
    name: str,
    version: str = CURRENT_VERSION,
    test: bool = False,
    directory_post_fix: Optional[str] = None,
    sub_dir: Optional[str] = None,
    custom_vep_annotation: Optional[str] = None,
    extension: str = "ht",
    path_post_fix: Optional[str] = None,
    temp: bool = False,
) -> Union[TableResource, str, ExpressionResource]:
    """
    Return path, TableResource, or ExpressionResource of requested constraint data.

    :param name: Name of the constraint data to retrieve.
    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :param test: Whether the Table is for testing purpose and only contains a subset of
        the data. Default is False.
    :param directory_post_fix: Postfix to append to the root path. Default is None.
    :param sub_dir: Subdirectory to append to the path. Default is None.
    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene"). Default is
        None.
    :param extension: File extension. Default is "ht".
    :param path_post_fix: Postfix to append to the file name. Default is None.
    :return: Path, TableResource, or ExpressionResource of the constraint data.
    """
    check_param_scope(
        version, custom_vep_annotation=custom_vep_annotation, extension=extension
    )

    if custom_vep_annotation:
        sub_dir = f"{sub_dir}/" if sub_dir else ""
        sub_dir = f"{sub_dir}{custom_vep_annotation}"

    path_post_fix = path_post_fix or ""
    if path_post_fix:
        path_post_fix = f".{path_post_fix}"

    root_dir = get_constraint_root(
        version=version,
        test=test,
        post_fix=directory_post_fix,
        sub_dir=sub_dir,
        temp=temp,
    )
    path = f"{root_dir}/gnomad.v{version}.{name}{path_post_fix}.{extension}"

    if extension == "ht":
        return TableResource(path)
    if extension in {"tsv", "tsv.bgz", "log"}:
        return path
    if extension == "he":
        return ExpressionResource(path)


def get_mutation_ht(**kwargs) -> TableResource:
    """
    Return mutation Table that includes the baseline mutation rate for each substitution and context.

    :return: Mutation rate Table.
    """
    return get_constraint_data("mutation_rate", sub_dir="mutation_rate", **kwargs)


def get_release_mutation_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Return TableResource for the release mutation rate Table.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of the release mutation rate Table.
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version)
    return TableResource(f"{root}/release/gnomad.v{version}.mutation_rate.ht")


def get_release_constraint_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Return TableResource for the release constraint metrics Table.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of the release constraint metrics Table.
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version)
    return TableResource(f"{root}/release/gnomad.v{version}.constraint_metrics.ht")


def get_release_constraint_tsv_path(version: str = CURRENT_VERSION) -> str:
    """
    Return path for the release constraint metrics TSV.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: Path of the release constraint metrics TSV.
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version)
    return f"{root}/release/gnomad.v{version}.constraint_metrics.tsv.bgz"


def get_annotated_context_ht(**kwargs) -> TableResource:
    """
    Return TableResource of annotated context Table.

    :return: TableResource of annotated context Table.
    """
    return get_constraint_data("annotated_context", temp=True, **kwargs)


def get_preprocessed_ht(**kwargs) -> TableResource:
    """
    Return TableResource of preprocessed genome, exomes, and context Table.

    :return: TableResource of processed context Table.
    """
    return get_constraint_data(
        "context.preprocessed", sub_dir="preprocessed_data", **kwargs
    )


def get_training_dataset(**kwargs) -> TableResource:
    """
    Return TableResource of training dataset with observed and possible variant count.

    :return: TableResource of training dataset.
    """
    return get_constraint_data("constraint_training", sub_dir="training_data", **kwargs)


def get_training_tsv_path(**kwargs) -> str:
    """
    Return tsv of training dataset with observed and possible variant count.

    :return: TSV path of training dataset.
    """
    return get_constraint_data(
        "constraint_training", sub_dir="training_data", extension="tsv.bgz", **kwargs
    )


def get_models(model_type: str, **kwargs) -> ExpressionResource:
    """
    Return path to a HailExpression that contains desired model type.

    :param model_type: The type of model. One of "plateau", "coverage". Default is None.
    :return: Path to the specified model.
    """
    check_param_scope(model_type=model_type)
    return get_constraint_data(model_type, sub_dir="models", extension="he", **kwargs)


def get_per_variant_expected_dataset(
    custom_vep_annotation: str = "transcript_consequences", **kwargs
) -> TableResource:
    """
    Return TableResource containing the expected variant counts and observed:expected ratio.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene").
    :return: Path of the model.
    """
    return get_constraint_data(
        "per_variant_expected",
        sub_dir="apply_models",
        custom_vep_annotation=custom_vep_annotation,
        **kwargs,
    )


def get_aggregated_per_variant_expected(
    custom_vep_annotation: str = "transcript_consequences", **kwargs
) -> TableResource:
    """
    Return TableResource containing the expected variant counts and observed:expected ratio.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene").
    :return: Path of the model.
    """
    return get_constraint_data(
        "per_variant_expected.aggregated",
        sub_dir="apply_models",
        custom_vep_annotation=custom_vep_annotation,
        **kwargs,
    )


def get_constraint_group_ht(custom_vep_annotation: str, **kwargs) -> TableResource:
    """
    Return TableResource of constraint group Table.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene").
    :return: TableResource of constraint group Table.
    """
    return get_constraint_data(
        "constraint_group",
        sub_dir="apply_models",
        custom_vep_annotation=custom_vep_annotation,
        **kwargs,
    )


def get_constraint_metrics_dataset(
    custom_vep_annotation: str = "transcript_consequences", **kwargs
) -> TableResource:
    """
    Return TableResource of pLI scores, observed:expected ratio, 90% confidence interval around the observed:expected ratio, and z scores.

    :param custom_vep_annotation: The VEP annotation used to customize the constraint
        model (one of "transcript_consequences" or "worst_csq_by_gene").
    :return: TableResource of constraint metrics.
    """
    return get_constraint_data(
        "constraint_metrics",
        sub_dir="metrics",
        custom_vep_annotation=custom_vep_annotation,
        **kwargs,
    )


def get_gene_quality_metrics_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Return TableResource of per-transcript gene quality metrics.

    Contains coverage and mapping quality metrics per transcript used for
    annotating the release constraint Table.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: TableResource of gene quality metrics.
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version)
    return TableResource(f"{root}/metrics/gnomad.v{version}.gene_quality_metrics.ht")


def get_release_downsampling_tsv_path(version: str = CURRENT_VERSION) -> str:
    """
    Return path for the release per-genetic-ancestry downsampling TSV.

    :param version: One of the release versions (`VERSIONS`). Default is
        `CURRENT_VERSION`.
    :return: Path of the release downsampling TSV.
    """
    check_param_scope(version=version)
    root = get_constraint_root(version=version)
    return f"{root}/release/gnomad.v{version}.constraint_metrics.downsampling.tsv.bgz"


def get_logging_path(name: str, **kwargs) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :return: Output log path.
    """
    return get_constraint_data(
        name, sub_dir="logging", extension="log", test=True, **kwargs
    )


def get_checkpoint_path(name: str, **kwargs) -> TableResource:
    """
    Create a checkpoint TableResource.

    :param name: Name of intermediate Table.
    :return: Output checkpoint TableResource.
    """
    return get_constraint_data(name, sub_dir="checkpoint_files", test=True, **kwargs)
