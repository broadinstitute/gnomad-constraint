"""Resource definitions for the promis3d pipeline."""

import logging

import hail as hl
from gnomad.resources.grch37.reference_data import gencode as grch37_gencode
from gnomad.resources.grch38.reference_data import gencode as grch38_gencode
from gnomad.resources.resource_utils import (
    BaseResource,
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("promis3d_pipeline")
logger.setLevel(logging.INFO)

VERSIONS = ["2.1.1", "4.1"]
"""Possible gnomAD versions for the promis3d pipeline."""

CURRENT_VERSION = "4.1"
"""Current gnomAD version for the promis3d pipeline."""

GENCODE_VERSION_MAP = {
    "2.1.1": "19",
    "4.1": "39",
}
"""GENCODE version map for each gnomAD version."""


def get_promis3d_root(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Get root path to promis3d resources.

    :param version: Version of promis3d resources to use.
    :param test: Whether to use a tmp path for testing.
    :return: Root path to promis3d resources.
    """
    return "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run"
    # return (
    #    f"gs://gnomad-tmp/gnomad_v{version}_testing/constraint/promis3d"
    #    if test
    #    else f"gs://gnomad/v{version}/constraint/promis3d"
    # )


def get_gencode_fasta(
    version: str = CURRENT_VERSION,
    name: str = "pc_transcripts",
) -> str:
    """
    Get GENCODE FASTA file path.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param name: Name of the type of GENCODE FASTA file to get. One of 'pc_transcripts',
        'pc_translations'. Default is 'pc_transcripts'.
    :return: GENCODE FASTA file path.
    """
    # TODO: Change this path when moved to a more permanent location.
    return (
        f"gs://gnomad-julia/promis3d/resources/"
        f"gencode.v{GENCODE_VERSION_MAP[version]}.{name}.fa.gz"
    )


def get_alpha_fold2_dir(version: str = CURRENT_VERSION) -> str:
    """
    Get AlphaFold2 directory path.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :return: AlphaFold2 directory path.
    """
    # TODO: Change this path when moved to a more permanent location and add any
    #  needed versioning.
    return "gs://gnomad-julia/alphafold2"


def get_gencode_seq_ht(
    version: str = CURRENT_VERSION,
    name: str = "pc_transcripts",
    test: bool = False,
) -> TableResource:
    """
    Get GENCODE sequences Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param name: Name of the type of GENCODE sequences to get. One of 'pc_transcripts',
        'pc_translations'. Default is 'pc_transcripts'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE sequences Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/"
        f"gencode_sequences.{name}.ht"
    )


def get_af2_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get alphafold2 Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE sequences Hail Table resource.
    """
    return TableResource(f"{get_promis3d_root('2.1.1', test)}/preprocessed_data/af2.ht")


def get_af2_dist_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get alphafold2 distance matrix Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: AlphaFold2 distance matrix Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root('2.1.1', test)}/preprocessed_data/af2_dist.ht"
    )


def get_af2_plddt_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get alphafold2 pLDDT Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: AlphaFold2 pLDDT Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root('2.1.1', test)}/preprocessed_data/af2_plddt.ht"
    )


def get_af2_pae_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get alphafold2 pAE Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: AlphaFold2 pAE Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root('2.1.1', test)}/preprocessed_data/af2_pae.ht"
    )


def get_gencode_translations_matched_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get GENCODE translations matched Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE translations matched Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/"
        f"gencode_translations_matched.ht"
    )


def get_gencode_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get GENCODE Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :return: GENCODE Hail Table resource.
    """
    if version == "2.1.1":
        gencode = grch37_gencode
    elif version == "4.1":
        gencode = grch38_gencode
    else:
        raise ValueError(f"Invalid version: {version}")

    return gencode.versions[f"v{GENCODE_VERSION_MAP[version]}"]


def get_gencode_pos_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get GENCODE positions Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: GENCODE positions Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/preprocessed_data/gencode_positions.ht"
    )


def get_obs_exp_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get observed/expected Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :return: Observed/expected Hail Table resource.
    """
    # TODO: Change this path when moved to a more permanent location.
    if version == "2.1.1":
        return TableResource(
            "gs://gnomad/v2.1.1/constraint/temp/gnomad.v2.1.1.per_base_expected.ht"
        )
    elif version == "4.1":
        return TableResource(
            "gs://gnomad/v4.1/constraint_coverage_corrected/apply_models/transcript_consequences/gnomad.v4.1.per_variant_expected.coverage_corrected.ht"
        )
    else:
        raise ValueError(f"Invalid version: {version}")


def get_greedy_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get Promis3D greedy algorithm Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Promis3D greedy algorithm Hail Table resource.
    """
    return TableResource(f"{get_promis3d_root(version, test)}/promis3D_greedy.ht")


def get_forward_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
    name: str = "",
) -> TableResource:
    """
    Get Promis3D forward algorithm Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Promis3D forward algorithm Hail Table resource.
    """
    name = f".{name}" if name else ""
    return TableResource(
        f"{get_promis3d_root(version, test)}/output/promis3D_forward{name}.ht"
    )


def get_forward_annotation_ht(
    name: str,
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get PROMIS3D forward algorithm annotation Hail Table resource.

    :param name: Annotation type ('per_variant', 'per_missense_variant', 'per_residue',
        or 'per_promis3d_region').
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: PROMIS3D annotated forward algorithm Hail Table resource.
    """
    return TableResource(
        f"{get_promis3d_root(version, test)}/promis3D_forward.{name}.annotated.5_29_25.ht"
    )


########################################################################################
# The following functions are for external resources.
########################################################################################
def get_kaplanis_variants() -> str:
    """
    Get Kaplanis annotated variants file path.

    :return: File path to Kaplanis annotated variants.
    """
    return "gs://gnomad-julia/kaplanis_variants_annotated_2024-05-15.txt"


def get_kaplanis_sig_variants() -> str:
    """
    Get Kaplanis significant variants file path.

    :return: File path to Kaplanis significant variants.
    """
    return "gs://gnomad-julia/promis3d/kaplanis_variants_sig.txt"


def get_kaplanis_variants_annotated_ht() -> TableResource:
    """
    Get processed Hail Table path for Kaplanis annotated de novo missense variants.

    This is the output of the `process_kaplanis_variants_ht` function, containing
    lifted GRCh38 loci and relevant annotations.

    :return: Hail Table resource path for processed Kaplanis missense variants.
    """
    return TableResource("gs://gnomad-julia/promis3d/kaplanis_variants_annotated.ht")


def get_interpro_annotations(version: str = CURRENT_VERSION) -> str:
    """
    Get Ensembl BioMart export InterPro annotations file path.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: File path to InterPro annotations.
    """
    return f"gs://gnomad/v{version}/constraint/resources/ensembl_biomart_export_interpro.txt"


def get_cosmis_score_tsv(model: str, version: str = CURRENT_VERSION) -> str:
    """
    Get COSMIS scores TSV path for a specified structure model.

    :param model: Structure model source ('alphafold', 'swiss_model', or 'pdb').
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: COSMIS scores TSV file path.
    """
    return f"gs://gnomad/v{version}/constraint/resources/cosmis_scores_{model}.tsv.gz"


def get_cosmis_score_ht(model: str, version: str = CURRENT_VERSION) -> TableResource:
    """
    Get COSMIS Hail Table resource for a specified structure model.

    :param model: Structure model source ('alphafold', 'swiss_model', or 'pdb').
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: COSMIS scores Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v{version}/constraint/resources/cosmis_scores_{model}.ht"
    )


def get_all_rmc_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get Hail Table resource with all regional missense constraint (RMC) scores.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: RMC Hail Table resource.
    """
    return TableResource(f"gs://gnomad/v{version}/constraint/resources/all_rmc.ht")


def get_context_preprocessed_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get preprocessed context Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Context preprocessed Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v{version}/constraint_coverage_corrected/preprocessed_data/gnomad.v{version}.context.preprocessed.ht"
    )


def get_constraint_metrics_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get coverage-corrected constraint metrics Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Constraint metrics Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v{version}/constraint_coverage_corrected/metrics/transcript_consequences/gnomad.v{version}.constraint_metrics.coverage_corrected.ht"
    )
