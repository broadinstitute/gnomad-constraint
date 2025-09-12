"""Resource definitions for the proemis3d pipeline."""

import logging
from typing import Optional

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
logger = logging.getLogger("proemis3d_pipeline")
logger.setLevel(logging.INFO)

VERSIONS = ["2.1.1", "4.1"]
"""Possible gnomAD versions for the proemis3d pipeline."""

CURRENT_VERSION = "4.1"
"""Current gnomAD version for the proemis3d pipeline."""

GENCODE_VERSION_MAP = {
    "2.1.1": "19",
    "4.1": "39",
}
"""GENCODE version map for each gnomAD version."""


def get_proemis3d_root(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Get root path to proemis3d resources.

    :param version: Version of proemis3d resources to use.
    :param test: Whether to use a tmp path for testing.
    :return: Root path to proemis3d resources.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/constraint/proemis3d"
        if test
        else f"gs://gnomad/v{version}/constraint/proemis3d"
    )


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
        f"gs://gnomad-julia/proemis3d/resources/"
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
        f"{get_proemis3d_root(version, test)}/preprocessed_data/"
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
    return TableResource(
        f"{get_proemis3d_root('2.1.1', test)}/preprocessed_data/af2.ht"
    )


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
        f"{get_proemis3d_root('2.1.1', test)}/preprocessed_data/af2_dist.ht"
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
        f"{get_proemis3d_root('2.1.1', test)}/preprocessed_data/af2_plddt.ht"
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
        f"{get_proemis3d_root('2.1.1', test)}/preprocessed_data/af2_pae.ht"
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
        f"{get_proemis3d_root(version, test)}/preprocessed_data/"
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
        f"{get_proemis3d_root(version, test)}/preprocessed_data/gencode_positions.ht"
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
    Get Proemis3D greedy algorithm Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Proemis3D greedy algorithm Hail Table resource.
    """
    return TableResource(f"{get_proemis3d_root(version, test)}/proemis3D_greedy.ht")


def get_forward_ht(
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get Proemis3D forward algorithm Hail Table resource.

    :param version: Version of gnomAD to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Proemis3D forward algorithm Hail Table resource.
    """
    return TableResource(f"{get_proemis3d_root(version, test)}/proemis3D_forward.ht")


def get_forward_annotation_ht(
    name: str,
    version: str = CURRENT_VERSION,
    test: bool = False,
) -> TableResource:
    """
    Get PROEMIS3D forward algorithm annotation Hail Table resource.

    :param name: Annotation type ('per_variant', 'per_missense_variant', 'per_residue',
        or 'per_proemis3d_region').
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: PROEMIS3D annotated forward algorithm Hail Table resource.
    """
    return TableResource(
        f"{get_proemis3d_root(version, test)}/proemis3D_forward.{name}.annotated.5_29_25.ht"
    )


def get_temp_all_possible_snvs_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get temp all possible SNVs Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Temp all possible SNVs Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad-tmp-4day/v{version}/constraint/proemis3d/all_possible_snvs.ht"
    )


def get_missense_viewer_input_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get missense viewer input Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Missense viewer input Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v{version}/constraint/proemis3d/missense_viewer_input.ht"
    )


########################################################################################
# The following functions are for external resources.
########################################################################################
def get_kaplanis_variants_tsv() -> str:
    """
    Get Kaplanis annotated variants file path.

    :return: File path to Kaplanis annotated variants.
    """
    return "gs://gnomad/v4.1/constraint/resources/variant_lists/tsv/kaplanis_variants_annotated_2024-05-15.txt"


def get_kaplanis_sig_variants_tsv() -> str:
    """
    Get Kaplanis significant variants file path.

    :return: File path to Kaplanis significant variants.
    """
    return "gs://gnomad/v4.1/constraint/resources/variant_lists/tsv/kaplanis_variants_sig.txt"


def get_kaplanis_variants_ht(
    liftover_to_grch38: bool = False,
    key_by_transcript: bool = False,
) -> TableResource:
    """
    Get processed Hail Table path for Kaplanis annotated de novo missense variants.

    This is the output of the `process_kaplanis_variants_ht` function, containing
    lifted GRCh38 loci and relevant annotations.

    :param liftover_to_grch38: Whether to liftover the variants to GRCh38. Default is
        False.
    :param key_by_transcript: Whether to key the table by transcript. Default is False.
    :return: Hail Table resource path for processed Kaplanis missense variants.
    """
    postfix = "liftover_to_grch38" if liftover_to_grch38 else ""
    postfix = f".{postfix}.keyed_by_transcript" if key_by_transcript else postfix
    return TableResource(
        f"gs://gnomad/v4.1/constraint/resources/variant_lists/ht/kaplanis_variants{postfix}.ht"
    )


def get_fu_variants_tsv() -> str:
    """
    Get Fu variants TSV file path.

    :return: Fu variants TSV file path.
    """
    return "gs://gnomad/v4.1/constraint/resources/variant_lists/tsv/fu_2022_supp20.txt"


def get_fu_variants_ht() -> TableResource:
    """
    Get processed Hail Table path for Fu annotated de novo missense variants.

    This is the output of the `process_fu_variants_ht` function, containing
    lifted GRCh38 loci and relevant annotations.

    :param key_by_transcript: Whether to key the table by transcript. Default is False.
    :return: Hail Table resource path for processed Fu missense variants.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/variant_lists/ht/fu_variants.ht"
    )


def get_interpro_annotations() -> str:
    """
    Get Ensembl BioMart export InterPro annotations file path.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: File path to InterPro annotations.
    """
    return "gs://gnomad/v4.1/constraint/resources/annotations/tsv/ensembl_biomart_export_interpro.txt"


def get_interpro_annotations_ht() -> TableResource:
    """
    Get InterPro annotations Hail Table resource.

    :return: InterPro annotations Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/annotations/ht/ensembl_biomart_export_interpro.ht"
    )


def get_cosmis_score_tsv(model: str) -> str:
    """
    Get COSMIS scores TSV path for a specified structure model.

    :param model: Structure model source ('alphafold', 'swiss_model', or 'pdb').
    :return: COSMIS scores TSV file path.
    """
    return f"gs://gnomad/v4.1/constraint/resources/3d_missense_methods/tsv/cosmis_scores_{model}.tsv.gz"


def get_cosmis_score_ht(model: str) -> TableResource:
    """
    Get COSMIS Hail Table resource for a specified structure model.

    :param model: Structure model source ('alphafold', 'swiss_model', or 'pdb').
    :return: COSMIS scores Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v4.1/constraint/resources/3d_missense_methods/ht/cosmis_scores_{model}.ht"
    )


def get_varity_tsv() -> str:
    """
    Get Varity TSV file path.

    :return: Varity TSV file path.
    """
    return "gs://gnomad/v4.1/constraint/resources/3d_missense_methods/tsv/varity_all_predictions.txt"


def get_varity_ht() -> TableResource:
    """
    Get Varity Hail Table resource.

    :return: Varity Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/3d_missense_methods/ht/varity_all_predictions.ht"
    )


def get_mtr3d_tsv() -> str:
    """
    Get MTR3D TSV file path.

    :return: MTR3D TSV file path.
    """
    return "gs://gnomad/v4.1/constraint/resources/3d_missense_methods/tsv/mtr_data.csv"


def get_mtr3d_ht() -> TableResource:
    """
    Get MTR3D Hail Table resource.

    :return: MTR3D Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/3d_missense_methods/ht/mtr3d_data.ht"
    )


def get_mtr_tsv() -> str:
    """
    Get MTR TSV file path.

    :return: MTR TSV file path.
    """
    return "gs://gnomad/v4.1/constraint/resources/annotations/tsv/full_MTR_scores.tsv"


def get_mtr_ht() -> TableResource:
    """
    Get MTR Hail Table resource.

    :return: MTR Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/annotations/ht/mtr_data.ht"
    )


def get_rmc_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get Hail Table resource with all regional missense constraint (RMC) scores.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: RMC Hail Table resource.
    """
    return TableResource(f"gs://gnomad/v{version}/constraint/resources/all_rmc.ht")


def get_rmc_browser_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get Hail Table resource with all regional missense constraint (RMC) scores.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: RMC Hail Table resource.
    """
    return TableResource(
        f"gs://regional_missense_constraint/constraint/{version}/2/rmc_browser.ht"
    )


def get_temp_processed_rmc_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get temp processed RMC Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Temp processed RMC Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad-tmp-4day/v{version}/constraint/resources/all_rmc.ht"
    )


def get_context_preprocessed_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get preprocessed context Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Context preprocessed Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad/v{version}/constraint_coverage_corrected/preprocessed_data/gnomad.v{version}.context.preprocessed.ht"
    )


def get_temp_context_preprocessed_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get temp preprocessed context Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Temp preprocessed context Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad-tmp-4day/v{version}/constraint_coverage_corrected/preprocessed_data/gnomad.v{version}.context.preprocessed.ht"
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


def get_temp_processed_constraint_ht(version: str = CURRENT_VERSION) -> TableResource:
    """
    Get temp processed constraint Hail Table resource.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: Temp processed constraint Hail Table resource.
    """
    return TableResource(
        f"gs://gnomad-tmp-4day/v{version}/constraint_coverage_corrected/gnomad.v{version}.constraint.ht"
    )


def get_clinvar_missense_ht() -> TableResource:
    """
    Get ClinVar missense Hail Table resource.

    :return: ClinVar missense Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/annotations/ht/clinvar_missense.ht"
    )


def get_revel_csv() -> str:
    """
    Get REVEL CSV file path.

    :return: REVEL CSV file path.
    """
    return "gs://gnomad-insilico/revel/revel-v1.3_all_chromosomes_with_transcript_ids.csv.bgz"


def get_insilico_annotations_ht(method: str) -> TableResource:
    """
    Get insilico annotations Hail Table resource.

    :param method: Insilico method to use. Must be one of 'cadd', 'phylop', or 'revel'.
    :return: Insilico annotations Hail Table resource.
    """
    if method == "revel":
        return TableResource(
            "gs://gnomad/v4.1/constraint/resources/annotations/ht/revel.ht"
        )
    else:
        if method not in {"cadd", "phylop"}:
            raise ValueError(
                f"Invalid method: {method}. Must be one of 'cadd' or 'phylop'."
            )

        return TableResource(
            f"gs://gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.{method}.grch38.ht"
        )


def get_genetics_gym_missense_scores_ht() -> TableResource:
    """
    Get Genetics Gym missense scores Hail Table resource.

    :return: Genetics Gym missense scores Hail Table resource.
    """
    return TableResource("gs://genetics-gym/vsm-tables/all-models-no-PAI3D.ht")


def get_processed_genetics_gym_missense_scores_ht() -> TableResource:
    """
    Get temp Genetics Gym missense scores Hail Table resource.

    :return: Temp Genetics Gym missense scores Hail Table resource.
    """
    return TableResource(
        "gs://gnomad/v4.1/constraint/resources/annotations/ht/all_missense_scores_percentile.with_uniprot.ht"
    )


def get_phaplo() -> ExpressionResource:
    """
    Get phaplo Hail expression resource.

    :return: Phaplo Hail expression resource.
    """
    return ExpressionResource(
        "gs://gnomad/v4.1/constraint/resources/gene_lists/he/phaplo_genes.he"
    )


def get_ptriplo() -> ExpressionResource:
    """
    Get ptriplo Hail expression resource.

    :return: Ptriplo Hail expression resource.
    """
    return ExpressionResource(
        "gs://gnomad/v4.1/constraint/resources/gene_lists/he/ptriplo_genes.he"
    )


def get_gnomad_de_novo_ht() -> TableResource:
    """
    Get gnomAD de novo Hail Table resource.

    :return: GnomAD de novo Hail Table resource.
    """
    return TableResource(
        "gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.de_novo.high_quality_coding.ht"
    )
