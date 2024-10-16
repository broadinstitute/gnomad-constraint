"""Script to perform the promis3d pipeline."""

import argparse
import logging

import hail as hl
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)

import gnomad_constraint.experimental.promis3d.resources as promis3d_res
from gnomad_constraint.experimental.promis3d.utils import (
    COLNAMES_TRANSLATIONS,
    convert_fasta_to_table,
    convert_gencode_transcripts_fasta_to_table,
    generate_codon_oe_table,
    get_gencode_positions,
    join_by_sequence,
    process_af2_structures,
    remove_multi_frag_uniprots,
    run_greedy,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("promis3d_pipeline")
logger.setLevel(logging.INFO)

TEST_TRANSCRIPT_ID = "ENST00000215754"
TEST_UNIPROT_ID = "P14174"
"""Transcript and UniProt IDs for testing."""


def get_promis3d_resources(
    version: str,
    overwrite: bool,
    test: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the promis3d pipeline.

    :param version: Version of promis3d resources to use.
    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :return: PipelineResourceCollection containing resources for all steps of the
        promis3D pipeline.
    """
    # Get glob for AlphaFold2 structures.
    af2_dir_path = promis3d_res.get_alpha_fold2_dir(version=version)
    if test:
        af2_dir_path = f"{af2_dir_path}/AF-{TEST_UNIPROT_ID}-*.cif.gz"
    else:
        af2_dir_path = f"{af2_dir_path}/*.cif.gz"

    # Initialize promis3D pipeline resource collection.
    promis3d_pipeline = PipelineResourceCollection(
        pipeline_name="promis3d",
        overwrite=overwrite,
        pipeline_resources={"AlphaFold2 directory": {"af2_dir_path": af2_dir_path}},
    )

    # Create resource collection for each step of the promis3D pipeline.
    gencode_transcipt = PipelineStepResourceCollection(
        "--convert-gencode-fastn-to-ht",
        input_resources={
            "GENCODE transcripts FASTA": {
                "gencode_transcipt_fasta_path": promis3d_res.get_gencode_fasta(
                    version=version, name="pc_transcripts"
                )
            },
        },
        output_resources={
            "gencode_transcipt_ht": promis3d_res.get_gencode_seq_ht(
                version=version, name="pc_transcripts", test=test
            ),
        },
    )
    gencode_translation = PipelineStepResourceCollection(
        "--convert-gencode-fasta-to-ht",
        input_resources={
            "GENCODE translations FASTA": {
                "gencode_translation_fasta_path": promis3d_res.get_gencode_fasta(
                    version=version, name="pc_translations"
                )
            },
        },
        output_resources={
            "gencode_translation_ht": promis3d_res.get_gencode_seq_ht(
                version=version, name="pc_translations", test=test
            ),
        },
    )
    read_af2_sequences = PipelineStepResourceCollection(
        "--read-af2-sequences",
        output_resources={"af2_ht": promis3d_res.get_af2_ht(version, test)},
    )
    compute_af2_distance_matrices = PipelineStepResourceCollection(
        "--compute-af2-distance-matrices",
        output_resources={"af2_dist_ht": promis3d_res.get_af2_dist_ht(version, test)},
    )
    gencode_alignment = PipelineStepResourceCollection(
        "--gencode-alignment",
        pipeline_input_steps=[gencode_translation, read_af2_sequences],
        output_resources={
            "matched_ht": promis3d_res.get_gencode_translations_matched_ht(
                version, test
            ),
        },
    )
    get_gencode_positions = PipelineStepResourceCollection(
        "--get-gencode-positions",
        pipeline_input_steps=[gencode_transcipt, gencode_alignment],
        add_input_resources={
            "GENCODE GTF": {"gencode_gtf_ht": promis3d_res.get_gencode_ht(version)}
        },
        output_resources={
            "gencode_pos_ht": promis3d_res.get_gencode_pos_ht(version, test),
        },
    )
    run_greedy = PipelineStepResourceCollection(
        "--run-greedy",
        pipeline_input_steps=[compute_af2_distance_matrices, get_gencode_positions],
        add_input_resources={
            "RMC OE Table": {"obs_exp_ht": promis3d_res.get_obs_exp_ht(version)}
        },
        output_resources={"greedy_ht": promis3d_res.get_greedy_ht(version, test)},
    )

    # Add all steps to the promis3D pipeline resource collection.
    promis3d_pipeline.add_steps(
        {
            "gencode_transcipt": gencode_transcipt,
            "gencode_translation": gencode_translation,
            "read_af2_sequences": read_af2_sequences,
            "compute_af2_distance_matrices": compute_af2_distance_matrices,
            "gencode_alignment": gencode_alignment,
            "get_gencode_positions": get_gencode_positions,
            "run_greedy": run_greedy,
        }
    )

    return promis3d_pipeline


def main(args):
    """Execute the Proemis 3D pipeline."""
    hl.init(
        log="/proemis_3d.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test = args.test
    overwrite = args.overwrite

    if version not in promis3d_res.VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = get_promis3d_resources(version, overwrite, test)

    if args.convert_gencode_fastn_to_ht:
        logger.info(
            "Importing and pre-process GENCODE transcripts FASTA file as a Hail Table."
        )
        res = resources.gencode_transcipt
        res.check_resource_existence()
        ht = convert_gencode_transcripts_fasta_to_table(
            res.gencode_transcipt_fasta_path
        )
        if test:
            ht = ht.filter(ht.enst == TEST_TRANSCRIPT_ID)
        ht = ht.checkpoint(res.gencode_transcipt_ht.path, overwrite=overwrite)
        ht.show()

    if args.convert_gencode_fasta_to_ht:
        logger.info(
            "Importing and pre-process GENCODE translations FASTA file as a Hail Table."
        )
        res = resources.gencode_translation
        res.check_resource_existence()
        ht = convert_fasta_to_table(
            res.gencode_translation_fasta_path, COLNAMES_TRANSLATIONS
        )
        if test:
            ht = ht.filter(ht.enst == TEST_TRANSCRIPT_ID)
        ht = ht.checkpoint(res.gencode_translation_ht.path, overwrite=overwrite)
        ht.show()

    if args.read_af2_sequences:
        logger.info(
            "Processing AlphaFold2 structures from a GCS bucket into a Hail Table."
        )
        res = resources.read_af2_sequences
        res.check_resource_existence()
        ht = process_af2_structures(resources.af2_dir_path)
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_ht.path, overwrite=overwrite)
        ht.show()

    if args.compute_af2_distance_matrices:
        logger.info("Computing distance matrices for AlphaFold2 structures.")
        res = resources.compute_af2_distance_matrices
        res.check_resource_existence()
        ht = process_af2_structures(resources.af2_dir_path, distance_matrix=True)
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_dist_ht.path, overwrite=overwrite)
        ht.show()

    if args.gencode_alignment:
        logger.info(
            "Joining the GENCODE translations and AlphaFold2 structures based on "
            "sequence."
        )
        res = resources.gencode_alignment
        res.check_resource_existence()
        ht = join_by_sequence(res.af2_ht.ht(), res.gencode_translation_ht.ht())
        ht = ht.checkpoint(res.matched_ht.path, overwrite=overwrite)
        ht.show()

    if args.get_gencode_positions:
        logger.info("Creating GENCODE positions Hail Table.")
        res = resources.get_gencode_positions
        res.check_resource_existence()
        ht = res.gencode_gtf_ht.ht()
        if test:
            ht.filter(ht.transcript_id == TEST_TRANSCRIPT_ID)
        ht = get_gencode_positions(
            res.gencode_transcipt_ht.ht(), res.matched_ht.ht(), ht
        )
        ht = ht.checkpoint(res.gencode_pos_ht.path, overwrite=overwrite)
        ht.show()

    if args.run_greedy:
        logger.info("Running greedy algorithm.")
        res = resources.run_greedy
        res.check_resource_existence()
        ht = res.obs_exp_ht.ht()
        if test:
            ht = ht.filter(ht.transcript == TEST_TRANSCRIPT_ID)
        ht = ht.filter(ht.annotation == "missense_variant")
        ht = ht.group_by("locus", "transcript").aggregate(
            obs=hl.agg.sum(ht.observed), exp=hl.agg.sum(ht.expected)
        )
        ht = generate_codon_oe_table(ht, res.gencode_pos_ht.ht())
        ht = ht.checkpoint(hl.utils.new_temp_file("codon_oe", "ht"))
        ht = run_greedy(res.af2_dist_ht.ht(), ht)
        ht = ht.checkpoint(res.greedy_ht.path, overwrite=overwrite)
        ht.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--version",
        help=(
            "Which version of the resource Tables will be used. Default is"
            f" {promis3d_res.CURRENT_VERSION}."
        ),
        type=str,
        default=promis3d_res.CURRENT_VERSION,
    )
    parser.add_argument(
        "--test",
        help="Whether to run a test instead of the full pipeline",
        action="store_true",
    )
    parser.add_argument(
        "--convert-gencode-fastn-to-ht",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--convert-gencode-fasta-to-ht",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--read-af2-sequences",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--compute-af2-distance-matrices",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--gencode-alignment",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--get-gencode-positions",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--run-greedy",
        help="",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
