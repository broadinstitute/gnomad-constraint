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
    join_by_sequence,
    process_af2_structures,
    remove_multi_frag_uniprots,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("promis3d_pipeline")
logger.setLevel(logging.INFO)


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
        promis3d pipeline.
    """
    # Initialize promis3d pipeline resource collection.
    promis3d_pipeline = PipelineResourceCollection(
        pipeline_name="promis3d",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the promis3d pipeline.
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
        input_resources={
            "AlphaFold2 directory": {
                "af2_dir_path": promis3d_res.get_alpha_fold2_dir(version=version),
            },
        },
        output_resources={
            "af2_ht": promis3d_res.get_af2_ht(version=version, test=test),
        },
    )
    gencode_alignment = PipelineStepResourceCollection(
        "--gencode-alignment",
        pipeline_input_steps=[gencode_translation, read_af2_sequences],
        output_resources={
            "matched_ht": promis3d_res.get_gencode_translations_matched_ht(
                version=version, test=test
            ),
        },
    )

    # Add all steps to the promis3d pipeline resource collection.
    promis3d_pipeline.add_steps(
        {
            "gencode_transcipt": gencode_transcipt,
            "gencode_translation": gencode_translation,
            "read_af2_sequences": read_af2_sequences,
            "gencode_alignment": gencode_alignment,
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
        ht = ht.checkpoint(res.gencode_translation_ht.path, overwrite=overwrite)
        ht.show()

    if args.read_af2_sequences:
        logger.info(
            "Processing AlphaFold2 structures from a GCS bucket into a Hail Table."
        )
        res = resources.read_af2_sequences
        res.check_resource_existence()
        ht = process_af2_structures(res.af2_dir_path)
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_ht.path, overwrite=overwrite)
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
        "--gencode-alignment",
        help="",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
