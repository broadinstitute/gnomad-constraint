"""Script to perform the proemis3d pipeline."""

import argparse
import logging

import hail as hl
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)

import gnomad_constraint.experimental.promis3d.resources as proemis3d_res
from gnomad_constraint.experimental.promis3d.constants import MIN_EXP_MIS
from gnomad_constraint.experimental.promis3d.utils import (
    COLNAMES_TRANSLATIONS,
    convert_fasta_to_table,
    convert_gencode_transcripts_fasta_to_table,
    create_missense_viewer_input_ht,
    create_per_proemis3d_region_ht_from_residue_ht,
    create_per_residue_ht_from_snv_ht,
    create_per_snv_combined_ht,
    determine_regions_with_min_oe_upper,
    generate_codon_oe_table,
    get_gencode_positions,
    join_by_sequence,
    process_af2_structures,
    remove_multi_frag_uniprots,
    run_forward,
    run_forward_no_catch_all,
    run_forward_no_catch_all_standardized,
    run_greedy,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("proemis3d_pipeline")
logger.setLevel(logging.INFO)

TEST_TRANSCRIPT_ID = "ENST00000372435"  # "ENST00000215754"
TEST_UNIPROT_ID = "P60891"  # "P14174"
"""Transcript and UniProt IDs for testing."""


def get_proemis3d_resources(
    version: str,
    overwrite: bool,
    test: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the proemis3d pipeline.

    :param version: Version of proemis3d resources to use.
    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :return: PipelineResourceCollection containing resources for all steps of the
        proemis3D pipeline.
    """
    # Get glob for AlphaFold2 structures.
    af2_dir_path = proemis3d_res.get_alpha_fold2_dir(version=version)
    if test:
        af2_struct_dir_path = f"{af2_dir_path}/AF-{TEST_UNIPROT_ID}-*.cif.gz"
    else:
        af2_struct_dir_path = f"{af2_dir_path}/*.cif.gz"

    # Get glob for AlphaFold2 confidence.
    if test:
        af2_conf_dir_path = (
            f"{af2_dir_path}/AF-{TEST_UNIPROT_ID}-*confidence_v4.json.gz"
        )
    else:
        af2_conf_dir_path = f"{af2_dir_path}/*confidence_v4.json.gz"

    # Get glob for AlphaFold2 pAE.
    if test:
        af2_pae_dir_path = (
            f"{af2_dir_path}/AF-{TEST_UNIPROT_ID}-*predicted_aligned_error_v4.json.gz"
        )
    else:
        af2_pae_dir_path = f"{af2_dir_path}/*predicted_aligned_error_v4.json.gz"

    # Initialize proemis3D pipeline resource collection.
    proemis3d_pipeline = PipelineResourceCollection(
        pipeline_name="proemis3d",
        overwrite=overwrite,
        pipeline_resources={
            "AlphaFold2 directory": {
                "af2_struct_dir_path": af2_struct_dir_path,
                "af2_conf_dir_path": af2_conf_dir_path,
                "af2_pae_dir_path": af2_pae_dir_path,
            }
        },
    )

    # Create resource collection for each step of the proemis3D pipeline.
    gencode_transcipt = PipelineStepResourceCollection(
        "--convert-gencode-fastn-to-ht",
        input_resources={
            "GENCODE transcripts FASTA": {
                "gencode_transcipt_fasta_path": proemis3d_res.get_gencode_fasta(
                    version=version, name="pc_transcripts"
                )
            },
        },
        output_resources={
            "gencode_transcipt_ht": proemis3d_res.get_gencode_seq_ht(
                version=version, name="pc_transcripts", test=test
            ),
        },
    )
    gencode_translation = PipelineStepResourceCollection(
        "--convert-gencode-fasta-to-ht",
        input_resources={
            "GENCODE translations FASTA": {
                "gencode_translation_fasta_path": proemis3d_res.get_gencode_fasta(
                    version=version, name="pc_translations"
                )
            },
        },
        output_resources={
            "gencode_translation_ht": proemis3d_res.get_gencode_seq_ht(
                version=version, name="pc_translations", test=test
            ),
        },
    )
    read_af2_sequences = PipelineStepResourceCollection(
        "--read-af2-sequences",
        output_resources={"af2_ht": proemis3d_res.get_af2_ht(version, test)},
    )
    compute_af2_distance_matrices = PipelineStepResourceCollection(
        "--compute-af2-distance-matrices",
        # output_resources={"af2_dist_ht": proemis3d_res.get_af2_dist_ht(version, test)},
        output_resources={"af2_dist_ht": proemis3d_res.get_af2_dist_ht(version)},
    )
    extract_af2_plddt = PipelineStepResourceCollection(
        "--extract-af2-plddt",
        output_resources={
            "af2_plddt_ht": proemis3d_res.get_af2_plddt_ht(version, test)
        },
    )
    extract_af2_pae = PipelineStepResourceCollection(
        "--extract-af2-pae",
        output_resources={"af2_pae_ht": proemis3d_res.get_af2_pae_ht(version, test)},
    )
    gencode_alignment = PipelineStepResourceCollection(
        "--gencode-alignment",
        pipeline_input_steps=[gencode_translation, read_af2_sequences],
        output_resources={
            "matched_ht": proemis3d_res.get_gencode_translations_matched_ht(
                version, test
            ),
        },
    )
    get_gencode_positions = PipelineStepResourceCollection(
        "--get-gencode-positions",
        pipeline_input_steps=[gencode_transcipt, gencode_alignment],
        add_input_resources={
            "GENCODE GTF": {"gencode_gtf_ht": proemis3d_res.get_gencode_ht(version)}
        },
        output_resources={
            # "gencode_pos_ht": promis3d_res.get_gencode_pos_ht(version, test),
            "gencode_pos_ht": proemis3d_res.get_gencode_pos_ht(version),
        },
    )
    run_greedy = PipelineStepResourceCollection(
        "--run-greedy",
        pipeline_input_steps=[compute_af2_distance_matrices, get_gencode_positions],
        add_input_resources={
            "RMC OE Table": {"obs_exp_ht": proemis3d_res.get_obs_exp_ht(version)}
        },
        output_resources={"greedy_ht": proemis3d_res.get_greedy_ht(version, test)},
    )
    run_forward = PipelineStepResourceCollection(
        "--run-forward",
        pipeline_input_steps=[compute_af2_distance_matrices, get_gencode_positions],
        add_input_resources={
            "RMC OE Table": {"obs_exp_ht": proemis3d_res.get_obs_exp_ht(version)}
        },
        output_resources={"forward_ht": proemis3d_res.get_forward_ht(version, test)},
    )
    # Add resources for per-variant, per-residue, and per-region HTs.
    write_per_variant = PipelineStepResourceCollection(
        "--write-per-variant",
        pipeline_input_steps=[
            gencode_transcipt,
            gencode_translation,
            gencode_alignment,
            compute_af2_distance_matrices,
            extract_af2_plddt,
            extract_af2_pae,
            run_forward,
        ],
        add_input_resources={
            "GENCODE GTF": {"gencode_gtf_ht": proemis3d_res.get_gencode_ht(version)}
        },
        output_resources={
            "per_variant_ht": proemis3d_res.get_forward_annotation_ht(
                "per_variant", version, test
            ),
        },
    )
    write_per_missense_variant = PipelineStepResourceCollection(
        "--write-per-missense-variant",
        pipeline_input_steps=[write_per_variant],
        output_resources={
            "per_missense_variant_ht": proemis3d_res.get_forward_annotation_ht(
                "per_missense_variant", version, test
            ),
        },
    )
    write_per_residue = PipelineStepResourceCollection(
        "--write-per-residue",
        pipeline_input_steps=[write_per_variant],
        output_resources={
            "per_residue_ht": proemis3d_res.get_forward_annotation_ht(
                "per_residue", version, test
            ),
        },
    )
    write_per_region = PipelineStepResourceCollection(
        "--write-per-region",
        pipeline_input_steps=[write_per_residue],
        output_resources={
            "per_region_ht": proemis3d_res.get_forward_annotation_ht(
                "per_region", version, test
            ),
        },
    )
    create_missense_viewer_input_ht = PipelineStepResourceCollection(
        "--create-missense-viewer-input-ht",
        pipeline_input_steps=[get_gencode_positions, run_forward],
        output_resources={
            "missense_viewer_input_ht": proemis3d_res.get_missense_viewer_input_ht(
                version
            ),
        },
    )

    # Add all steps to the proemis3D pipeline resource collection.
    proemis3d_pipeline.add_steps(
        {
            "gencode_transcipt": gencode_transcipt,
            "gencode_translation": gencode_translation,
            "read_af2_sequences": read_af2_sequences,
            "compute_af2_distance_matrices": compute_af2_distance_matrices,
            "extract_af2_plddt": extract_af2_plddt,
            "extract_af2_pae": extract_af2_pae,
            "gencode_alignment": gencode_alignment,
            "get_gencode_positions": get_gencode_positions,
            "run_greedy": run_greedy,
            "run_forward": run_forward,
            "write_per_variant": write_per_variant,
            "write_per_missense_variant": write_per_missense_variant,
            "write_per_residue": write_per_residue,
            "write_per_region": write_per_region,
        }
    )

    return proemis3d_pipeline


def main(args):
    """Execute the Proemis 3D pipeline."""
    hl.init(
        log="/proemis_3d.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    version = args.version
    test = args.test
    overwrite = args.overwrite

    if version not in proemis3d_res.VERSIONS:
        raise ValueError("The requested version of resource Tables is not available.")

    # Construct resources with paths for intermediate Tables generated in the pipeline.
    resources = get_proemis3d_resources(version, overwrite, test)

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
            res.gencode_translation_fasta_path, COLNAMES_TRANSLATIONS[version]
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
        ht = process_af2_structures(resources.af2_struct_dir_path, mode="sequence")
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_ht.path, overwrite=overwrite)
        ht.show()

    if args.compute_af2_distance_matrices:
        logger.info("Computing distance matrices for AlphaFold2 structures.")
        res = resources.compute_af2_distance_matrices
        res.check_resource_existence()
        ht = process_af2_structures(resources.af2_struct_dir_path, mode="distance")
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_dist_ht.path, overwrite=overwrite)
        ht.show()

    if args.extract_af2_plddt:
        logger.info("Extracting pLDDT scores from AlphaFold2 structures.")
        res = resources.extract_af2_plddt
        res.check_resource_existence()
        ht = process_af2_structures(resources.af2_conf_dir_path, mode="plddt")
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_plddt_ht.path, overwrite=overwrite)
        ht.show()

    if args.extract_af2_pae:
        logger.info("Extracting pAE scores from AlphaFold2 structures.")
        res = resources.extract_af2_pae
        res.check_resource_existence()
        ht = process_af2_structures(resources.af2_pae_dir_path, mode="pae")
        ht = remove_multi_frag_uniprots(ht)
        ht = ht.checkpoint(res.af2_pae_ht.path, overwrite=overwrite)
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

    if True:
        # if args.run_greedy or args.run_forward:
        # logger.info("Preparing to run greedy and/or forward algorithms.")
        # if args.run_greedy:
        #    res = resources.run_greedy
        # res.check_resource_existence()
        # if args.run_forward:
        #    res = resources.run_forward
        # res.check_resource_existence()

        res = resources.run_forward

        # Use new shuffle method for apply models to prevent shuffle errors.
        hl._set_flags(use_new_shuffle="1")

        # ht = res.obs_exp_ht.ht()
        # gencode_pos_ht = hl.read_table("gs://gnomad/v4.1/constraint/promis3d/preprocessed_data/gencode_positions.ht")#res.gencode_pos_ht.ht()
        # af2_dist_ht = hl.read_table("gs://gnomad/v2.1.1/constraint/promis3d/preprocessed_data/af2_dist.ht") #res.af2_dist_ht.ht()
        # test_transcripts = hl.set(["ENST00000644876", "ENST00000233146", "ENST00000371953", "ENST00000347132", "ENST00000357654"])
        # ht = ht.filter(test_transcripts.contains(ht.transcript))
        # gencode_pos_ht = gencode_pos_ht.filter(test_transcripts.contains(gencode_pos_ht.enst)).cache()
        # uniprot_ids = hl.set(gencode_pos_ht.aggregate(hl.agg.collect_as_set(gencode_pos_ht.uniprot_id)))
        # af2_dist_ht = af2_dist_ht.filter(uniprot_ids.contains(af2_dist_ht.uniprot_id))

        # if test:
        #    genes_for_testing_ht = hl.import_table(
        #        "gs://gnomad/v4.1/constraint/resources/genes_for_testing.txt"
        #    ).cache()
        #    test_transcripts = hl.set(genes_for_testing_ht.transcript_id.collect())
        #    ht = ht.filter(test_transcripts.contains(ht.transcript))
        #    gencode_pos_ht = gencode_pos_ht.filter(test_transcripts.contains(gencode_pos_ht.enst)).cache()
        #    uniprot_ids = hl.set(gencode_pos_ht.aggregate(hl.agg.collect_as_set(gencode_pos_ht.uniprot_id)))
        #    af2_dist_ht = af2_dist_ht.filter(uniprot_ids.contains(af2_dist_ht.uniprot_id))
        #    #ht = ht.filter(ht.transcript == TEST_TRANSCRIPT_ID)

        # ht = ht.filter(ht.annotation == "missense_variant")

        # if version == "2.1.1":
        #    ht = ht.group_by("locus", "transcript").aggregate(
        #        obs=hl.agg.sum(ht.observed), exp=hl.agg.sum(ht.expected)
        #    )
        # elif version == "4.1":
        # ht = ht.group_by("locus", "transcript").aggregate(
        #        obs=hl.agg.sum(ht.calibrate_mu.observed_variants[0]),
        #        exp=hl.agg.sum(ht.expected_variants[0]),
        # )
        # else:
        #    raise ValueError(
        #        "The requested version of the resource Tables is not available."
        #    )
        # ht = generate_codon_oe_table(ht, gencode_pos_ht)
        # ht = ht.repartition(5000).checkpoint(hl.utils.new_temp_file("codon_oe", "ht"))
        # ht = ht.repartition(1).checkpoint(
        #    "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/codon_oe.ht"
        # )
        # af2_ht = (
        #    af2_dist_ht
        #    #.repartition(5000)
        #    .repartition(50)
        #    .repartition(1)
        #    .checkpoint(
        #        "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/af2_dist.ht",
        #        overwrite=True
        #    )
        # )
        # af2_ht = hl.read_table("gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/af2_dist.ht")
        af2_ht = hl.read_table(
            "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/af2_dist.ht"
        )
        # ht = determine_regions_with_min_oe_upper(
        #    af2_ht, ht, min_exp_mis=args.min_exp_mis
        # )
        # ht = ht.repartition(5000).checkpoint(
        # ht = ht.repartition(50).checkpoint(
        #    hl.utils.new_temp_file("sort_regions_by_oe", "ht")
        # )
        # ht = hl.read_table("gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/sort_regions_by_oe.ht")
        # af2_ht = hl.read_table(
        #    "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/af2_dist.ht"
        # )
        # ht = hl.read_table(
        #    "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/codon_oe.ht"
        # )
        ht = hl.read_table(
            "gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/codon_oe.ht"
        )

        plddt_out = ""
        if args.plddt_cutoff is not None:
            # plddt_ht = promis3d_res.get_af2_plddt_ht("2.1.1", args.test).ht()
            plddt_ht = hl.read_table(
                "gs://gnomad/v2.1.1/constraint/promis3d/preprocessed_data/af2_plddt.ht"
            )
            # print(plddt_ht.filter(plddt_ht.uniprot_id == "A0A024R2K8").collect())
            ht = ht.annotate(
                oe_by_transcript=ht.oe_by_transcript.map(
                    lambda x: x.annotate(
                        oe=hl.zip(x.oe, plddt_ht[ht.uniprot_id].plddt).map(
                            lambda y: hl.or_missing(y[1] >= args.plddt_cutoff, y[0])
                        )
                    )
                )
            )
            # print(ht.filter(ht.uniprot_id == "A0A024R2K8").collect())
            plddt_out = f".plddt_cutoff_{args.plddt_cutoff}"

        # if args.run_greedy:
        #    logger.info("Running greedy algorithm.")
        #    res = resources.run_greedy
        #    greedy_ht = run_greedy(ht)
        #    greedy_ht = greedy_ht.checkpoint(res.greedy_ht.path, overwrite=overwrite)
        #    greedy_ht.show()

        min_exp_mis_out = (
            ""
            if args.min_exp_mis == MIN_EXP_MIS
            else f".min_exp_mis_{args.min_exp_mis}"
        )

        if args.run_forward_original:
            # logger.info("Running forward algorithm.")
            # res = resources.run_forward
            ht = determine_regions_with_min_oe_upper(
                af2_ht, ht, min_exp_mis=args.min_exp_mis, oe_upper_method="chisq"
            )
            # ht = ht.repartition(200).checkpoint(
            #    f"gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/sort_regions_by_oe.min_exp_mis_{args.min_exp_mis}.chisq.ht",
            #    _read_if_exists=True,
            # )
            ht = ht.repartition(1).checkpoint(
                f"gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/sort_regions_by_oe.min_exp_mis_{args.min_exp_mis}.chisq.ht",
                _read_if_exists=True,
            )
            output_path = promis3d_res.get_forward_ht(
                name=f"oe_upper_chisq{min_exp_mis_out}{plddt_out}"
            ).path
            forward_ht = run_forward(
                ht, min_exp_mis=args.min_exp_mis, oe_upper_method="chisq"
            )
            # forward_ht = forward_ht.checkpoint(res.forward_ht.path, overwrite=overwrite)
            forward_ht = forward_ht.checkpoint(output_path, overwrite=overwrite)
            forward_ht.show()

        ht = determine_regions_with_min_oe_upper(
            af2_ht, ht, min_exp_mis=args.min_exp_mis, oe_upper_method="gamma"
        )
        # ht = ht.repartition(200).checkpoint(
        #    f"gs://gnomad/v4.1/constraint/promis3d/test_gene_set_run/sort_regions_by_oe.min_exp_mis_{args.min_exp_mis}{plddt_out}.gamma.ht",
        #    #_read_if_exists=True,
        #    overwrite=True,
        # )
        ht = ht.repartition(1).checkpoint(
            f"gs://gnomad/v4.1/constraint/promis3d/test_gene_set_2_run/sort_regions_by_oe.min_exp_mis_{args.min_exp_mis}{plddt_out}.gamma.ht",
            _read_if_exists=True,
            # overwrite=True,
        )
        ht.show(5)

        if args.run_forward:
            # logger.info("Running forward algorithm.")
            # res = resources.run_forward
            output_path = promis3d_res.get_forward_ht(
                name=f"oe_upper_gamma{min_exp_mis_out}{plddt_out}"
            ).path
            forward_ht = run_forward(
                ht, min_exp_mis=args.min_exp_mis, oe_upper_method="gamma"
            )
            plddt_ht = hl.read_table(
                "gs://gnomad/v2.1.1/constraint/promis3d/preprocessed_data/af2_plddt.ht"
            )
            plddt_ht = plddt_ht.annotate(plddt=hl.enumerate(plddt_ht.plddt))
            plddt_ht = plddt_ht.explode("plddt")
            plddt_ht = plddt_ht.annotate(
                residue_index=plddt_ht.plddt[0],
                plddt=plddt_ht.plddt[1],
            )
            plddt_ht = plddt_ht.key_by("uniprot_id", "residue_index")
            forward_ht = forward_ht.annotate(
                plddt=plddt_ht[forward_ht.uniprot_id, forward_ht.residue_index].plddt
            )
            # forward_ht = forward_ht.checkpoint(res.forward_ht.path, overwrite=overwrite)
            forward_ht = forward_ht.checkpoint(output_path, overwrite=overwrite)
            forward_ht.show()

        if args.run_forward_no_catch_all:
            output_path = promis3d_res.get_forward_ht(
                name=f"oe_upper_gamma.no_catch_all{min_exp_mis_out}{plddt_out}"
            ).path
            forward_ht = run_forward_no_catch_all(ht, min_exp_mis=args.min_exp_mis)
            forward_ht = forward_ht.checkpoint(output_path, overwrite=overwrite)
            forward_ht.show()

        if args.run_forward_no_catch_all_standardized:
            output_path = promis3d_res.get_forward_ht(
                name=f"oe_upper_gamma.no_catch_all_standardized{min_exp_mis_out}{plddt_out}"
            ).path
            forward_ht = run_forward_no_catch_all_standardized(
                ht, min_exp_mis=args.min_exp_mis
            )
            forward_ht = forward_ht.checkpoint(output_path, overwrite=overwrite)
            forward_ht.show()

    if args.write_per_variant:
        logger.info("Creating per-variant annotated Hail Table.")
        res = resources.write_per_variant
        res.check_resource_existence()

        # all_snv_temp_path = hl.utils.new_temp_file("all_snv", "ht")
        # ht = generate_all_possible_snvs_from_gencode_positions(
        #    res.gencode_transcipt_ht.ht(),
        #    res.gencode_translation_ht.ht().repartition(1000),
        #    res.gencode_gtf_ht.ht(),
        #    res.matched_ht.ht(),
        # ).checkpoint(all_snv_temp_path, overwrite=overwrite)
        # partition_intervals = ht._calculate_new_partitions(args.all_snv_n_partitions)
        # ht = hl.read_table(
        #    all_snv_temp_path, _intervals=partition_intervals
        # ).checkpoint(
        #    proemis3d_res.get_temp_all_possible_snvs_ht().path, overwrite=overwrite
        # )

        ht = hl.read_table(
            "gs://gnomad-tmp-4day/v4.1/constraint/proemis3d/all_possible_snvs.ht"
        )

        ht = create_per_snv_combined_ht(
            ht,
            res.forward_ht.ht(),
            res.af2_plddt_ht.ht(),
            res.af2_pae_ht.ht(),
            res.af2_dist_ht.ht(),
        )
        ht = ht.checkpoint(res.per_variant_ht.path, overwrite=overwrite)
        ht.describe()

    if args.write_per_missense_variant:
        logger.info("Filtering per-variant annotated Hail Table to missense variants.")
        res = resources.write_per_missense_variant
        res.check_resource_existence()
        ht = res.per_variant_ht.ht()
        ht = ht.filter(
            hl.any(
                ht.variant_level_annotations.transcript_consequences.map(
                    lambda x: x == "missense_variant"
                )
            )
        )
        ht = ht.checkpoint(res.per_missense_variant_ht.path, overwrite=overwrite)
        ht.show()

    if args.write_per_residue:
        logger.info("Creating per-residue annotated Hail Table.")
        res = resources.write_per_residue
        res.check_resource_existence()
        ht = create_per_residue_ht_from_snv_ht(res.per_variant_ht.ht())
        ht = ht.checkpoint(res.per_residue_ht.path, overwrite=overwrite)
        ht.show()

    if args.write_per_region:
        logger.info("Creating per-region annotated Hail Table.")
        res = resources.write_per_region
        res.check_resource_existence()
        ht = create_per_proemis3d_region_ht_from_residue_ht(res.per_residue_ht.ht())
        ht = ht.checkpoint(res.per_region_ht.path, overwrite=overwrite)
        ht.show()

    if args.create_missense_viewer_input_ht:
        logger.info("Creating missense viewer input Hail Table.")
        res = resources.create_missense_viewer_input_ht
        res.check_resource_existence()
        ht = create_missense_viewer_input_ht(
            res.gencode_pos_ht.ht(), res.forward_ht.ht()
        )
        ht = ht.checkpoint(res.missense_viewer_input_ht.path, overwrite=overwrite)
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
            f" {proemis3d_res.CURRENT_VERSION}."
        ),
        type=str,
        default=proemis3d_res.CURRENT_VERSION,
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
        "--extract-af2-plddt",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--extract-af2-pae",
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
    parser.add_argument(
        "--run-forward",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--run-forward-original",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--run-forward-no-catch-all",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--run-forward-no-catch-all-standardized",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--min-exp-mis",
        help=(
            "Minimum expected number of missense variants to consider for the greedy "
            f"algorithm. Default is {MIN_EXP_MIS}."
        ),
        type=int,
        default=MIN_EXP_MIS,
    )
    parser.add_argument(
        "--plddt-cutoff",
        help="Minimum pLDDT cutoff to filter on.",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--write-per-variant",
        action="store_true",
        help="Generate per-variant annotated HT",
    )
    parser.add_argument(
        "--all-snv-n-partitions",
        help="Number of partitions to use for the all possible SNVs Hail Table.",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--write-per-missense-variant",
        action="store_true",
        help="Generate per-variant annotated HT",
    )
    parser.add_argument(
        "--write-per-residue",
        action="store_true",
        help="Generate per-residue HT from per-variant",
    )
    parser.add_argument(
        "--write-per-region",
        action="store_true",
        help="Generate per-region HT from per-residue",
    )
    parser.add_argument(
        "--create-missense-viewer-input-ht",
        action="store_true",
        help="Create missense viewer input HT",
    )

    args = parser.parse_args()
    main(args)
