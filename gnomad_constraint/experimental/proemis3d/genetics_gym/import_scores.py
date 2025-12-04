from typing import Callable, Dict, Any, Optional, List
import hail as hl
import argparse
import logging

import gnomad_constraint.experimental.proemis3d.genetics_gym.constants as ct
from gnomad_constraint.experimental.proemis3d.genetics_gym.constants import (
    SCORE_KEY_GROUPS,
    KEY_GROUPS,
    SCORE_FIELDS,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("genetics_gym_missense_scores_import")
logger.setLevel(logging.INFO)


def import_score(
    path: str,
    import_func: Optional[Callable[[str], hl.Table]] = None, 
    import_args: Optional[Dict[str, Any]] = None,
    parse_func: Optional[Callable[[hl.Table], hl.Table]] = None,
    parse_args: Optional[Dict[str, Any]] = None,
    rename_fields: Optional[Dict[str, str]] = None,
    key_by_fields: Optional[List[str]] = None,
    select_fields: Optional[List[str]] = None,
    naive_coalesce_partitions: Optional[int] = None,
    repartition_partitions: Optional[int] = None,
) -> hl.Table:
    """
    Import a score from a path and parse it into a Hail Table.

    .. note::

        The actions related to the parameters are performed in the order they are 
        listed in the function definition. Therefore, parse functions are run before 
        renaming fields, renaming fields are run before keying by fields, and keying by 
        fields are run before selecting fields.

    :param path: The path to the score.
    :param import_func: The function to use to import the score.
    :param import_args: The arguments to pass to the import function.
    :param parse_func: The function to use to parse the score.
    :param parse_args: The arguments to pass to the parse function.
    :param rename_fields: The fields to rename.
    :param key_by_fields: The fields to key by.
    :param select_fields: The fields to select.
    :param naive_coalesce_partitions: The number of partitions to coalesce to.
    :param repartition_partitions: The number of partitions to repartition to.
    :return: A Hail Table with the score.
    """
    if import_func is None:
        import_func = hl.import_table

    ht = import_func(path, **(import_args or {}))

    if parse_func is not None:
        ht = parse_func(ht, **(parse_args or {}))

    if rename_fields is not None:
        ht = ht.rename(rename_fields)

    if key_by_fields is not None:
        # Use new shuffle method to prevent shuffle errors.
        hl._set_flags(use_new_shuffle="1")
        ht = ht.cache().key_by(*key_by_fields)
        hl._set_flags(use_new_shuffle=None)

    if select_fields is not None:
        ht = ht.select(*select_fields)

    if repartition_partitions is not None:
        ht = ht.repartition(repartition_partitions, shuffle=True)
    elif naive_coalesce_partitions is not None:
        ht = ht.naive_coalesce(naive_coalesce_partitions)

    return ht


def parse_misfit(ht: hl.Table, mf_mapping_path: str) -> hl.Table:
    """
    Parse the Misfit score.

    `ht` input format:

        Uniprot_position    AA_alt    MisFit_D    MisFit_S       filename
        1                   A         0.678089    0.000430874    gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/MisFit_UP000005640_9606_HUMAN_S_gene.txt.gz
        1                   C         0.99952     0.00085987     gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/MisFit_UP000005640_9606_HUMAN_S_gene.txt.gz
        1                   D         0.942231    0.000760302    gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/MisFit_UP000005640_9606_HUMAN_S_gene.txt.gz

    `mf_mapping_path` input format:

        UniprotID    GeneID             TranscriptID       Symbol    ProteinID          UniprotEntryName    Strand   Chrom    Length    MisFit_sgene_mis       MisFit_sgene_ptv      logit_s_mean_mis    logit_s_sd_mis    logit_s_mean_ptv     logit_s_sd_ptv
        Q99972       ENSG00000034971    ENST00000037502    MYOC	     ENSP00000037502	MYOC_HUMAN          -        1        504       0.00581318221561312    0.0028308776672929    -5.141797           1.2674515         -5.86433362960815    2.04042434692383
        Q13395       ENSG00000059588    ENST00000040877    TARBP1	 ENSP00000040877	TARB1_HUMAN         -        1        1621      0.00281245358145288    0.0042722406797111    -5.8708816          1.5425556         -5.45133543014526    0.815543234348297
        Q15836       ENSG00000049245    ENST00000054666    VAMP3	 ENSP00000054666	VAMP3_HUMAN         +        1        100       0.014188575286456      0.0020731259137392    -4.241028           1.0843155         -6.17662239074707    1.54537832736969

    :param ht: The Hail Table to parse.
    :param mf_mapping_path: The path to the Misfit mapping table.
    :return: A Hail Table with the Misfit score.
    """
    mf_mapping_ht = hl.import_table(mf_mapping_path)
    mf_mapping_ht = mf_mapping_ht.select(
        **{
            ct.UNIPROT_ID_FIELD: mf_mapping_ht.UniprotID,
            ct.ENSEMBL_TRANSCRIPT_ID_FIELD: mf_mapping_ht.TranscriptID,
            ct.ENSEMBL_GENE_ID_FIELD: mf_mapping_ht.GeneID,
            ct.ENSEMBL_PROTEIN_ID_FIELD: mf_mapping_ht.ProteinID,
            ct.GENE_SYMBOL_FIELD: mf_mapping_ht.Symbol,
        }
    )
    mf_mapping_ht = mf_mapping_ht.key_by(ct.UNIPROT_ID_FIELD)
    
    ht = ht.transmute(**{ct.UNIPROT_ID_FIELD: ht.filename.split('/')[-1].split('\.')[0]})
    ht = ht.annotate(**mf_mapping_ht[ht[ct.UNIPROT_ID_FIELD]])

    return ht


# TODO: Confirm that numbering is all on the same base.
def parse_popeve(ht: hl.Table) -> hl.Table:
    """
    Parse the popEVE score.

    Table input format:
    
        mutant,gap frequency,popEVE,popped EVE,popped ESM-1v,EVE,ESM-1v
        M1A,0.751,-3.308,-2.907,-3.709,6.216,-3.254
        M1C,0.751,-3.473,-2.886,-4.06,6.078,-4.285
        M1D,0.751,-3.412,-2.896,-3.928,6.141,-3.826

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the popEVE score.
    """
    ht = ht.transmute(
        **{
            ct.AA_POS_FIELD: hl.int(ht.mutant[1:-1]),
            ct.AA_REF_FIELD: ht.mutant[0],
            ct.AA_ALT_FIELD: ht.mutant[-1],
        },
        NP=ht.filename.split('/')[-1][:-4],
    )

    return ht


def parse_rasp(ht: hl.Table) -> hl.Table:
    """
    Parse the RaSP score.

    `ht` input format:

        ----------------------------------------
        File Type: Table
            Partitions: 23391
            Rows: 227747600
            Empty partitions: 0
            Min(rows/partition): 283
            Max(rows/partition): 54767
            Median(rows/partition): 7443
            Mean(rows/partition): 9736
            StdDev(rows/partition): 7560
        ----------------------------------------
        Global fields:

        ----------------------------------------
        Row fields:
            'uniprot': str
            'pdbid': str
            'variant': str
            'wt_nlf': float64
            'mt_nlf': float64
            'score_ml_fermi': float64
            'score_ml': float64
            'b_factors': float64
        ----------------------------------------
        Key: ['uniprot', 'variant']
        ----------------------------------------

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the RaSP score.
    """
    ht = ht.annotate(
        **{
            ct.AA_POS_FIELD: hl.int(ht.variant[1:-1]),
            ct.AA_REF_FIELD: ht.variant[0],
            ct.AA_ALT_FIELD: ht.variant[-1],
            ct.UNIPROT_ID_FIELD: ht.uniprot,
        },
        rasp_score = ht.score_ml,
    )

    return ht


def parse_revel(ht: hl.Table) -> hl.Table:
    """
    Parse the REVEL score.

    Table input format:

        chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
        1,35142,35142,G,A,T,M,0.027,ENST00000417324
        1,35142,35142,G,C,T,R,0.035,ENST00000417324
        1,35142,35142,G,T,T,K,0.043,ENST00000417324

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the REVEL score.
    """
    ht = ht.filter(ht.grch38_pos != '.')
    ht = ht.annotate(
        locus=hl.locus(
            'chr' + ht.chr, hl.int(ht.grch38_pos), reference_genome='GRCh38'
        ),
        alleles=[ht.ref, ht.alt], 
        **{
            ct.AA_REF_FIELD: ht.aaref, 
            ct.AA_ALT_FIELD: ht.aaalt, 
            ct.ENSEMBL_TRANSCRIPT_ID_FIELD: ht.Ensembl_transcriptid.split(';'), 
        },
        revel=hl.float(ht.REVEL),
    )
    ht = ht.explode(ct.ENSEMBL_TRANSCRIPT_ID_FIELD)
    
    return ht


def parse_cpt(ht: hl.Table) -> hl.Table:
    """
    Parse the CPT score.

    Table input format:

        mutant,CPT1_score
        M1A,0.266730204583536
        M1R,0.2823860405594356
        M1N,0.4940127001804201

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the CPT score.
    """
    ht = ht.annotate(
        **{
            ct.AA_POS_FIELD: hl.int(ht.mutant[1:-1]),
            ct.AA_REF_FIELD: ht.mutant[0],
            ct.AA_ALT_FIELD: ht.mutant[-1],
        },
        cpt1_score=hl.float(ht.CPT1_score), 
        uniprot_entry=ht.filename.split('/')[-1].split('\.')[0]
    )

    return ht

# TODO: What is the difference between UNIPROT_ID_FIELD and UNIPROT_ISOFORM_FIELD and
# is it OK to use the uniprot ID to map to the isoform?
def parse_proteinmpnn(ht: hl.Table, af2_path: str) -> hl.Table:
    """
    Parse the ProteinMPNN score.

    `ht` input format:

        ----------------------------------------
        File Type: Table
            Partitions: 156
            Rows: 221282586
            Empty partitions: 0
            Min(rows/partition): 689324
            Max(rows/partition): 2689336
            Median(rows/partition): 1386583.0
            Mean(rows/partition): 1418478
            StdDev(rows/partition): 331306
        ----------------------------------------
        Global fields:

        ----------------------------------------
        Row fields:
            'AF_structure_name': str
            'uniprot': str
            'aa_pos': int32
            'aa_ref': str
            'aa_alt': str
            'ref_log_p': float64
            'alt_log_p': float64
            'proteinmpnn_llr': float64
        ----------------------------------------
        Key: ['uniprot', 'aa_pos', 'aa_ref', 'aa_alt']
        ----------------------------------------

    `af2_path` input format:

        ----------------------------------------
        File Type: Table
            Partitions: 1
            Rows: 19812
            Empty partitions: 0
        ----------------------------------------
        Global fields:

        ----------------------------------------
        Row fields:
            'uniprot_isoform': str
            'uniprot_id': str
        ----------------------------------------
        Key: ['uniprot_id']
        ----------------------------------------

    :param ht: The Hail Table to parse.
    :param af2_path: The path to the AlphaFold2 mapping of Uniprot IDs to isoforms.
    :return: A Hail Table with the ProteinMPNN score.
    """
    af2_ht = hl.read_table(af2_path)
    ht = ht.annotate(**{ct.UNIPROT_ISOFORM_FIELD: af2_ht[ht.uniprot].uniprot_isoform})

    return ht


def parse_pai3d(ht: hl.Table) -> hl.Table:
    """
    Parse the PAI3D score.

    Input format:

        chr,pos,non_flipped_ref,non_flipped_alt,gene_name,ref_aa,alt_aa,score_PAI3D
        chr1,69134,A,G,ENST00000335137.4,E,G,0.48945943988859647
        chr1,69101,A,G,ENST00000335137.4,E,G,0.6441008591651916
        chr1,69916,A,G,ENST00000335137.4,N,D,0.8164566327631475

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the PAI3D score.
    """
    ht = ht.annotate(
        locus=hl.locus(ht.chr, hl.int(ht.pos), 'GRCh38'),
        alleles=[ht.non_flipped_ref, ht.non_flipped_alt],
        **{ct.ENSEMBL_TRANSCRIPT_ID_FIELD: ht.gene_name.split('\.')[0]},
        score_PAI3D=hl.float(ht.score_PAI3D),
    )

    return ht


def parse_polyphen(ht: hl.Table) -> hl.Table:
    """
    Parse the Polyphen score.

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the Polyphen score.
    """
    csqs = ht.vep.transcript_consequences
    ht = ht.select(
        _=csqs.map(
            lambda x: x.select(
                **{
                    ct.ENSEMBL_TRANSCRIPT_ID_FIELD: x.transcript_id,
                    ct.ENSEMBL_GENE_ID_FIELD: x.gene_id,
                    ct.GENE_SYMBOL_FIELD: x.gene_symbol,
                    ct.POLYPHEN_SCORE_FIELD: x.polyphen_score,
                }
            ),
        )
    ).explode("_")
    ht = ht.annotate(**ht._)
    ht = ht.filter(
        ht[ct.ENSEMBL_TRANSCRIPT_ID_FIELD].startswith('ENST') 
        & hl.is_defined(ht.polyphen_score)
    )

    return ht


def parse_cadd(ht: hl.Table) -> hl.Table:
    """
    Parse the CADD score.

    Input format:

        ##CADD GRCh38-v1.7 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health at Charite - Universitatsmedizin Berlin 2013-2023. All rights reserved.
        #Chrom    Pos      Ref    Alt    RawScore    PHRED
        1         10001    T      A      0.767991	 7.993
        1         10001    T      C      0.799933	 8.285
        1         10001    T      G      0.777613	 8.080

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the CADD score.
    """
    ht = ht.annotate(
        locus=hl.locus(
            'chr' + ht['#Chrom'], hl.int(ht['Pos']), reference_genome='GRCh38'
        ),
        alleles=[ht['Ref'], ht['Alt']],
        cadd_score=hl.float(ht['RawScore']),
    )
    ht = ht.filter(hl.is_defined(ht.cadd_score))

    return ht


def parse_gpn_msa(ht: hl.Table) -> hl.Table:
    """
    Parse the GPN-MSA score.

    Input format:

        1    10065    C    A    -0.39
        1    10065    C    G    -2.76
        1    10065    C    T    -0.47

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the GPN-MSA score.
    """
    ht = ht.annotate(
        locus=hl.locus('chr'+ht.f0, hl.int(ht.f1), reference_genome='GRCh38'),
        alleles=[ht.f2, ht.f3], 
        gpn_msa_score=hl.float(ht.f4),
    )

    return ht


def parse_am_isoforms(ht: hl.Table) -> hl.Table:
    """
    Parse the AlphaMissense isoforms score.

    Input format:

        # Copyright 2023 DeepMind Technologies Limited
        #
        # Licensed under CC BY-NC-SA 4.0 license
        transcript_id         protein_variant    am_pathogenicity    am_class
        ENST00000000442.11    M1A                0.2808              likely_benign
        ENST00000000442.11    M1R                0.3341              likely_benign
        ENST00000000442.11    M1N                0.5455              ambiguous

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the AlphaMissense isoforms score.
    """
    ht = ht.annotate(
        **{
            ct.AA_POS_FIELD: hl.int(ht.protein_variant[1:-1]),
            ct.AA_REF_FIELD: ht.protein_variant[0],
            ct.AA_ALT_FIELD: ht.protein_variant[-1],
            ct.ENSEMBL_TRANSCRIPT_ID_FIELD: ht.transcript_id.split('\.')[0],
        },
        AM_score=hl.float(ht.am_pathogenicity), 
        enst_orig=ht.transcript_id,
    )

    return ht


def parse_am_canonical(ht: hl.Table) -> hl.Table:
    """
    Parse the AlphaMissense canonical score.

    Input format:
    
        # Copyright 2023 DeepMind Technologies Limited
        #
        # Licensed under CC BY-NC-SA 4.0 license
        #CHROM    POS      REF    ALT    genome    uniprot_id    transcript_id        protein_variant    am_pathogenicity    am_class
        chr1      69094    G      T      hg38      Q8NH21        ENST00000335137.4    V2L                0.2937              likely_benign
        chr1      69094    G      C      hg38      Q8NH21        ENST00000335137.4    V2L                0.2937              likely_benign
        chr1      69094    G      A      hg38      Q8NH21        ENST00000335137.4    V2M                0.3296              likely_benign

    :param ht: The Hail Table to parse.
    :return: A Hail Table with the AlphaMissense canonical score.
    """
    ht = ht.annotate(
        locus=hl.locus(ht.f0, hl.int(ht.f1), reference_genome='GRCh38'),
        alleles=[ht.f2, ht.f3], 
        **{
            ct.UNIPROT_ID_FIELD: ht.f5,
            ct.ENSEMBL_TRANSCRIPT_ID_FIELD: ht.f6.split('\.')[0],
        },
        AM_score=hl.float(ht.f8), 
        enst_orig=ht.f6, 
        consequence=ht.f9, 
        genome=ht.f4,
    )

    return ht


def get_scores_config() -> Dict[str, Dict[str, Any]]:
    """
    Get the scores configuration.

    :return: A dictionary of scores configuration.
    """
    return {
        "esm1b": {
            "path": ct.ESM1B_PATH,
            "import_args": {
                "types": {
                    "aa_pos": hl.tint,
                    "aa_ref": hl.tstr,
                    "aa_alt": hl.tstr,
                    "esm_score": hl.tfloat,
                    "uniprot": hl.tstr,
                },
                "force": True,
                "missing": "",
            },
            "rename_fields": {
                "uniprot": ct.UNIPROT_ID_FIELD, "esm_score": ct.ESM1B_SCORE_FIELD
            },
        },
        "rasp": {
            "path": ct.RASP_PATH,
            "import_func": hl.read_table,
            "parse_func": parse_rasp,
        },
        "proteinmpnn": {
            "path": ct.PROTEINMPNN_PATH,
            "import_func": hl.read_table,
            "parse_func": parse_proteinmpnn,
            "parse_args": {"af2_path": ct.AF2_UNIPROT_ISOFORM_MAPPING_PATH},
            "rename_fields": {"uniprot": ct.UNIPROT_ID_FIELD},
            "select_fields": [ct.UNIPROT_ID_FIELD, *SCORE_FIELDS["proteinmpnn"]],
            "repartition_partitions": 1000,
        }, 
        "misfit": {
            "path": ct.MISFIT_PATH, 
            "import_args": {
                "source_file_field": "filename",
                "types": {
                    "Uniprot_position": hl.tint,
                    "AA_alt": hl.tstr,
                    "MisFit_D": hl.tfloat,
                    "MisFit_S": hl.tfloat,
                },
                "force": True,
                "missing": "",
            },
            "parse_func": parse_misfit,
            "parse_args": {"mf_mapping_path": ct.MISFIT_MAPPING_PATH},
            "rename_fields": {"Uniprot_position": ct.AA_POS_FIELD, "AA_alt": ct.AA_ALT_FIELD},
        },
        "am_isos": {
            "path": ct.AM_ISOFORMS_PATH, 
            "import_args": {
                "force_bgz": True, 
                "no_header": False, 
                "comment": '#',
                "min_partitions": 1000,
            }, 
            "parse_func": parse_am_isoforms,
            "select_fields": [*SCORE_FIELDS["am_isos"], 'enst_orig'],
        },
        "popeve": {
            "path": ct.POPEVE_PATH,
            "import_args": {
                "delimiter": ",",
                "source_file_field": "filename",
                "types": {
                    "popEVE": hl.tfloat,
                    "popped EVE": hl.tfloat,
                    "popped ESM-1v": hl.tfloat,
                    "EVE": hl.tfloat,
                    "ESM-1v": hl.tfloat,
                },
                "force": True,
                "missing": "",
            },
            "parse_func": parse_popeve,
            "rename_fields": {
                "popped EVE": "popped_EVE", 
                "popped ESM-1v": "popped_ESM_1v", 
                "ESM-1v": ct.ESM_1V_SCORE_FIELD,
            },
            "select_fields": [
                *SCORE_FIELDS["popeve"],
                "NP",
                "popped_EVE",
                "popped_ESM_1v",
            ],
        },
        "cpt": {
            "path": ct.CPT_PATH, 
            "import_args": {
                "delimiter": ",", 
                "force": True, 
                "source_file_field": "filename",
            },
            "parse_func": parse_cpt,
        },
        "am_canonical": {
            "path": ct.AM_HG38_PATH, 
            "import_args": {
                "force_bgz": True, 
                "no_header": True, 
                "comment": '#',
                "min_partitions": 1000,
            }, 
            "parse_func": parse_am_canonical,
            "select_fields": [
                *SCORE_FIELDS["am_canonical"], 'enst_orig', 'consequence', 'genome'
            ],
        },
        "pai3d": {
            "path": ct.PRIMATEAI3D_PATH, 
            "import_func": hl.import_csv,
            "import_args": {"min_partitions": 1000},
            "parse_func": parse_pai3d,
            "rename_fields": {"ref_aa": ct.AA_REF_FIELD, "alt_aa": ct.AA_ALT_FIELD},
            "select_fields":[
                ct.AA_REF_FIELD, ct.AA_ALT_FIELD, *SCORE_FIELDS["pai3d"]
            ],
        },
        "revel": {
            "path": ct.REVEL_PATH, 
            "import_args": {"delimiter": ',', "min_partitions": 1000},
            "parse_func": parse_revel,
        },
        "polyphen": {
            "path": ct.CONTEXT_VEP_ANNOTATED_PATH,
            "import_func": hl.read_table,
            "parse_func": parse_polyphen,
            "select_fields": [
                ct.ENSEMBL_GENE_ID_FIELD, 
                ct.GENE_SYMBOL_FIELD, 
                *SCORE_FIELDS["polyphen"],
            ],
        },
        "cadd": {
            "path": ct.CADD_PATH, 
            "import_args": {"filter": '##', "force_bgz": True, "min_partitions": 1000},
            "parse_func": parse_cadd,
        },
        "gpn_msa": {
            "path": ct.GPN_MSA_PATH, 
            "import_args": {
                "force_bgz": True, "no_header": True, "min_partitions": 1000
            }, 
            "parse_func": parse_gpn_msa,
        },
        # TODO: Not tested because I don't have access to this bucket, but why not use
        # The public resource and lift over to GRCh38?
        #"mpc": {
        #    "path": ct.MPC_PATH,
        #    "import_func": hl.read_table,
        #},
    }


def main(args):
    """Import missense scores for Genetics Gym."""
    hl.init(
        log="/genetics_gym_missense_scores_import.log",
        #tmp_dir="gs://trisha-tmp",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    overwrite = args.overwrite
    resources = args.resources

    scores_config = get_scores_config()

    if args.import_all:
        resources = scores_config.keys()

    for r in resources:
        resource_info = scores_config[r]
        logger.info(f"Importing {r}...")
        ht = import_score(
            resource_info["path"], 
            import_func=resource_info.get("import_func"), 
            import_args=resource_info.get("import_args"), 
            parse_func=resource_info.get("parse_func"), 
            parse_args=resource_info.get("parse_args"), 
            rename_fields=resource_info.get("rename_fields"), 
            key_by_fields=resource_info.get("key_by_fields", KEY_GROUPS[SCORE_KEY_GROUPS[r]]), 
            select_fields=resource_info.get("select_fields", SCORE_FIELDS[r]), 
            naive_coalesce_partitions=resource_info.get("naive_coalesce_partitions", 1000),
            repartition_partitions=resource_info.get("repartition_partitions"),
        )
        ht.write(
            #f"gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/{r}.ht", 
            f"gs://gnomad-julia/genetics_gym/imported_vsm/{r}.ht", 
            overwrite=overwrite,
        )
        logger.info(f"Finished importing {r}.\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--import-all",
        help="Import all data.",
        action="store_true",
    )
    parser.add_argument(
        "--resources",
        choices=list(get_scores_config().keys()),
        nargs="+",
        help="Resource to import. Choices are:\n\n"
        + "\n".join([f"- {key}" for key in get_scores_config().keys()]),
        default=None,
    )

    args = parser.parse_args()
    main(args)
