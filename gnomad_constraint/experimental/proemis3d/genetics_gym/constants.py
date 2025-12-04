ENSEMBL_TRANSCRIPT_ID_FIELD = "enst"
"""The Ensembl transcript ID field."""
ENSEMBL_GENE_ID_FIELD = "ensg"
"""The Ensembl gene ID field."""
ENSEMBL_PROTEIN_ID_FIELD = "ensp"
"""The Ensembl protein ID field."""
GENE_SYMBOL_FIELD = "gene_symbol"
"""The gene symbol field."""

AA_POS_FIELD = "aa_pos"
"""The amino acid position field."""
AA_REF_FIELD = "aa_ref"
"""The amino acid reference field."""
AA_ALT_FIELD = "aa_alt"
"""The amino acid alternative field."""

UNIPROT_ID_FIELD = "uniprot_id"
"""The UniProt ID field."""
UNIPROT_ISOFORM_FIELD = "uniprot_isoform"
"""The UniProt isoform field."""

VARIANT_KEY_FIELDS = ["locus", "alleles"]
"""The variant key fields."""
VARIANT_TRANSCRIPT_KEY_FIELDS = [*VARIANT_KEY_FIELDS, ENSEMBL_TRANSCRIPT_ID_FIELD]
"""The variant transcript key fields."""
AA_SUBSTITUTION_KEY_FIELDS = [AA_POS_FIELD, AA_REF_FIELD, AA_ALT_FIELD]
"""The amino acid substitution key fields."""
AA_SUBSTITUTION_UNIPROT_KEY_FIELDS = [*AA_SUBSTITUTION_KEY_FIELDS, UNIPROT_ID_FIELD]
"""The amino acid substitution UniProt key fields."""
AA_SUBSTITUTION_TRANSCRIPT_KEY_FIELDS = [
    *AA_SUBSTITUTION_KEY_FIELDS, ENSEMBL_TRANSCRIPT_ID_FIELD
]
"""The amino acid substitution UniProt transcript key fields."""
VARIANT_UNIPROT_TRANSCRIPT_KEY_FIELDS = [
    *VARIANT_KEY_FIELDS, UNIPROT_ID_FIELD, ENSEMBL_TRANSCRIPT_ID_FIELD
]
"""The variant UniProt transcript key fields."""
AA_ALT_UNIPROT_TRANSCRIPT_KEY_FIELDS = [
    AA_POS_FIELD, 
    AA_ALT_FIELD, 
    UNIPROT_ID_FIELD,
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
LINKER_KEY_FIELDS = [
    *VARIANT_KEY_FIELDS, 
    *AA_SUBSTITUTION_UNIPROT_KEY_FIELDS, 
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
"""The linker key fields."""

KEY_GROUPS = {
    "variant": VARIANT_KEY_FIELDS,
    "variant_transcript": VARIANT_TRANSCRIPT_KEY_FIELDS,
    "variant_uniprot_transcript": VARIANT_UNIPROT_TRANSCRIPT_KEY_FIELDS,
    "aa_substitution": AA_SUBSTITUTION_KEY_FIELDS,
    "aa_substitution_transcript": AA_SUBSTITUTION_TRANSCRIPT_KEY_FIELDS,
    "aa_substitution_uniprot": AA_SUBSTITUTION_UNIPROT_KEY_FIELDS,
    "aa_substitution_alt_uniprot_transcript": AA_ALT_UNIPROT_TRANSCRIPT_KEY_FIELDS,
}
"""The table key groups."""

ESM1B_SCORE_FIELD = "esm1b"
MISFIT_D_SCORE_FIELD = "MisFit_D"
MISFIT_S_SCORE_FIELD = "MisFit_S"
POPEVE_SCORE_FIELD = "popEVE"
EVE_SCORE_FIELD = "EVE"
ESM_1V_SCORE_FIELD = "ESM_1v"
MPC_SCORE_FIELD = "mpc"
RASP_SCORE_FIELD = "rasp_score"
REVEL_SCORE_FIELD = "revel"
CPT1_SCORE_FIELD = "cpt1_score"
PROTEINMPNN_LLR_SCORE_FIELD = "proteinmpnn_llr"
PAI3D_SCORE_FIELD = "score_PAI3D"
POLYPHEN_SCORE_FIELD = "polyphen_score"
CADD_SCORE_FIELD = "cadd_score"
GPN_MSA_SCORE_FIELD = "gpn_msa_score"
AM_SCORE_FIELD = "AM_score"
AM_CANONICAL_SCORE_FIELD = "AM_score"

SCORE_FIELDS = {
    "esm1b": [ESM1B_SCORE_FIELD],
    "misfit": [MISFIT_D_SCORE_FIELD, MISFIT_S_SCORE_FIELD],
    "popeve": [POPEVE_SCORE_FIELD, EVE_SCORE_FIELD, ESM_1V_SCORE_FIELD],
    "mpc": [MPC_SCORE_FIELD],
    "rasp": [RASP_SCORE_FIELD], 
    "revel": [REVEL_SCORE_FIELD], 
    "cpt": [CPT1_SCORE_FIELD],
    "proteinmpnn": [PROTEINMPNN_LLR_SCORE_FIELD], 
    "pai3d": [PAI3D_SCORE_FIELD], 
    "polyphen": [POLYPHEN_SCORE_FIELD],
    "cadd": [CADD_SCORE_FIELD],
    "gpn_msa": [GPN_MSA_SCORE_FIELD],
    "am_isos": [AM_SCORE_FIELD],
    "am_canonical": [AM_SCORE_FIELD],
}
"""The score fields in each score resource."""

SCORE_KEY_GROUPS = {
    "esm1b": "aa_substitution_uniprot",
    "misfit": "aa_substitution_alt_uniprot_transcript",
    "popeve": "aa_substitution",
    #"mpc": "variant", # TODO: Uncomment this when MPC scores are processed.
    "rasp": "aa_substitution_uniprot", 
    "revel": "variant_transcript", 
    "cpt": "aa_substitution",
    "proteinmpnn": "aa_substitution_uniprot", # *ct.AA_SUBSTITUTION_KEY_FIELDS, ct.UNIPROT_ISOFORM_FIELD
    "pai3d": "variant_transcript", 
    "polyphen": "variant_transcript",
    "cadd": "variant",
    "gpn_msa": "variant",
    "am_isos": "aa_substitution_transcript",
    "am_canonical": "variant_uniprot_transcript",
}
    
HIGHER_IS_LESS_DELETERIOUS = {
    "esm1b": {ESM1B_SCORE_FIELD: True},
    "misfit":{MISFIT_D_SCORE_FIELD: False, MISFIT_S_SCORE_FIELD: False},
    "popeve":{
        POPEVE_SCORE_FIELD: True, 
        EVE_SCORE_FIELD: False, 
        ESM_1V_SCORE_FIELD: True,
    },
    "mpc": {MPC_SCORE_FIELD: False},
    "rasp": {RASP_SCORE_FIELD: False},
    "revel": {REVEL_SCORE_FIELD: False},
    "cpt": {CPT1_SCORE_FIELD: False},
    "proteinmpnn": {PROTEINMPNN_LLR_SCORE_FIELD: True},
    "pai3d": {PAI3D_SCORE_FIELD: False},
    "polyphen": {POLYPHEN_SCORE_FIELD: False},
    "cadd": {CADD_SCORE_FIELD: False},
    "gpn_msa": {GPN_MSA_SCORE_FIELD: False},
    "am_isos": {AM_SCORE_FIELD: False},
    "am_canonical": {AM_SCORE_FIELD: False},
}
"""Whether higher scores are less deleterious for each score."""

ESM1B_PATH = "gs://missense-scoring/esm1b-650m-brandes/proc/*.txt.bgz"
MISFIT_PATH = "gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/*.txt.gz"
MISFIT_MAPPING_PATH = "gs://missense-scoring/misfit/raw/geneset_s_gene.txt"
POPEVE_PATH = "gs://missense-scoring/popEVE_ukbb_20250312/*.csv"
MPC_PATH = "gs://asc-v17/mpc_gnomad_2.1.1/mpc_grch38_deduped_with_outliers_2024-04-30.ht"
RASP_PATH = "gs://nnfc-fdp-konrad-public/RaSP/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.ht"
REVEL_PATH = "gs://missense-scoring/revel_with_transcript_ids"
CPT_PATH = "gs://missense-scoring/cpt_all/*.csv.gz"
PROTEINMPNN_PATH = "gs://missense-scoring/mpnn-outputs/AF_total_variants.ht"
AF2_UNIPROT_ISOFORM_MAPPING_PATH = "gs://genetics-gym-not-public/Emily/af2db_uniprot_final.ht"
PRIMATEAI3D_PATH = "gs://missense-scoring/primate_ai3d/PrimateAI-3D_scores.csv.bgz"
CONTEXT_VEP_ANNOTATED_PATH = "gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht"
CADD_PATH = "gs://genetics-gym-not-public/Trisha/whole_genome_SNVs.tsv.gz"
GPN_MSA_PATH = "gs://missense-scoring/GPN-MSA/scores.tsv.bgz"
AM_ISOFORMS_PATH = "gs://dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz"
AM_HG38_PATH = "gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz"
