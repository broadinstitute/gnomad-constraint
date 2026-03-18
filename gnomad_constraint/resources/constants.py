"""Constants used across the constraint pipeline and release formatting."""

# ---------------------------------------------------------------------------
# Pipeline configuration
# ---------------------------------------------------------------------------

EXTENSIONS = ["ht", "tsv", "tsv.bgz", "he", "log"]
"""Valid file extensions for constraint pipeline resources."""

VERSIONS = ["2.1.1", "4.0", "4.1", "4.1.1"]
"""Supported gnomAD constraint pipeline versions."""

CURRENT_VERSION = "4.1.1"
"""Current default gnomAD constraint pipeline version."""

DATA_TYPES = ["context", "exomes", "genomes"]
"""Data types used in the constraint pipeline."""

MODEL_TYPES = ["plateau", "coverage"]
"""Model types used for constraint calibration."""

GENOMIC_REGIONS = ["autosome_par", "chrx_nonpar", "chry_nonpar"]
"""Genomic regions used to partition constraint calculations."""

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

SFS_BIN_CUTOFFS = (0, 1e-6, 2e-6, 4e-6, 2e-5, 5e-5, 5e-4, 5e-3, 0.5)
"""Allele frequency upper bounds defining site frequency spectrum bins.

Variants with missing frequency are assigned bin 0. Otherwise, each variant is
assigned the index of the first cutoff its AF falls at or below.
"""

COVERAGE_CUTOFF = 40
"""
Minimum median exome coverage differentiating high coverage sites from low coverage sites.

Low coverage sites require an extra calibration when computing the proportion of expected variation.
"""

CLASSIC_LOF_ANNOTATIONS = (
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
)
"""Classic loss-of-function VEP annotations."""

MU_GROUPING = ("context", "ref", "alt", "methylation_level")
"""
Annotations used to group variants for the mutation rate calculation.
"""

CALIBRATION_GROUPING = ("genomic_region", "build_model", "cpg", "exomes_coverage")
"""
Annotations used to group variants for the mutation rate calibration.
"""

AGGREGATE_SUM_FIELDS = (
    "mu_snp",
    "mu",
    "observed_variants",
    "possible_variants",
    "predicted_proportion_observed",
    "coverage_correction",
    "expected_variants",
)
"""
Fields to sum (or array sum) when aggregating the expected counts Table.
"""

MUTATION_TYPE_FIELDS = (
    "cpg",
    "transition",
    "mutation_type",
    "mutation_type_model",
)
"""
Fields added by `annotate_mutation_type`.
"""

# ---------------------------------------------------------------------------
# Frequency metadata
# ---------------------------------------------------------------------------

ADJ_FREQ_META = {"group": "adj"}
"""Frequency metadata key for the adjusted allele frequency group."""

# ---------------------------------------------------------------------------
# GENCODE field renames
# ---------------------------------------------------------------------------

GENCODE_FIELD_RENAMES = {
    "transcript_id_version": "transcript_version",
    "level": "transcript_level",
}
"""GENCODE field renames applied when preparing the release Table."""

# ---------------------------------------------------------------------------
# Release format constants
# ---------------------------------------------------------------------------

RELEASE_KEY_ORDER = ["gene", "gene_id", "transcript", "canonical", "mane_select"]
"""Key fields for the release Table, in display order."""

RELEASE_SCALAR_FIELDS = [
    "transcript_version",
    "transcript_type",
    "transcript_level",
    "chromosome",
    "start_position",
    "end_position",
    "cds_length",
    "num_coding_exons",
]
"""Scalar transcript/gene annotation fields included in the flat release Table."""

RELEASE_TOP_LEVEL_ANNOTATIONS = RELEASE_SCALAR_FIELDS + [
    "gene_quality_metrics",
    "gene_flags",
    "constraint_flags",
]
"""Non-key, non-constraint-group row fields included in the release Table."""

RELEASE_GROUP_NAMES = ["syn", "mis", "lof_hc_lc", "lof"]
"""Constraint group names exposed in the release Table, in display order."""

RELEASE_LOF_FIELDS = ["pLI", "pNull", "pRec"]
"""LoF-specific fields appended to the ``lof`` constraint group in release format."""

RELEASE_CI_FIELDS = ["lower", "upper"]
"""OE confidence interval sub-fields included in the release Table."""

RELEASE_RANK_FIELDS = [
    "upper_rank",
    "upper_bin_percentile",
    "upper_bin_decile",
    "upper_bin_sextile",
]
"""Rank and bin sub-fields appended to the CI struct for ranked groups."""

RELEASE_CI_FIELDS_WITH_RANK = RELEASE_CI_FIELDS + RELEASE_RANK_FIELDS
"""CI fields including rank annotations, used for groups in RELEASE_GROUPS_WITH_RANK."""

RELEASE_GROUPS_WITH_RANK = ["lof_hc", "lof_hc_lc"]
"""Constraint groups for which rank and bin annotations are included."""

RELEASE_GROUPS_WITH_PLI = ["lof_hc", "lof_hc_lc"]
"""Constraint groups that include pLI/pNull/pRec in release format."""

RELEASE_GROUP_RENAMES = {"lof_hc": "lof"}
"""Internal constraint group names that are renamed for public release."""

RELEASE_PIPELINE_PARAM_GLOBALS = [
    (
        "calculate_mu_globals",
        "calculate_mu_params",
        ["freq_meta", "genetic_ancestry_groups", "downsampling_idx"],
    ),
    (
        "build_models_globals",
        "build_models_params",
        ["synonymous_transcript_filter_field", "skip_coverage_model"],
    ),
    (
        "apply_models_globals",
        "apply_models_params",
        ["skip_coverage_model", "groupings"],
    ),
]
"""Pipeline parameter globals: (internal name, release name, fields to drop)."""

RELEASE_CG_RENAME = {
    "mu_snp": "mu",
    "possible_variants": "possible",
    "observed_variants": "obs",
    "expected_variants": "exp",
}
"""Field renames applied to the constraint-group structs in release format."""

RELEASE_CG_SELECT = [
    "mu",
    "possible",
    "obs",
    "exp",
    "oe",
    "z_raw",
    "z_score",
    "oe_ci",
    "gen_anc_obs",
    "gen_anc_exp",
]
"""Fields selected from the release constraint-group struct, in display order."""

# ---------------------------------------------------------------------------
# Constraint percentile threshold computation and annotation
# ---------------------------------------------------------------------------

CONSTRAINT_SCORE_CAP = 2.0
"""OE upper CI value above which scores are considered unconstrained and capped."""

PLI_EXPECTED_VALUES = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}
"""Expected o/e values for the pLI model (null, recessive, loss-of-function intolerant)."""

CONSTRAINT_METRICS = ["lof", "mis", "syn"]
"""Constraint metrics for which percentile thresholds are computed."""

CONSTRAINT_GRANULARITIES = {
    "percentile": list(range(1, 100)),
    "decile": list(range(1, 10)),
    "sextile": list(range(1, 6)),
}
"""
Granularities for percentile binning.

Keys are granularity names; values are boundary bin labels (1-indexed).
Quantile probabilities are bin / (max_bin + 1).
"""
