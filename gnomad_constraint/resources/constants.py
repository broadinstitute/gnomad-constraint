"""Constants used across the constraint pipeline and release formatting."""

# ---------------------------------------------------------------------------
# Pipeline configuration
# ---------------------------------------------------------------------------

EXTENSIONS = ["ht", "tsv", "tsv.bgz", "he", "log"]
"""Valid file extensions for constraint pipeline resources."""

VERSIONS = ["4.0", "4.1", "4.1.1", "4.1.2", "4.1.3"]
"""Supported gnomAD constraint pipeline versions."""

CURRENT_VERSION = "4.1.3"
"""
Default gnomAD constraint pipeline version: the next, unreleased version.

Keep it ahead of the latest release, so a run that doesn't pass ``--version`` never
writes to the paths of published tables.
"""

SITES_VERSION_MAP = {
    "4.0": "4.0",
    "4.1": "4.1",
    "4.1.1": "4.1",
    "4.1.2": "4.1",
    "4.1.3": "4.1",
}
"""Map from constraint pipeline version to gnomAD sites release version."""

VEP_VERSION_MAP = {
    "4.0": "105",
    "4.1": "105",
    "4.1.1": "105",
    "4.1.2": "105",
    "4.1.3": "105",
}
"""Map from constraint pipeline version to the VEP version of its context Table."""

GENCODE_VERSION_MAP = {"105": "v39"}
"""Map from VEP version to the GENCODE version it annotates with."""

DATA_TYPES = ["context", "exomes", "genomes"]
"""Data types used in the constraint pipeline."""

MODEL_TYPES = ["plateau", "coverage"]
"""Model types used for constraint calibration."""

VEP_ANNOTATION = "transcript_consequences"
"""
VEP annotation the models are applied to, one row per transcript consequence.

The outputs of the apply-models and metrics steps are in a subdirectory with this
name, left from when the annotation was an option.
"""

SFS_BIN_CUTOFFS = (0, 1e-6, 2e-6, 4e-6, 2e-5, 5e-5, 5e-4, 5e-3, 0.5)
"""Allele frequency upper bounds defining site frequency spectrum bins.

Variants with missing frequency are assigned bin 0. Otherwise, each variant is
assigned the index of the first cutoff its AF falls at or below.
"""

COVERAGE_CUTOFF = 90
"""
Exome AN percent at or above which a site counts as high coverage.

Low coverage sites require an extra calibration when computing the proportion of expected
variation. The value is on the default ``AN_percent`` coverage metric.
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
# Constraint percentile threshold computation and annotation
# ---------------------------------------------------------------------------

PLI_EXPECTED_VALUES = {"Null": 1.0, "Rec": 0.706, "LI": 0.207}
"""Expected o/e values for the pLI model (null, recessive, loss-of-function intolerant)."""

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
"""LoF-specific fields appended to the LoF groups in ``RELEASE_GROUPS_WITH_PLI`` in release format."""

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

RELEASE_GROUPS_WITH_RANK = ["lof_hc"]
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
    (
        "regional_rr_globals",
        "regional_rr_params",
        ["per_variant_apply_models_globals", "per_variant_row_count"],
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
# Regional mutation rate (rr)
# ---------------------------------------------------------------------------

GENOMIC_FEATURES = [
    "dist2telo",
    "dist2cent",
    "GC_content",
    "LCR",
    "SINE",
    "LINE",
    "recomb_male",
    "recomb_female",
    "met_sperm",
    "Nucleosome",
    "cDNM_maternal_05M",
    "cDNM_paternal_05M",
    "CpG_island",
]
"""
Regional genomic features used to model the regional mutation rate (rr).

Order matches the columns of the published Gnocchi feature matrix
(``genomic_features13_genome_1kb.txt``).
"""

GENOMIC_FEATURE_SCALES = ["1k", "10k", "100k", "1M"]
"""Flanking window scales at which each regional genomic feature is computed."""

GENOMIC_FEATURE_COLUMNS = [
    f"{f}_{s}" for s in GENOMIC_FEATURE_SCALES for f in GENOMIC_FEATURES
]
"""Feature matrix columns (``{feature}_{scale}``), in published column order."""

GENOMIC_FEATURE_MISSING = ["NA", "nan", "NaN", ""]
"""Spellings of a missing value in the published genomic feature matrix."""

RR_CONTEXTS = sorted(f"{b1}{c}{b3}" for b1 in "ACGT" for c in "AC" for b3 in "ACGT")
"""Strand-collapsed trinucleotide contexts (central base A or C) modeled for rr."""

CPG_CONTEXTS = ["ACG", "CCG", "GCG", "TCG"]
"""Strand-collapsed CpG trinucleotide contexts."""

METHYLATION_CORRELATED_FEATURES = [
    "GC_content",
    "SINE",
    "met_sperm",
    "Nucleosome",
    "CpG_island",
]
"""
Features correlated with methylation that are left out of the CpG context rr models.

The mutation rates already model methylation at CpG sites.
"""

RR_CONTEXT_FEATURES = {
    context: [
        f"{f}_{s}"
        for s in GENOMIC_FEATURE_SCALES
        for f in GENOMIC_FEATURES
        if not (context in CPG_CONTEXTS and f in METHYLATION_CORRELATED_FEATURES)
    ]
    for context in RR_CONTEXTS
}
"""
Features of each context's rr model.

Every feature at every scale, without the methylation-correlated features for CpG
contexts.
"""

RR_FEATURE_SELECTION_ALL = "all"
"""
``--rr-feature-selection`` value that models every feature in ``RR_CONTEXT_FEATURES``.

This is the default, and the ridge penalty then decides how much each feature counts.
"""

RR_BONFERRONI_CUTOFF = "bonferroni"
"""
``--rr-feature-selection`` value for Gnocchi's Bonferroni feature selection.

A context's feature is selected when its univariate p-value is at most
``RR_FEATURE_SELECTION_ALPHA / len(GENOMIC_FEATURE_SCALES) / n``, where ``n`` is the
number of features the context's model can use: 13, or 8 in CpG contexts, which leave
out the methylation-correlated features. Only the smallest selected scale of each
feature is kept.
"""

RR_FEATURE_SELECTION_ALPHA = 0.05
"""Family-wise alpha of the Bonferroni rr feature selection."""

RR_RIDGE_PENALTIES = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
"""
Per-site ridge penalties tried by ``expected_models.py --tune-rr-penalty``.

The step picks the penalty with the lowest held-out log loss of models fit on odd
chromosomes and scored on even ones, and vice versa, summed over contexts. In the
analysis behind the rr, this chose 0.3 from this range.

Each of those models is fit on about half the training sites, and
``--train-rr-models`` refits on all of them with the chosen penalty. The penalty is
per site (the fit minimizes the mean log loss plus the penalty term), so the final
fit is shrunk about as much as the models that were scored.
"""

RR_RIDGE_PENALTY = 0.3
"""
Ridge penalty used when ``--tune-rr-penalty`` has no held-out sites to score.

In a ``--test`` run only chr20 has features (chrX and chrY have none), so there is
no other chromosome to hold out.
"""

RR_BACKGROUND_FRACTION = 0.05
"""Fraction of non-mutated synonymous sites sampled as the rr training background."""

RR_BACKGROUND_SEED = 0
"""Seed used when sampling the non-mutated rr training background."""

RR_GENE_FLAGS_VERSION = "4.1.2"
"""Version of the constraint metrics Table whose gene flags pick the genes left out of rr training."""

RR_TRAINING_N_PARTITIONS = 100
"""
Number of partitions of the rr site features Table, and so of the training Table.

The sampled sites keep the partitioning of the preprocessed Table (about 5,000
partitions for about 500,000 rows), so the Table is coalesced before it is written.
"""
