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

SITES_VERSION_MAP = {
    "2.1.1": "2.1.1",
    "4.0": "4.0",
    "4.1": "4.1",
    "4.1.1": "4.1",
}
"""Map from constraint pipeline version to gnomAD sites release version."""

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
# Constraint percentile threshold computation and annotation
# ---------------------------------------------------------------------------

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
        [],
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
    "unadj",
    "gen_anc_obs",
    "gen_anc_exp",
]
"""Fields selected from the release constraint-group struct, in display order."""

RELEASE_UNADJ_FIELDS = ["exp", "oe", "oe_ci", "z_raw", "z_score"]
"""
Fields in the release ``unadj`` struct (metrics without the regional rr adjustment).
"""

# ---------------------------------------------------------------------------
# Regional mutation rate (rr) adjustment
# ---------------------------------------------------------------------------

UNADJ_EXPECTED_FIELD = "unadj_expected"
"""
Per-variant and aggregated field holding the expected variant count before the
regional rr adjustment is applied.
"""

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

RR_SCREEN_FEATURE_ORDER = [
    "dist2telo",
    "dist2cent",
    "LCR",
    "SINE",
    "LINE",
    "GC_content",
    "recomb_male",
    "recomb_female",
    "met_sperm",
    "Nucleosome",
    "CpG_island",
    "cDNM_maternal_05M",
    "cDNM_paternal_05M",
]
"""
Order in which features are screened for rr feature selection.

This is the legacy Gnocchi screening order (it differs from the published column
order), which determines the column order of the rr models.
"""

GENOMIC_FEATURE_SCALES = ["1k", "10k", "100k", "1M"]
"""Flanking window scales at which each regional genomic feature is computed."""

GENOMIC_FEATURE_COLUMNS = [
    f"{f}_{s}" for s in GENOMIC_FEATURE_SCALES for f in GENOMIC_FEATURES
]
"""Feature matrix columns (``{feature}_{scale}``), in published column order."""

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
Features correlated with methylation that are dropped from CpG context rr models.
"""

RR_FEATURE_SELECTION_ALPHA = 0.05
"""Family-wise alpha used by the legacy Bonferroni rr feature selection."""

RR_BONFERRONI_CUTOFF = "bonferroni"
"""
Feature selection cutoff token for the legacy Bonferroni rule.

A feature is selected when its univariate p-value is at most
``RR_FEATURE_SELECTION_ALPHA / len(GENOMIC_FEATURE_SCALES) / n``, where ``n`` is 13
for non-CpG contexts and 8 (13 minus the methylation-correlated features) for CpG
contexts.
"""

RR_BACKGROUND_FRACTION = 0.05
"""Fraction of non-mutated synonymous sites sampled as the rr training background."""

RR_BACKGROUND_SEED = 0
"""Seed used when sampling the non-mutated rr training background."""

RR_CLIP_BOUNDS = (0.2, 5.0)
"""Default lower and upper bounds used to clip rr when safeguards are applied."""

RR_XY_MASKED_FEATURES = [
    "recomb_male",
    "recomb_female",
    "met_sperm",
    "Nucleosome",
    "cDNM_maternal_05M",
    "cDNM_paternal_05M",
]
"""
Features masked (pinned to the genome-average window) on chrX and chrY windows.

These features have no source data covering chrX/chrY (or the source track is not
available), so their values are missing in the X/Y feature matrix.
"""

GENOMIC_FEATURE_COVERAGE_TRACKS = {
    "LCR": "btu356_LCR-hs38.bed",
    "SINE": "repeatMasker_hg38.SINE.bed",
    "LINE": "repeatMasker_hg38.LINE.bed",
    "CpG_island": "cpgIslandExt.bed",
}
"""
Source tracks (in the genomic feature track directory) for coverage features.

Each feature is the percent of the window covered by the (merged) track intervals.
"""

GENOMIC_FEATURE_DISTANCE_TRACKS = {
    "dist2telo": "telomeres_hg38.bed",
    "dist2cent": "centromeres_hg38.bed",
}
"""
Source tracks (in the genomic feature track directory) for distance features.

Each feature is the distance from the window midpoint to the nearest track
interval, as a percent of the chromosome length.
"""

GENOMIC_FEATURE_WINDOW_BEDS = {
    "1k": "hg38.chrom.1kb.bed",
    "10k": "hg38.chrom.1kb_flnk_10k.bed",
    "100k": "hg38.chrom.1kb_flnk_100k.bed",
    "1M": "hg38.chrom.1kb_flnk_1M.bed",
}
"""Window (and flanking window) BEDs, per feature scale, keyed by 1kb element ID."""

GENOMIC_FEATURE_MID_BEDS = {
    "1k": "hg38.chrom.1kb.mid.bed",
    "10k": "hg38.chrom.1kb_flnk_10k.mid.bed",
    "100k": "hg38.chrom.1kb_flnk_100k.mid.bed",
    "1M": "hg38.chrom.1kb_flnk_1M.mid.bed",
}
"""Window midpoint BEDs, per feature scale, used for the distance features."""

RR_STATSMODELS_VERSION = "0.14.4"
"""
statsmodels version used to fit the legacy syn_adj_r rr models.

``Logit.fit_regularized`` results differ slightly (~1e-5 relative rr) with
statsmodels >= 0.14.5, so this version (which requires scipy < 1.16) is needed to
exactly reproduce the legacy models.
"""

RR_FOLDS = ("odd", "even")
"""
Chromosome-split folds used for out-of-fold rr evaluation.

Autosomes are assigned by chromosome number parity, chrX to "odd", and chrY to
"even".
"""
