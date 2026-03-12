"""
Constraint pipeline for DD de novo variants — obs/exp only.

Computes observed and expected variant counts for DD case de novo variants
using the gnomAD constraint framework. Stops after per-transcript obs/exp
aggregation (does not compute pLI, z-scores, or confidence intervals).

Uses the v4.1 coverage-corrected preprocessed context table, which already
has all constraint annotations (trimer context, mutation type, methylation,
coverage, VEP). No additional preprocessing is needed.

Uses imports from:
  - gnomad.utils.constraint (jg/add_functions_for_per_base_constraint branch)
"""

import logging

import hail as hl
from bokeh.io import output_file, save
from bokeh.models import TabPanel, Tabs

# --- gnomad_methods (jg/add_functions_for_per_base_constraint branch) ---
from gnomad.utils.constraint import (
    annotate_exploded_vep_for_constraint_groupings,
    annotate_with_mu,
    apply_plateau_models,
    build_plateau_models,
    count_observed_and_possible_by_group,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_test_dd")
logger.setLevel(logging.INFO)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HIGH_COVERAGE_CUTOFF = 90

# GCS paths (coverage-corrected v4.1 resources)
DD_HT_PATH = "gs://gnomad-julia/constraint/dd.ht"
PREPROCESSED_CONTEXT_PATH = (
    "gs://gnomad/v4.1/constraint_coverage_corrected/"
    "preprocessed_data/gnomad.v4.1.context.preprocessed.ht"
)
MUTATION_RATE_PATH = (
    "gs://gnomad/v4.1/constraint_coverage_corrected/"
    "mutation_rate/gnomad.v4.1.mutation_rate.ht"
)
CHECKPOINT_ROOT = "gs://gnomad-tmp/julia/constraint"


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    """Run DD de novo constraint analysis (obs/exp only)."""
    hl.init(tmp_dir="gs://gnomad-tmp-4day")

    # -------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------
    logger.info("Loading preprocessed context, mutation rate, and DD variant tables...")
    context_ht = hl.read_table(PREPROCESSED_CONTEXT_PATH)
    mutation_ht = hl.read_table(MUTATION_RATE_PATH)
    dd_ht = hl.read_table(DD_HT_PATH)

    dd_ht.show()
    print(f"DD positive variants: {dd_ht.filter(dd_ht.is_pos_dd).count()}")

    # -------------------------------------------------------------------
    # Annotate DD status onto preprocessed context
    # -------------------------------------------------------------------
    logger.info("Joining DD variant status onto preprocessed context...")
    context_ht = context_ht.annotate(**dd_ht[context_ht.key])
    context_ht.filter(hl.is_defined(context_ht.is_pos_dd)).show()

    # -------------------------------------------------------------------
    # Explode VEP and checkpoint
    # -------------------------------------------------------------------
    logger.info("Exploding VEP transcript consequences...")
    ht = context_ht.filter(
        hl.is_defined(context_ht.exomes_coverage)
        & hl.or_else(hl.len(context_ht.filters.exomes) == 0, True)
    )

    ht, grouping = annotate_exploded_vep_for_constraint_groupings(
        ht,
        vep_annotation="transcript_consequences",
        include_canonical_group=True,
        include_mane_select_group=True,
    )
    ht = ht.checkpoint(
        f"{CHECKPOINT_ROOT}/dd_exploded_vep.gnomad_pos.v2.ht",
        _read_if_exists=True,
    )
    print("Grouping:", grouping)
    ht.show()

    # -------------------------------------------------------------------
    # Count observed and possible variants
    # -------------------------------------------------------------------
    logger.info("Counting observed and possible variants...")
    hl._set_flags(use_new_shuffle="1")

    ht = hl.read_table(f"{CHECKPOINT_ROOT}/dd_exploded_vep.gnomad_pos.v2.ht")

    # DD variant is "observed" if is_pos_dd is True; all sites are "possible".
    observed_expr = hl.array([hl.int(hl.or_else(ht.is_pos_dd, False))])
    possible_expr = hl.int(1)

    keys = (
        "context", "ref", "alt", "methylation_level", "exomes_coverage",
    ) + grouping
    additional_grouping = list(keys) + ["cpg", "mutation_type"]

    ht = count_observed_and_possible_by_group(
        ht,
        possible_expr=possible_expr,
        observed_expr=observed_expr,
        additional_grouping=additional_grouping,
        partition_hint=1000,
    )
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/julia/constraint/"
        "dd_observed_possible.gnomad_pos.intermediate.exomes_coverage.v2.ht",
        _read_if_exists=True,
    )
    ht = ht.key_by(*keys)

    # Annotate with mutation rate.
    ht = annotate_with_mu(ht, mutation_ht)

    ht = ht.checkpoint(
        f"{CHECKPOINT_ROOT}/dd_observed_possible.gnomad_pos.exomes_coverage.v2.ht",
        _read_if_exists=True,
    )
    ht.show()

    hl._set_flags(use_new_shuffle=None)

    # -------------------------------------------------------------------
    # Build plateau models (all sites + high coverage)
    # -------------------------------------------------------------------
    logger.info("Building plateau models...")
    ht = hl.read_table(
        f"{CHECKPOINT_ROOT}/dd_observed_possible.gnomad_pos.exomes_coverage.v2.ht"
    )

    # Filter to synonymous MANE Select Ensembl transcripts for training.
    ht_syn = ht.filter(
        ht.mane_select
        & (ht.annotation == "synonymous_variant")
        & ht.transcript.startswith("ENST")
    )
    high_cov_ht = ht_syn.filter(ht_syn.exomes_coverage >= HIGH_COVERAGE_CUTOFF)

    # Group by mutation rate keys for model building.
    model_keys = ("context", "ref", "alt", "methylation_level", "mu_snp")
    mu_type_fields = ("cpg", "mutation_type")

    def aggregate_for_models(t: hl.Table) -> hl.Table:
        """Aggregate training data by mutation rate keys."""
        return (
            t.group_by(*(model_keys + mu_type_fields))
            .aggregate(
                observed_variants=hl.agg.array_sum(t.observed_variants),
                possible_variants=hl.agg.sum(t.possible_variants),
            )
            .key_by(*model_keys)
        )

    ht_syn_grouped = aggregate_for_models(ht_syn)
    ht_syn_grouped = ht_syn_grouped.checkpoint(
        f"{CHECKPOINT_ROOT}/dd_observed_possible.syn.gnomad_pos.ht",
        _read_if_exists=True,
    )
    print("All sites training HT")
    ht_syn_grouped.show()

    high_cov_grouped = aggregate_for_models(high_cov_ht)
    high_cov_grouped = high_cov_grouped.checkpoint(
        f"{CHECKPOINT_ROOT}/dd_observed_possible.syn.gnomad_pos.high_cov.ht",
        _read_if_exists=True,
    )
    print("High coverage training HT")
    high_cov_grouped.show()

    # Build plateau models using build_plateau_models from gnomad_methods.
    plateau_models = ht_syn_grouped.aggregate(
        build_plateau_models(
            mu_snp_expr=ht_syn_grouped.mu_snp,
            observed_variants_expr=ht_syn_grouped.observed_variants[0],
            possible_variants_expr=ht_syn_grouped.possible_variants,
            cpg_expr=ht_syn_grouped.cpg,
        )
    )
    print("Plateau models (all sites):", plateau_models)

    high_cov_plateau_models = high_cov_grouped.aggregate(
        build_plateau_models(
            mu_snp_expr=high_cov_grouped.mu_snp,
            observed_variants_expr=high_cov_grouped.observed_variants[0],
            possible_variants_expr=high_cov_grouped.possible_variants,
            cpg_expr=high_cov_grouped.cpg,
        )
    )
    print("Plateau models (high coverage):", high_cov_plateau_models)

    # -------------------------------------------------------------------
    # Scatter plots: mutation rate vs observed/possible
    # -------------------------------------------------------------------
    plot = Tabs(
        tabs=[
            TabPanel(
                child=hl.plot.scatter(
                    ht_syn_grouped.mu_snp,
                    ht_syn_grouped.observed_variants[0]
                    / hl.float(ht_syn_grouped.possible_variants),
                    label=ht_syn_grouped.mutation_type,
                    xlabel="mu_snp",
                    ylabel="obs/pos",
                ),
                title="All sites",
            ),
            TabPanel(
                child=hl.plot.scatter(
                    high_cov_grouped.mu_snp,
                    high_cov_grouped.observed_variants[0]
                    / hl.float(high_cov_grouped.possible_variants),
                    label=high_cov_grouped.mutation_type,
                    xlabel="mu_snp",
                    ylabel="obs/pos",
                ),
                title="High coverage",
            ),
        ]
    )
    plot_path = "/tmp/dd_mu_vs_obs_pos.html"
    output_file(plot_path)
    save(plot)
    gcs_plot_path = f"{CHECKPOINT_ROOT}/dd_mu_vs_obs_pos.html"
    hl.hadoop_copy(f"file://{plot_path}", gcs_plot_path)
    logger.info("Scatter plot saved to %s", gcs_plot_path)

    # -------------------------------------------------------------------
    # Apply models per-variant
    # -------------------------------------------------------------------
    logger.info("Applying models per-variant...")
    hl._set_flags(use_new_shuffle="1")

    ht = hl.read_table(
        f"{CHECKPOINT_ROOT}/dd_observed_possible.gnomad_pos.exomes_coverage.v2.ht"
    )

    # Convert observed_variants from array to scalar.
    ht = ht.annotate(observed_variants=ht.observed_variants[0])

    # Apply plateau models per-variant using apply_plateau_models.
    cpg_key = hl.struct(cpg=ht.cpg)
    all_ppo = apply_plateau_models(
        ht.mu_snp, hl.literal(plateau_models)[cpg_key]
    )
    high_ppo = apply_plateau_models(
        ht.mu_snp, hl.literal(high_cov_plateau_models)[cpg_key]
    )

    # No coverage correction (cov_corr = 1) for DD analysis.
    ht = ht.annotate(
        all_sites_expected=all_ppo * ht.possible_variants,
        high_cov_expected=high_ppo * ht.possible_variants,
    )

    ht = ht.naive_coalesce(1000).checkpoint(
        f"{CHECKPOINT_ROOT}/dd_per_variant.gnomad_pos.v2.ht",
        _read_if_exists=True,
    )
    ht.show()

    # -------------------------------------------------------------------
    # Aggregate per-transcript obs/exp
    # -------------------------------------------------------------------
    logger.info("Aggregating per-transcript obs/exp...")
    ht = hl.read_table(f"{CHECKPOINT_ROOT}/dd_per_variant.gnomad_pos.v2.ht")

    vep_grouping = (
        "annotation", "modifier", "gene", "gene_id",
        "transcript", "canonical", "mane_select",
    )
    ht = (
        ht.group_by(*vep_grouping)
        .partition_hint(1000)
        .aggregate(
            all_sites=hl.struct(
                observed_variants=hl.agg.sum(ht.observed_variants),
                possible_variants=hl.agg.sum(ht.possible_variants),
                expected_variants=hl.agg.sum(ht.all_sites_expected),
            ),
            high_cov=hl.struct(
                observed_variants=hl.agg.sum(ht.observed_variants),
                possible_variants=hl.agg.sum(ht.possible_variants),
                expected_variants=hl.agg.sum(ht.high_cov_expected),
            ),
        )
    )

    # Compute obs/exp ratio.
    ht = ht.annotate(
        all_sites=ht.all_sites.annotate(
            obs_exp=ht.all_sites.observed_variants / ht.all_sites.expected_variants
        ),
        high_cov=ht.high_cov.annotate(
            obs_exp=ht.high_cov.observed_variants / ht.high_cov.expected_variants
        ),
    )

    ht = ht.naive_coalesce(1000).checkpoint(
        f"{CHECKPOINT_ROOT}/dd_apply_model.gnomad_pos.v2.ht",
        _read_if_exists=True,
    )
    ht.show()

    hl._set_flags(use_new_shuffle=None)


if __name__ == "__main__":
    main()
