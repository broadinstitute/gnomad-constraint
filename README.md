# gnomad-constraint

Genic constraint analysis pipeline for gnomAD. Computes observed/expected ratios, pLI scores, z-scores, and confidence intervals for LoF, missense, and synonymous variants at the gene/transcript level. Current version: **v4.1.1** (GRCh38). Historical version 2.1.1 (GRCh37) is also supported.

An overview of the pipeline and functions that are used can be found in [/gnomad_constraint/flowchart/constraint_pipeline_v4.pdf](./gnomad_constraint/flowchart/constraint_pipeline_v4.pdf). Note that many functions are imported from the [gnomad_methods](https://github.com/broadinstitute/gnomad_methods) repo.

The `gnomad_constraint/experimental/proemis3d/` directory contains the ProEmis3D project for regional missense constraint visualization.

## Project Structure

| Directory | Purpose |
|-----------|---------|
| `gnomad_constraint/pipeline/constraint_pipeline.py` | Main constraint pipeline |
| `gnomad_constraint/pipeline/constraint_pipeline_complex.py` | Complex region constraint pipeline |
| `gnomad_constraint/utils/constraint.py` | Core utility functions (preprocessing, model building, metrics) |
| `gnomad_constraint/utils/constraint_complex.py` | Complex region constraint utilities |
| `gnomad_constraint/resources/resource_utils.py` | Resource paths and `TableResource` definitions |
| `gnomad_constraint/resources/constants.py` | Pipeline constants |
| `gnomad_constraint/experimental/proemis3d/` | ProEmis3D regional missense constraint |
| `gnomad_constraint/plots/` | R and Python plotting scripts |

## Constraint Pipeline Steps

The main pipeline (`constraint_pipeline.py`) has these steps (each a CLI flag):

| Step | CLI Flag | Function |
|------|----------|----------|
| 1 | `--prepare-context-ht` | Annotate VEP context with methylation, coverage, GERP |
| 2 | `--preprocess-data` | Add VEP context annotations to exome/genome tables, prepare for constraint |
| 3 | `--calculate-gerp-cutoffs` | Optional: compute GERP percentile cutoffs |
| 4 | `--calculate-mutation-rate` | Compute baseline mutation rate per substitution/context |
| 5 | `--create-training-set` | Count observed + possible variants at synonymous sites |
| 6 | `--build-models` | Build plateau and coverage regression models |
| 7a | `--apply-models-per-variant` | Per-variant: apply models to compute expected counts per variant |
| 7b | `--apply-models-aggregated` | Aggregated: aggregate counts first, then apply models (alternative to 7a+8+9) |
| 8 | `--aggregate-per-variant-expected` | Aggregate per-variant expected counts by transcript/consequence |
| 9 | `--aggregate-by-constraint-groups` | Group aggregated counts into constraint groups (lof, mis, syn) |
| 10 | `--compute-gene-quality-metrics` | Compute per-transcript gene quality metrics (coverage, MQ, segdup, LCR) |
| 11 | `--compute-constraint-metrics` | Compute pLI, z-scores, o/e with CIs |
| 12 | `--prepare-release` | Format constraint metrics for public release |
| 13 | `--export-release-tsv` | Export release table to TSV |
| 14 | `--export-release-downsampling-tsv` | Export downsampling constraint metrics to TSV |

Steps 7a/8/9 (per-variant path) and 7b (aggregated path) are alternative ways to go from models to constraint groups. The per-variant path applies coverage correction per-variant then aggregates; the aggregated path aggregates counts first then applies models. Use `--use-aggregated-expected` with `--aggregate-by-constraint-groups` to read from the aggregated path output instead of the per-variant path.

### Key Resource Paths (v4.1.1)

```
gs://gnomad/v4.1.1/constraint/                              # Production root
gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint/            # Test root

# Key outputs:
.../preprocessed_data/annotated_context.ht
.../preprocessed_data/gnomad.v4.1.1.{context|exomes|genomes}.preprocessed.{region}.ht
.../mutation_rate/gnomad.v4.1.1.mutation_rate.ht
.../training_data/gnomad.v4.1.1.constraint_training.{region}.ht
.../models/gnomad.v4.1.1.{plateau|coverage}.{region}.he
.../apply_models/transcript_consequences/gnomad.v4.1.1.per_variant_expected.{region}.ht
.../apply_models/transcript_consequences/gnomad.v4.1.1.aggregated_expected.{region}.ht
.../apply_models/transcript_consequences/gnomad.v4.1.1.apply.{region}.ht
.../constraint_groups/transcript_consequences/gnomad.v4.1.1.constraint_groups.{region}.ht
.../metrics/gnomad.v4.1.1.constraint_metrics.ht
```

### constraint_metrics Table

Keyed by `(gene, transcript, canonical)` (and optionally `mane_select`, `gene_id`).

Output struct per annotation category (`lof`, `mis`, `syn`):
- `.obs` — observed variant count
- `.exp` — expected variant count
- `.oe` — observed/expected ratio
- `.oe_ci` — 90% CI around o/e
- `.z_raw` — raw z-score
- `.pLI` — probability of loss-of-function intolerance (LoF only)
- `.pNull`, `.pRec` — null/recessive probabilities (LoF only)

## Missense Score Percentile Analysis

`gnomad_constraint/plots/determine_missense_score_percentiles.py` computes per-percentile depletion of missense variants binned by missense prediction scores.

### Scores analyzed

ProteinMPNN, ESM, REVEL, RASP, AM, MisFit, PolyPhen, CPT1, popEVE, EVE, MPC, CADD, GPN-MSA

### Pipeline steps (CLI flags)

| Step | Flag |
|------|------|
| 1 | `--preprocess-scores` |
| 2 | `--compute-percentiles` |
| 3 | `--annotate-constraint-data` |
| 4 | `--aggregate-by-transcript` |
| 5 | `--compute-cumulative` |
| 6 | `--export-percentile-summary` |
| 7 | `--export-matched-plof-summary` |

Step 7 computes matched pLoF o/e per missense percentile bin. It computes adj_r-corrected gene-level pLoF from the per-SNV table (`--constraint-ht-path`), not from the pre-computed constraint_metrics table (which lacks adj_r for pLoF).

## Dependencies

- **hail** — distributed genomics framework
- **numpy**, **pandas**, **scipy** — numerical/statistical
- **gnomad** (gnomad_methods) — shared gnomAD utilities
- **gnomad_qc** — gnomAD QC pipeline resources
