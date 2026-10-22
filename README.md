# gnomad-constraint

Genic constraint analysis pipeline for gnomAD. Computes observed/expected ratios, pLI scores, z-scores, and confidence intervals for LoF, missense, and synonymous variants at the gene/transcript level. Current version: **v4.1.1** (GRCh38). Version 2.1.1 (GRCh37) is no longer supported by this pipeline; commit 39928d1 is the last version of the script that supports it.

An interactive map of the pipeline is published at <https://broadinstitute.github.io/gnomad-constraint/>: every step with its flags, the Tables it reads and writes, the function calls behind it (in this repo and in [gnomad_methods](https://github.com/broadinstitute/gnomad_methods)), and the schema, row count, and first rows of each Table. The page is [gnomad_constraint/flowchart/constraint_pipeline.html](./gnomad_constraint/flowchart/constraint_pipeline.html); see [Updating the flowchart](#updating-the-flowchart).

The `gnomad_constraint/experimental/proemis3d/` directory contains the ProEmis3D project for regional missense constraint visualization.

## Project Structure

| Directory | Purpose |
|-----------|---------|
| `gnomad_constraint/pipeline/constraint_pipeline.py` | Main constraint pipeline |
| `gnomad_constraint/utils/constraint.py` | Core utility functions (preprocessing, model building, metrics) |
| `gnomad_constraint/resources/resource_utils.py` | Resource paths and `TableResource` definitions |
| `gnomad_constraint/resources/constants.py` | Pipeline constants |
| `gnomad_constraint/experimental/proemis3d/` | ProEmis3D regional missense constraint |
| `gnomad_constraint/plots/` | R and Python plotting scripts |
| `gnomad_constraint/flowchart/` | Interactive pipeline flowchart and the script that generates it |

## Constraint Pipeline Steps

The main pipeline (`constraint_pipeline.py`) has these steps (each a CLI flag):

| Step | CLI Flag | Function |
|------|----------|----------|
| 1 | `--prepare-context-ht` | Annotate VEP context with methylation, coverage, GERP |
| 2 | `--preprocess-data` | Add observed/possible, mutation rate, and model-group annotations to the annotated context Table |
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
gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint/            # Test root (--test)

# Key outputs:
.../preprocessed_data/gnomad.v4.1.1.annotated_context.ht
.../preprocessed_data/gnomad.v4.1.1.context.preprocessed.ht
.../mutation_rate/gnomad.v4.1.1.mutation_rate.ht
.../training_data/gnomad.v4.1.1.constraint_training.ht
.../models/gnomad.v4.1.1.{plateau|coverage}.he
.../apply_models/transcript_consequences/gnomad.v4.1.1.per_variant_expected.ht
.../apply_models/transcript_consequences/gnomad.v4.1.1.per_variant_expected.aggregated.ht
.../apply_models/transcript_consequences/gnomad.v4.1.1.aggregated_expected.ht
.../apply_models/transcript_consequences/gnomad.v4.1.1.constraint_group.ht
.../metrics/gnomad.v4.1.1.gene_quality_metrics.ht
.../metrics/transcript_consequences/gnomad.v4.1.1.constraint_metrics.ht
.../release/gnomad.v4.1.1.constraint_metrics.{ht|tsv.bgz}
```

All outputs, including the release files, go under the test root with `--test`. See [gnomad_constraint/resources/README.md](gnomad_constraint/resources/README.md) for the full layout and the schema of every output.

### Release constraint metrics Table

Keyed by `(gene, gene_id, transcript, canonical, mane_select)`. One struct per constraint group (`syn`, `mis`, `lof_hc_lc`, `lof`):
- `.obs` — observed variant count
- `.exp` — expected variant count
- `.oe` — observed/expected ratio
- `.oe_ci` — 90% gamma CI around o/e (for `lof`, the LOEUF, with rank and bin fields)
- `.z_raw`, `.z_score` — raw and standardized z-scores
- `.pLI`, `.pNull`, `.pRec` — pLI / null / recessive probabilities, computed from each group's own counts (`lof_hc_lc` and `lof` only)

The internal `metrics/.../constraint_metrics.ht` stores the same metrics in a `constraint_groups` array; see the resources README.

## Updating the flowchart

[generate_flowchart.py](gnomad_constraint/flowchart/generate_flowchart.py) fills three generated blocks in the page:

- **Code**: steps, flags, input and output resources, and the function calls behind each step, read from `constraint_pipeline.py`, `resource_utils.py`, and the utility modules of this repo and gnomad_methods. Prose the code cannot supply (what each step does to the data, pitfalls, Table summaries) lives in [notes.py](gnomad_constraint/flowchart/notes.py).
- **Notes**: `notes.py`, checked against the code so a renamed step, function, flag, or resource fails instead of going stale.
- **Tables**: the schema, row count, globals, and first rows of every pipeline input and output, read on Dataproc.

After changing the pipeline, commit, then run from the repository root with gnomad_methods and gnomad_qc importable:

```bash
python -m gnomad_constraint.flowchart.generate_flowchart --extract
```

`--extract` records the commit of both repositories so every source link is a permalink, and refuses to read files with uncommitted changes (`--allow-dirty` overrides this, and the page then shows a warning). CI runs `--check`, which fails when the page no longer matches the code or `notes.py`.

After rerunning pipeline steps, refresh the Table panels (build the single zip as in [CLAUDE.md](CLAUDE.md)):

```bash
hailctl dataproc submit <cluster> \
  <repo>/gnomad_constraint/flowchart/generate_flowchart.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --collect --out gs://gnomad-tmp-4day/constraint_flowchart/tables.json
gsutil cp gs://gnomad-tmp-4day/constraint_flowchart/tables.json .
python -m gnomad_constraint.flowchart.generate_flowchart --inject-tables tables.json
```

[pages.yml](.github/workflows/pages.yml) publishes the page to GitHub Pages whenever it changes on `main`.

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
