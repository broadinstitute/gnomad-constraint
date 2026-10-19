# gnomad-constraint

Genic constraint analysis pipeline for gnomAD. Computes observed/expected ratios, pLI scores, z-scores, and confidence intervals for LoF, missense, and synonymous variants at the gene/transcript level. Current version: **v4.1.1** (GRCh38). Historical version 2.1.1 (GRCh37) is also supported.

An overview of the pipeline and functions that are used can be found in [/gnomad_constraint/flowchart/constraint_pipeline_v4.pdf](./gnomad_constraint/flowchart/constraint_pipeline_v4.pdf). Note that many functions are imported from the [gnomad_methods](https://github.com/broadinstitute/gnomad_methods) repo.

The `gnomad_constraint/experimental/proemis3d/` directory contains the ProEmis3D project for regional missense constraint visualization.

## Project Structure

| Directory | Purpose |
|-----------|---------|
| `gnomad_constraint/pipeline/constraint_pipeline.py` | Main constraint pipeline |
| `gnomad_constraint/pipeline/constraint_pipeline_complex.py` | Complex region constraint pipeline |
| `gnomad_constraint/pipeline/regional_rate_pipeline.py` | Regional mutation rate (rr) models and per-window rr |
| `gnomad_constraint/pipeline/prepare_syn_dnms.py` | Combine DNM cohorts into the synonymous DNM training Table |
| `gnomad_constraint/pipeline/build_xy_features.py` | Build regional genomic features for chrX/chrY windows |
| `gnomad_constraint/utils/constraint.py` | Core utility functions (preprocessing, model building, metrics) |
| `gnomad_constraint/utils/constraint_complex.py` | Complex region constraint utilities |
| `gnomad_constraint/utils/regional_rate.py` | Regional mutation rate (rr) utilities |
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
| 8 | `--aggregate-per-variant-expected` | Apply the regional rr adjustment (default) and aggregate per-variant expected counts by transcript/consequence |
| 9 | `--aggregate-by-constraint-groups` | Group aggregated counts into constraint groups (lof, mis, syn) |
| 10 | `--compute-gene-quality-metrics` | Compute per-transcript gene quality metrics (coverage, MQ, segdup, LCR) |
| 11 | `--compute-constraint-metrics` | Compute pLI, z-scores, o/e with CIs |
| 12 | `--prepare-release` | Format constraint metrics for public release |
| 13 | `--export-release-tsv` | Export release table to TSV |
| 14 | `--export-release-downsampling-tsv` | Export downsampling constraint metrics to TSV |

Steps 7a/8/9 (per-variant path) and 7b (aggregated path) are alternative ways to go from models to constraint groups. The per-variant path applies coverage correction per-variant then aggregates; the aggregated path aggregates counts first then applies models. Use `--use-aggregated-expected` with `--aggregate-by-constraint-groups` to read from the aggregated path output instead of the per-variant path. The regional rr adjustment (see below) is only available on the per-variant path; use `--skip-regional-rr` with the aggregated path.

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
- `.unadj` — `exp`, `oe`, `oe_ci`, `z_raw`, and `z_score` without the regional rr adjustment

## Regional Mutation Rate (rr) Adjustment

Expected counts are multiplied per variant by a regional mutation rate factor (rr) of the variant's 1kb window and strand-collapsed trinucleotide context, a port of the [Gnocchi](https://github.com/atgu/gnomad_nc_constraint) regional adjustment retrained on synonymous de novo mutations (DNMs). For each context, a logistic model separates synonymous DNM sites from non-mutated synonymous sites using 13 regional genomic features at 4 flanking scales (features are selected with a univariate screen), and `rr = P(window) / P(average window)`.

rr is applied **after model calibration** in `--aggregate-per-variant-expected`, using the existing per-variant expected Table (`.../apply_models/transcript_consequences/gnomad.v4.1.1.per_variant_expected.ht`); the plateau and coverage models are not re-fit or re-applied. By default, `exp` and all downstream metrics (oe, CIs, z-scores, LOEUF bins, pLI) use the rr-adjusted expected counts, and the unadjusted values are kept alongside (`unadj_expected`, `unadj_oe_info`, and the release `unadj` struct). rr is 1 for sites without a feature window.

### Run order

1. `build_xy_features.py --validate` then `--build-xy-features` (chrX/chrY features; copies and extracts the published `gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13.tar.gz` source tracks to `--tracks-cache-dir`, or use an already extracted copy with `--tracks-dir`).
2. `prepare_syn_dnms.py --prepare-syn-dnms` (combine 8 DNM cohorts, keep synonymous DNMs, always deduplicate by locus/alleles). The ASD and DD DNMs come from the Satterstrom 2026 `satterstrom_2026_asd` and `satterstrom_2026_ddd` files (trios shared with the ASD cohort are removed from the DD file), which replace the older DDD + ASD DNM Table.
3. `regional_rate_pipeline.py`:

| Step | CLI Flag | Function |
|------|----------|----------|
| 1 | `--import-genomic-features` | Import the 1kb feature matrices (`--include-xy-features`, `--extra-features-paths`) |
| 2 | `--create-rr-training-set` | Synonymous DNM sites + sampled non-mutated synonymous background |
| 3 | `--train-rr-models` | Univariate feature screen and per-context models per `--feature-pval-cutoffs` (`--cv-chrom-split` for held-out metrics) |
| 4 | `--compute-rr` | rr per window/context (`--apply-rr-safeguards`, `--mask-features`, `--out-of-fold`, `--mask-ablation`) |
| 5 | `--evaluate-rr-cutoffs` | Synonymous calibration (dispersion, outliers, o/e by LOEUF decile) of out-of-fold rr per cutoff (LOEUF deciles from `--loeuf-ht-path`, default the v4.1.2 constraint metrics Table) |
| 6 | `--compare-rr-variants` | Synonymous calibration and per-context mean rr of the rr variants in `--rr-variants` (e.g. legacy `syn_adj_r` vs. the new rr with or without clipping, normalization, masking, or out-of-fold models) |

Synonymous calibration is evaluated on Ensembl MANE Select transcripts only (the RefSeq MANE Select transcripts would count each gene twice), and the rr normalization factors are computed on the same sites.

4. `constraint_pipeline.py --aggregate-per-variant-expected` onward (rr applied by default; `--regional-rr-type {calibrated,raw}`, `--regional-rr-feature-cutoff`, `--regional-rr-corrected-window-coords`, `--regional-rr-xy-masked`, `--regional-rr-ht-path`, `--skip-regional-rr`).

### Options

- **Safeguards** (`--apply-rr-safeguards`): clip rr to `[--rr-min, --rr-max]` (default `[0.2, 5.0]`) and divide rr by its per-context expected-weighted mean over exome synonymous MANE Select sites. Without this flag, rr is written as trained, which reproduces the legacy `syn_adj_r` Tables exactly (with legacy window coordinates, the `bonferroni` cutoff, and `statsmodels==0.14.4`, which requires `scipy<1.16`).
- **Window coordinates** (`--corrected-window-coords`): use `[start + 1, end]` for the BED windows instead of the legacy `[start, end)` (1bp offset).
- **chrX/chrY masking** (`--mask-features`): on chrX/chrY windows only, pin the features without X/Y source data (`--masked-features`, default `recomb_male`, `recomb_female`, `met_sperm`, `Nucleosome`, `cDNM_maternal_05M`, `cDNM_paternal_05M`) to the genome-average window. Autosomes are never masked.
- **Feature selection cutoff** (`--feature-pval-cutoffs`): `bonferroni` (legacy: `0.05/4/13`, `0.05/4/8` for CpG contexts) and/or raw p-value cutoffs. Compare cutoffs with `--train-rr-models --cv-chrom-split`, `--compute-rr --out-of-fold`, and `--evaluate-rr-cutoffs`.
- **rr variants** (`--rr-variants`): each spec is a source (`legacy` for the legacy `syn_adj_r`, `legacy_noncoding` for the legacy non-coding `adj_r`, `path=<rr Table>`, or `cutoff=<token>` for rr built from the trained models) followed by comma-separated modifiers (`clip`, `norm`, and for `cutoff=` also `mask` and `oof`), e.g. `cutoff=bonferroni,clip,norm,mask`. The default list adds one change of the new code at a time to the legacy `syn_adj_r` baseline.

### Output roots

All scripts accept `--output-root` (e.g. `gs://gnomad-tmp/rp`) to write outputs under `{output_root}/gnomad_v{version}[_testing]/constraint/...` instead of the production root. In `constraint_pipeline.py`, `--output-root` only applies to `--aggregate-per-variant-expected` and downstream steps (including release and export); the upstream production Tables are read-only inputs. Hail temporary files are written to `gs://gnomad-tmp-4day`.

### Dataproc

The rr models need `statsmodels` on the cluster. To reproduce the legacy `syn_adj_r` exactly, start the cluster with `--packages "statsmodels==0.14.4,scipy==1.15.3"`, then submit with the single-zip `--pyfiles` workaround (see `CLAUDE.md`), e.g.:

```bash
hailctl dataproc submit <CLUSTER> gnomad_constraint/pipeline/regional_rate_pipeline.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --output-root gs://gnomad-tmp/rp --create-rr-training-set
```

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
- **statsmodels**, **scikit-learn** — regional rr models (`statsmodels==0.14.4` with `scipy<1.16` to reproduce the legacy `syn_adj_r` exactly)
