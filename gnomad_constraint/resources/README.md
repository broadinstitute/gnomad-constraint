# gnomad_constraint resources

This package wires every Hail Table produced by the constraint pipeline to a
canonical GCS path. The functions live in [resource_utils.py](resource_utils.py)
and the constants they consume live in [constants.py](constants.py).

This document is the schema reference for every **pipeline output** of the
constraint pipeline and of `import_data.py` and `regional_rate_pipeline.py`. The
schemas were first captured with
`ht.describe()` from the v4.1.1 tables on GCS and have been updated by hand as the
pipeline changed, so the released v4.1.1 tables differ from them in places. For upstream **input**
resources (VEP context, sites,
coverage, methylation, GERP, GENCODE, `all_sites_an`) see
[`gnomad.resources.grch38`](https://github.com/broadinstitute/gnomad_methods/tree/main/gnomad/resources/grch38).

> **Scope.** This README documents v4.1.1. Older versions (`4.0`, `4.1`)
> follow the same general layout but partition some artifacts by
> genomic region (e.g. `*.autosome_par.ht`, `*.chrx_nonpar.ht`,
> `*.chry_nonpar.ht`). v4.1.1 keeps all regions in one Table and labels each
> site with a `genomic_region` annotation (`autosome_or_par`, `chrx_nonpar`, or
> `chry_nonpar`).

---

## Path conventions

Built by [`get_constraint_root`](resource_utils.py) and
[`get_constraint_data`](resource_utils.py):

```
gs://gnomad/v{version}/constraint/
├── preprocessed_data/
│   ├── gnomad.v{version}.annotated_context.ht
│   └── gnomad.v{version}.context.preprocessed.ht
├── mutation_rate/
│   └── gnomad.v{version}.mutation_rate.ht
├── training_data/
│   ├── gnomad.v{version}.constraint_training.ht
│   └── gnomad.v{version}.constraint_training.tsv.bgz
├── models/
│   ├── gnomad.v{version}.plateau.he
│   └── gnomad.v{version}.coverage.he
├── imported_data/
│   ├── gnomad.v{version}.dnm.ht
│   └── gnomad.v{version}.genomic_features13_1kb.ht
├── regional_rate/
│   ├── gnomad.v{version}.rr_training.ht
│   ├── gnomad.v{version}.rr_penalty_tuning.ht
│   ├── gnomad.v{version}.rr_models.ht
│   └── gnomad.v{version}.syn_dnm_rr_1kb.ht
├── apply_models/transcript_consequences/
│   ├── gnomad.v{version}.per_variant_expected.ht
│   ├── gnomad.v{version}.per_variant_expected.aggregated.ht
│   ├── gnomad.v{version}.aggregated_expected.ht
│   └── gnomad.v{version}.constraint_group.ht
├── metrics/
│   ├── gnomad.v{version}.gene_quality_metrics.ht
│   └── transcript_consequences/
│       ├── gnomad.v{version}.constraint_metrics_pre_rank.ht
│       └── gnomad.v{version}.constraint_metrics.ht
└── release/
    ├── gnomad.v{version}.constraint_metrics.ht
    ├── gnomad.v{version}.constraint_metrics.tsv.bgz
    ├── gnomad.v{version}.constraint_metrics.downsampling.tsv.bgz
    ├── gnomad.v{version}.loeuf_percentile_thresholds.tsv
    ├── gnomad.v{version}.mutation_rate.ht
    └── gnomad.v{version}.mutation_rate.tsv
```

Tests use `gs://gnomad-tmp/gnomad_v{version}_testing/constraint/...` (`test=True`)
and intermediate checkpoints use `gs://gnomad-tmp/gnomad_v{version}/constraint/...`
(`temp=True`). These rules apply to every output, including the release files.
`--directory-post-fix` changes `constraint/` to `constraint_{post_fix}/`, and
`--path-post-fix` adds `.{post_fix}` before the file extension of training, model,
apply-models, metrics, and constraint-metrics release outputs. `import_data.py` takes
neither post-fix. `regional_rate_pipeline.py` applies both to its own outputs and to
the per-variant expected Table it normalizes rr on, and
`constraint_pipeline.py --apply-regional-rr` reads the rr Table with both, so run the
two scripts with the same post-fixes.

---

## Shared structures

Several tables carry identical **pipeline-parameter globals** that capture the
arguments used to produce them. They are written by
[`calculate_mu_by_downsampling`](../utils/constraint.py), [`build_models`](../utils/constraint.py),
and the apply-models steps:

### `calculate_mu_globals` — `struct`
Parameters passed to `calculate_mu_by_downsampling`.
- `freq_meta: array<dict<str,str>>` — frequency metadata entries (one per
  freq-array index): `{group: adj}` for the global adj entry plus one
  `{downsampling, group, pop}` entry per downsampling × genetic-ancestry-group.
- `ac_cutoff: int32` — variants with `AC > ac_cutoff` are excluded from mu
  calculation (typical value: 5).
- `min_cov: int32` — minimum mean genome coverage filter for mu sites (15).
- `max_cov: int32` — maximum mean genome coverage filter (60).
- `gerp_lower_cutoff: float64` — lower GERP bound for the mu site set (−3.9885).
- `gerp_upper_cutoff: float64` — upper GERP bound (2.6607).
- `genetic_ancestry_groups: array<str>` — gen-anc groups used for downsampling
  (`global, afr, amr, eas, nfe, sas`).
- `downsampling_level: int32` — the downsampling size used for the canonical
  mu calculation (1000).
- `downsampling_idx: int32` — index into the parallel freq arrays whose entry
  was used as `mu_snp` (the global×1000 entry).
- `most_severe_consequence: array<str>` — VEP `most_severe_consequence` values
  retained as putatively neutral for mu fitting
  (`intron_variant`, `intergenic_variant`).

### `build_models_globals` — `struct`
Parameters passed to the plateau/coverage model fit.
- `synonymous_transcript_filter_field: str` — transcript filter used when
  selecting synonymous training variants (e.g. `canonical`).
- `low_cov_cutoff: int32` — coverage below this is treated as the low-coverage
  regime (19).
- `high_cov_cutoff: int32` — coverage at or above this is the high-coverage
  regime that the plateau model is fit to (90).
- `upper_cov_cutoff: int32` — upper bound for the coverage model fit (nullable).
- `skip_coverage_model: bool` — whether the coverage model fit was skipped.

### `apply_models_globals` — `struct`
Carries the previous two plus the fitted models (after `--apply-models-*`).
- `low_cov_cutoff`, `high_cov_cutoff`, `skip_coverage_model` — as above.
- `plateau_models: dict<struct{cpg: bool, genomic_region: str}, array<array<float64>>>`
  — fitted plateau slope/intercept per (cpg, genomic_region). Each value is
  an array parallel to `mutation_rate.mu` (one fit per downsampling).
- `coverage_model: array<float64>` — `[intercept, slope]` of the coverage
  correction regression.
- `log10_coverage: bool` — whether the coverage model is in log10 space.
- `groupings: tuple(str×7)` or `array<str>` — the field names used to group
  variants when applying models (annotation, modifier, gene, gene_id,
  transcript, canonical, mane_select). Stored as a tuple on the per-variant
  path, as an array on the aggregated path.

### Other shared globals
- `exomes_freq_meta: array<dict<str,str>>` — frequency metadata of the exomes
  freq array (one entry per `(downsampling, gen_anc)`).
- `genetic_ancestry_groups: array<str>` — copy of the gen-anc groups list.
- `downsamplings: array<int32>` — ordered list of downsampling sizes the freq
  array covers (e.g. `[10, 20, ..., 1000, ...]`).
- `max_af: float64` — AF cap for the "observed" variant filter.

---

# Pipeline output tables

Each section below lists the **canonical path**, the **key**, the **globals**
and the **row fields** with sub-struct expansion. Producing function is the
`--<step>` flag from [constraint_pipeline.py](../pipeline/constraint_pipeline.py).

---

## 1. `annotated_context.ht`

The fully annotated universe of all possible single-nucleotide substitutions
in the genome, produced by `--prepare-context-ht` /
[`prepare_ht_for_constraint_calculations`](../utils/constraint.py).

- **Path:** `gs://gnomad/v{version}/constraint/preprocessed_data/gnomad.v{version}.annotated_context.ht`
- **Resource fn:** [`get_annotated_context_ht`](resource_utils.py)
- **Key:** `locus`, `alleles`

### Globals
- `grange: array<int32>` — methylation-level grange bin edges (10 bins).
- `vep_help: str` — captured `vep --help` text from the run that annotated VEP.
- `vep_config: str` — JSON-serialized VEP runner config (cache version, plugins,
  command).
- `version: str` — gnomAD version string.
- `an_globals: struct{exomes, genomes}` — AN strata metadata (from
  [`all_sites_an`](https://github.com/broadinstitute/gnomad_methods)) for each data type:
  - `strata_sample_count: array<int32>` — sample count for each strata index.
  - `strata_meta: array<dict<str,str>>` — strata key (e.g. `{group: adj, gen_anc: afr}`).
- `freq_globals: struct{exomes, genomes}` — frequency strata metadata:
  - `exomes.freq_meta_sample_count: array<int32>`
  - `exomes.freq_meta: array<dict<str,str>>`
  - `genomes.freq_meta: array<dict<str,str>>`

### Rows
- `locus: locus<GRCh38>`
- `alleles: array<str>` — `[ref, alt]` of one SNV.
- `context: str` — trinucleotide context (ref-strand-collapsed via
  [`collapse_strand`](../utils/constraint.py): rows with G/T ref are reverse-complemented and
  `was_flipped` is set true).
- `vep: struct` — minimal VEP output kept for constraint:
  - `most_severe_consequence: str`
  - `transcript_consequences: array<struct>` — per-transcript VEP entries with
    `transcript_id, gene_id, gene_symbol, biotype, most_severe_consequence,
    mane_select, canonical, lof, lof_flags, sift_score, polyphen_score,
    domains, uniprot_isoform, amino_acids, codons`.
- `ref: str`, `alt: str` — ref-strand-collapsed alleles.
- `was_flipped: bool` — `true` if the original ref was G or T and
  ref/alt/context were reverse-complemented.
- `transition: bool` — purine↔purine or pyrimidine↔pyrimidine.
- `cpg: bool` — variant is in a CpG dinucleotide context.
- `mutation_type: str` — coarse mutation type (e.g. `CpG`, `non-CpG transition`,
  `transversion`).
- `mutation_type_model: str` — mutation-type label used as the model group.
- `methylation_level: int32` — discretized methylation level (0–`len(grange)-1`).
- `gerp: float64` — GERP RS score at the locus.
- `coverage: struct{exomes, genomes}` — coverage from gnomAD coverage HT:
  - `<dt>.mean: float64`, `<dt>.median_approx: int32`.
- `AN: struct{exomes: int64, genomes: int64}` — global-adj allele number from
  [`all_sites_an`](https://github.com/broadinstitute/gnomad_methods).
- `freq: struct{exomes, genomes}` — per-downsampling allele freqs from the
  sites table. Each entry is `array<struct{AC, AF, AN, homozygote_count}>`
  parallel to `freq_globals.<dt>.freq_meta`. (Note: exomes
  `homozygote_count` is `int64`, genomes `int32`.)
- `filters: struct{exomes: set<str>, genomes: set<str>}` — site filters
  inherited from the corresponding release sites table.
- `genomic_region: str` — `autosome_or_par`, `chrx_nonpar`, or `chry_nonpar`.
- `adj_r: float64` — per-context regional depletion correction (see
  `adj_r.ht`).
- `syn_adj_r: float64` — synonymous-DNM variant of `adj_r` (see `syn_adj_r.ht`).
- `sfs_bin: int32` — site-frequency-spectrum bin assigned via
  [`annotate_sfs_bin`](../utils/constraint.py) from `SFS_BIN_CUTOFFS` (0 if AF
  is missing; otherwise index of first cutoff `af_expr <= cutoff` is true).

---

## 2. `context.preprocessed.ht`

Context joined with exomes/genomes frequency, coverage, AN, and the
pre-computed mu inputs (`compute_mu`, `calibrate_mu` structs). Produced by the
`preprocess_data` step / [`preprocess_data`](../utils/constraint.py) — the
shared upstream input for nearly every downstream step.

- **Path:** `gs://gnomad/v{version}/constraint/preprocessed_data/gnomad.v{version}.context.preprocessed.ht`
- **Resource fn:** [`get_preprocessed_ht`](resource_utils.py)
- **Key:** `locus`, `alleles`

### Globals
`calculate_mu_globals`, `build_models_globals`, `apply_models_globals`,
`exomes_freq_meta`, `genetic_ancestry_groups`, `downsamplings`, `max_af` —
see [Shared structures](#shared-structures).

### Rows
Carries every annotated_context row field *except* `freq` (which is split
into `compute_mu.genomes_freq` / `calibrate_mu.exomes_freq`), plus:
- `exomes_coverage: int32` — exome coverage from `get_exome_coverage_expr`,
  set by `--exome-coverage-metric`: the percent of exome samples with a
  genotype (`AN_percent`, the default), the median coverage (`median`), or the
  allele number (`AN`). It isn't clamped; the model coverage cutoffs are
  applied to it when assigning `build_model` and `apply_model`.
- `compute_mu: struct` — fields used to fit the mutation-rate model on the
  *genome* SNVs:
  - `genomes_freq: array<struct{AC, AF, AN, homozygote_count}>` — genomes freq
    array (parallel to `calculate_mu_globals.freq_meta`). No later step reads
    it; it is kept for analyses that start from this Table.
  - `observed_variants: array<int32>` — observed counts (one per freq entry).
  - `possible_variants: int32` — count of possible SNVs at this site (always 1
    for context rows; used by aggregation).
- `calibrate_mu: struct` — fields used to *apply* mu (i.e. calibrate by
  comparison to *exome* SNVs):
  - `exomes_freq: array<struct{AC, AF, AN, homozygote_count}>` — exomes freq
    array (parallel to `exomes_freq_meta`), missing for variants that fail the
    exome filters. No later step reads it; it is kept for analyses that start
    from this Table.
  - `observed_variants: array<int32>`
  - `possible_variants: int32`
  - `build_model: struct{high_or_low_coverage: str, model_group: struct{cpg: bool, genomic_region: str}}`
    — which plateau model bucket this row was fit into.
  - `apply_model: struct{high_or_low_coverage: str, model_group: struct{cpg: bool, genomic_region: str}}`
    — which plateau model bucket this row will be evaluated against.

---

## 3. `mutation_rate.ht`

Per-context mutation rate, produced by `--calculate-mutation-rate` /
[`calculate_mu_by_downsampling`](../utils/constraint.py).

- **Path:** `gs://gnomad/v{version}/constraint/mutation_rate/gnomad.v{version}.mutation_rate.ht`
- **Resource fn:** [`get_mutation_ht`](resource_utils.py)
- **Key:** `context`, `ref`, `alt`, `methylation_level` (the `MU_GROUPING`)

### Globals
`calculate_mu_globals`, `build_models_globals`, `apply_models_globals`,
`exomes_freq_meta`, `genetic_ancestry_groups`, `downsamplings`, `max_af` —
see [Shared structures](#shared-structures).

### Rows
- `context: str`, `ref: str`, `alt: str`, `methylation_level: int32` — key.
- `observed_variants: array<int64>` — observed counts per freq-meta entry.
- `possible_variants: int64` — count of possible SNVs of this trinucleotide
  context/substitution.
- `proportion_observed: array<float64>` — `observed / possible` per freq entry.
- `mu: array<float64>` — scaled mutation rate per freq entry.
- `mu_snp: float64` — scalar mutation rate at the canonical
  `(global, downsampling=downsampling_level)` index (typically global×1000).
- `transition: bool`, `cpg: bool`, `mutation_type: str`, `mutation_type_model: str`
  — mutation-type annotations propagated from
  [`annotate_mutation_type`](../utils/constraint.py).

---

## 4. `constraint_training.ht`

Training set for plateau / coverage model fit. Produced by
`--create-training-set` from synonymous high-coverage SNVs.

- **Path:** `gs://gnomad/v{version}/constraint/training_data/gnomad.v{version}.constraint_training.ht`
- **Resource fn:** [`get_training_dataset`](resource_utils.py)
- **Key:** `context, ref, alt, methylation_level, cpg, transition, mutation_type, mutation_type_model, genomic_region, build_model, exomes_coverage`

### Globals
Same shared pipeline parameter globals as `mutation_rate.ht`.

### Rows
- `context`, `ref`, `alt`, `methylation_level` — substitution / methylation key.
- `cpg`, `transition`, `mutation_type`, `mutation_type_model` — mutation-type annotations.
- `genomic_region: str` — `autosome_or_par` / `chrx_nonpar` / `chry_nonpar`.
- `build_model: struct{high_or_low_coverage: str, model_group: struct{cpg: bool, genomic_region: str}}`
  — the (coverage-band × cpg × region) bucket whose plateau model this row trains.
- `exomes_coverage: int32` — exome coverage, as in §2; rows are grouped by it.
- `observed_variants: array<int64>` — aggregated observed SNV counts per freq entry.
- `possible_variants: int64` — aggregated possible SNV count.
- `mu_snp: float64` — joined-in mu_snp for this (context, ref, alt, methylation).

A `tsv.bgz` mirror of this table is written to the same directory.

---

## 5. `plateau.he` and `coverage.he`

Fitted Hail Expressions produced by `--build-models` /
[`build_models`](gnomad.utils.constraint:build_models). These are *not* tables — they are
pickled Hail expressions reread by `apply_models_globals`.

- **Paths:**
  - `gs://gnomad/v{version}/constraint/models/gnomad.v{version}.plateau.he`
  - `gs://gnomad/v{version}/constraint/models/gnomad.v{version}.coverage.he`
- **Resource fn:** [`get_models`](resource_utils.py)

### Types and meaning
- `plateau.he : dict<struct{cpg: bool, genomic_region: str}, array<array<float64>>>`
  — for each (cpg, genomic_region) bucket, a list of `[intercept, slope]`
  pairs parallel to `exomes_freq_meta` (one entry per exome frequency index).
  The global-adj fit is at the index of `{group: adj}` in `exomes_freq_meta`;
  `calculate_mu_globals.downsampling_idx` instead indexes the genome mutation rates.
- `coverage.he : struct{coefficients: array<float64>, log10_coverage: bool}` —
  `coefficients` is the `[intercept, slope]` of the linear regression mapping
  exome coverage (`--exome-coverage-metric`), or its log10 when
  `log10_coverage` is True, to the multiplicative coverage correction in the
  low-coverage regime. The apply steps use the saved scale. Models built before
  the scale was saved, including v4.1.1, are only the `array<float64>` of
  coefficients, and are applied on the scale set by
  `--use-logarithmic-coverage-model`. For these, the apply steps raise an error
  if the model's correction at the high coverage cutoff is further from 1 on
  that scale than on the other.

These expressions are read back into `apply_models_globals.plateau_models`,
`apply_models_globals.coverage_model` (the coefficients), and
`apply_models_globals.log10_coverage` on every downstream apply step.

---

## 6. `per_variant_expected.ht`

Per-SNV expected counts, produced by `--apply-models-per-variant`. One row per
locus×allele pair, exploded through `transcript_consequences`. It keeps every
preprocessed field and every consequence, including the non-coding ones that
`--aggregate-per-variant-expected` drops, because downstream analyses read this
Table directly.

- **Path:** `gs://gnomad/v{version}/constraint/apply_models/transcript_consequences/gnomad.v{version}.per_variant_expected.ht`
- **Resource fn:** [`get_per_variant_expected_dataset`](resource_utils.py)
- **Key:** `locus`, `alleles`

### Globals
`calculate_mu_globals`, `build_models_globals`, `apply_models_globals` (now
populated with `plateau_models` and `coverage_model`), `exomes_freq_meta`,
`genetic_ancestry_groups`, `downsamplings`, `max_af`.

### Rows
All preprocessed-context fields (`locus, alleles, context, ref, alt,
was_flipped, transition, cpg, mutation_type, mutation_type_model,
methylation_level, gerp, coverage, AN, filters, genomic_region, adj_r,
syn_adj_r, sfs_bin, exomes_coverage, compute_mu, calibrate_mu`) plus the
applied-model fields:
- `annotation: str` — VEP `most_severe_consequence` for this transcript row.
- `modifier: str` — finer-grained consequence sub-class (e.g. LOFTEE HC/LC,
  `missense_variant` modifier, etc.).
- `gene: str`, `gene_id: str`, `transcript: str` — VEP transcript identity.
- `canonical: bool`, `mane_select: bool` — transcript flags.
- `mu_snp: float64` — joined per-context mu (scalar).
- `mu: float64` — `mu_snp * possible_variants * coverage_correction` (the
  per-variant mu mass, with the coverage correction).
- `predicted_proportion_observed: array<float64>` — plateau-applied predicted
  observed proportion per freq entry (one element per `freq_meta` row).
- `expected_variants: array<float64>` — `predicted_proportion_observed *
  possible * coverage_correction` per freq entry.
- `coverage_correction: float64` — coverage-model multiplicative correction
  evaluated at this row's `exomes_coverage`.

---

## 7. `per_variant_expected.aggregated.ht`

Sum of `per_variant_expected.ht` over variants within each
`(annotation, modifier, gene, gene_id, transcript, canonical, mane_select)`
tuple. Produced by `--aggregate-per-variant-expected`.

- **Path:** `gs://gnomad/v{version}/constraint/apply_models/transcript_consequences/gnomad.v{version}.per_variant_expected.aggregated.ht`
- **Resource fn:** [`get_aggregated_per_variant_expected`](resource_utils.py)
- **Key:** `annotation, modifier, gene, gene_id, transcript, canonical, mane_select`

### Globals
Identical pipeline-parameter globals as `per_variant_expected.ht`.

### Rows
- `annotation, modifier, gene, gene_id, transcript, canonical, mane_select` — key.
- `mu_snp: float64` — summed `mu_snp * possible_variants` (i.e. total per-context
  mu mass over the bucket; despite the name, this is no longer a per-SNV rate).
- `mu: float64` — summed per-variant `mu` (`mu_snp * possible *
  coverage_correction` over the bucket).
- `observed_variants: array<int64>` — summed observed counts per freq entry.
- `possible_variants: int64` — total possible SNV count in the bucket.
- `predicted_proportion_observed: array<float64>` — `Σ predicted_proportion_observed
  · possible` (note: weighted by `possible`, divide by `possible_variants` to
  recover the bucket-mean fraction).
- `coverage_correction: float64` — `Σ coverage_correction · possible`.
- `expected_variants: array<float64>` — summed expected counts per freq entry.

---

## 8. `aggregated_expected.ht`

Alternative path: aggregate *first*, then apply models. Produced by
`--apply-models-aggregated` and has the same schema as
`per_variant_expected.aggregated.ht`.

- **Path:** `gs://gnomad/v{version}/constraint/apply_models/transcript_consequences/gnomad.v{version}.aggregated_expected.ht`
- **Resource fn:** [`get_aggregated_expected`](resource_utils.py)
- **Key:** identical to §7.

### Globals
Same as §7, except `apply_models_globals.groupings` is `array<str>` (rather
than `tuple(str×7)`) because aggregation happens before per-variant grouping.

### Rows
Identical to §7. Use this table to compare aggregate-before-apply vs.
apply-then-aggregate (§7) results; the canonical downstream input is §7.

---

## 9. `constraint_group.ht`

Pre-metric per-transcript table organized into constraint groups (`syn`, `mis`,
`lof_classic`, `lof_hc_lc`, `lof_hc`, plus any additional groupings). Produced by
`--aggregate-by-constraint-groups` /
[`aggregate_by_constraint_groups`](../utils/constraint.py).

- **Path:** `gs://gnomad/v{version}/constraint/apply_models/transcript_consequences/gnomad.v{version}.constraint_group.ht`
- **Resource fn:** [`get_constraint_group_ht`](resource_utils.py)
- **Key:** `gene, gene_id, transcript, canonical, mane_select`

### Globals
Pipeline-parameter globals (same as §7) plus:
- `constraint_group_meta: array<dict<str,str>>` — one dict per
  `constraint_groups[i]` describing the group's filter: `{csq_set: syn}`,
  `{csq_set: mis}`, `{lof: classic}` (stop gained, splice donor, and splice
  acceptor variants), `{lof: hc_lc}` (LOFTEE HC or LC), and `{lof: hc}`
  (LOFTEE HC). Position in this array matches position in the row-level
  `constraint_groups` array.

### Rows
- `gene, gene_id, transcript, canonical, mane_select` — key.
- `constraint_groups: array<struct>` — one element per metadata entry in
  `constraint_group_meta`, each carrying summed counts:
  - `mu_snp: float64` — `Σ mu_snp · possible` over variants in this group.
  - `mu: float64` — `Σ mu` (i.e. `Σ mu_snp · possible · coverage_correction`).
  - `possible_variants: int64` — total possible SNVs.
  - `coverage_correction: float64` — `Σ coverage_correction · possible`.
  - `oe_info: array<struct{observed_variants: int64, predicted_proportion_observed: float64, expected_variants: float64}>`
    — one entry per freq-meta index. `observed_variants` is the summed
    observed count, `expected_variants` is the summed predicted expected
    count, and `predicted_proportion_observed` is the `Σ ppo·possible` sum.
- `no_variants: bool` — true if the `syn`, `mis`, and `lof_hc` groups have zero
  observed variants at the global-adj index. Rows with no observed and no
  expected variants in any group are dropped.

---

## 10. `constraint_metrics.ht` (internal)

Per-transcript metrics: pLI (LoF groups), OE confidence intervals (two estimators),
raw and standardized z-scores, percentile/decile/sextile bins. Produced by
`--compute-constraint-metrics`. This is the *internal* metrics table; the
public release table (§12) is a renamed/flattened projection.

- **Path:** `gs://gnomad/v{version}/constraint/metrics/transcript_consequences/gnomad.v{version}.constraint_metrics.ht`
- **Resource fn:** [`get_constraint_metrics_dataset`](resource_utils.py)
- **Key:** `gene, gene_id, transcript, canonical, mane_select`

### Globals
Pipeline-parameter globals (as §9) plus:
- `constraint_group_meta` — see §9.
- `sd_raw_z: array<float64>` — per-constraint-group standard deviation of the
  raw z-statistic across all transcripts (parallel to `constraint_group_meta`).
  Used to scale `z_raw → z_score`.
- `percentile_thresholds: struct` — observed quantile cutoffs of the
  `oe_ci.upper` distribution. One sub-struct per metric:
  - `syn`, `mis`, `lof: struct{percentile: array<float64>, decile: array<float64>, sextile: array<float64>}`
    — each granularity's boundary values (e.g.
    `percentile_thresholds.lof.decile[0]` = 10th-percentile of LOEUF upper).

### Rows
- `gene, gene_id, transcript, canonical, mane_select` — key.
- `constraint_groups: array<struct>` — one element per group (parallel to
  `constraint_group_meta`). Each element extends §9's struct:
  - `mu_snp`, `mu`, `possible_variants`, `coverage_correction` — as in §9.
  - `oe_info: array<struct>` — extended with metric fields:
    - `observed_variants: int64`, `predicted_proportion_observed: float64`,
      `expected_variants: float64`
    - `oe: float64` — `observed / expected` ratio.
    - `oe_ci_discretized_poisson: struct{lower: float64, upper: float64}` — OE
      CI using a discretized-Poisson estimator.
    - `oe_ci_gamma: struct{lower: float64, upper: float64}` — OE CI using a
      Gamma estimator (the LOEUF estimator).
    - `z_raw: float64` — `(expected - observed) / sqrt(expected)`, so it is
      positive when there are fewer variants than expected.
  - `flags: set<str>` — per-group QC flags: `no_exp_{syn,mis,lof}` when the
    expected count is missing or 0, and `outlier_{syn,mis,lof}` when `z_raw`
    is outside the group's outlier thresholds. Every LoF group uses the `lof`
    name and thresholds.
  - `z_score: float64` — `z_raw / sd_raw_z[group_idx]`, taken at the
    canonical (global-adj) freq index.
  - `pNull: float64`, `pRec: float64`, `pLI: float64` — pNull / pRec / pLI
    scores computed from this group's global-adj observed / expected via
    [`compute_pli`](https://github.com/broadinstitute/gnomad_methods). Defined
    only for the groups in `RELEASE_GROUPS_WITH_PLI` (`lof_hc`, `lof_hc_lc`);
    missing for every other group.
  - `oe_ci_discretized_poisson_rank: struct` and `oe_ci_gamma_rank: struct` —
    rank/bin annotations for the `upper` bound of each CI:
    - `upper_rank: int64` — 0-based ascending rank within the group.
    - `upper_bin_percentile: int32` — percentile bin (0–99).
    - `upper_bin_decile: int32` — decile bin (0–9).
    - `upper_bin_sextile: int32` — sextile bin (0–5).
- `no_variants: bool` — as §9.
- `constraint_flags: set<str>` — transcript-level QC flags: the union of the
  `syn`, `mis`, and `lof_hc` groups' `flags` (the `no_exp_*` and `outlier_*`
  flags), plus `no_variants` when `no_variants` is true. Flagged transcripts are
  left out of the percentile thresholds.
- `constraint_bins: struct{percentile, decile, sextile}` — convenience bins
  per metric, where each is `struct{syn: int32, mis: int32, lof: int32}`.
- `gene_quality_metrics: struct` — joined-in from §11:
  - `exome_prop_bp_AN90: float64` — fraction of transcript CDS bp with
    exome AN ≥ 90% of max.
  - `exome_mean_AS_MQ: float64` — mean AS_MQ over the transcript.
  - `exome_prop_segdup: float64` — fraction of transcript CDS in segmental
    duplications.
  - `exome_prop_LCR: float64` — fraction in low-complexity regions.
- `gene_flags: set<str>` — joined-in from §11 (`low_exome_coverage`,
  `low_exome_mapping_quality`).
- `level: str`, `transcript_type: str`, `chromosome: str`,
  `start_position: int32`, `end_position: int32`,
  `gene_id_version: str`, `transcript_id_version: str`,
  `cds_length: int64`, `num_coding_exons: int64` — GENCODE annotations.

---

## 11. `gene_quality_metrics.ht`

Per-transcript coverage / mapping-quality and region-overlap metrics. Produced
by `--compute-gene-quality-metrics`.

- **Path:** `gs://gnomad/v{version}/constraint/metrics/gnomad.v{version}.gene_quality_metrics.ht`
- **Resource fn:** [`get_gene_quality_metrics_ht`](resource_utils.py)
- **Key:** `transcript`

### Globals
*(none)*

### Rows
- `transcript: str` — key.
- `gene_quality_metrics: struct{exome_prop_bp_AN90, exome_mean_AS_MQ, exome_prop_segdup, exome_prop_LCR}`
  — see §10 for the individual fields.
- `gene_flags: set<str>` — transcript-level flags derived from the above:
  `low_exome_coverage` when `exome_prop_bp_AN90` < 0.1, and
  `low_exome_mapping_quality` when `exome_mean_AS_MQ` < 50.

---

## 12. `release/constraint_metrics.ht`

Public-release flattened constraint table. One row per
`(gene, gene_id, transcript, canonical, mane_select)`. Produced by
`--prepare-release` from §10 by renaming/selecting fields per the
`RELEASE_*` constants in [constants.py](constants.py).

- **Path:** `gs://gnomad/v{version}/constraint/release/gnomad.v{version}.constraint_metrics.ht`
- **Resource fn:** [`get_release_constraint_ht`](resource_utils.py)
- **Key:** `gene, gene_id, transcript, canonical, mane_select`
- A flat **TSV** (`constraint_metrics.tsv.bgz`), a **downsampling TSV**
  (`constraint_metrics.downsampling.tsv.bgz`), and a **LOEUF percentile
  thresholds TSV** (`loeuf_percentile_thresholds.tsv`) are exported alongside.

### Globals
- `version: str` — gnomAD release version (e.g. `4.1.1`).
- `calculate_mu_params: struct` — release-cleaned subset of
  `calculate_mu_globals` (drops `freq_meta`, `genetic_ancestry_groups`,
  `downsampling_idx`):
  - `ac_cutoff, min_cov, max_cov, gerp_lower_cutoff, gerp_upper_cutoff,
    downsampling_level, most_severe_consequence`.
- `build_models_params: struct{low_cov_cutoff, high_cov_cutoff, upper_cov_cutoff}`
  — release-cleaned `build_models_globals` (drops
  `synonymous_transcript_filter_field`, `skip_coverage_model`).
- `apply_models_params: struct{low_cov_cutoff, high_cov_cutoff, plateau_models, coverage_model, log10_coverage}`
  — release-cleaned `apply_models_globals` (drops `skip_coverage_model`,
  `groupings`).
- `downsamplings: struct{global, afr, amr, eas, nfe, sas}` — each is
  `array<int32>` of downsampling sizes present for that gen-anc group
  (parallel to the per-row `gen_anc_obs.<group>` / `gen_anc_exp.<group>`
  arrays).
- `max_af: float64` — AF cap for observed.
- `sd_raw_z: struct{syn, mis, lof_hc_lc, lof}` — per-release-group `sd_raw_z`
  (re-keyed from positional array to named struct via `RELEASE_GROUP_RENAMES`,
  which maps internal `lof_hc → lof` for release).
- `loeuf_percentile_thresholds: struct{percentile, decile, sextile}` — LOEUF
  upper-bound thresholds at each granularity (copied from
  `percentile_thresholds.lof` in §10).

### Rows
Key + GENCODE annotations + per-group structs:

- `gene, gene_id, transcript, canonical, mane_select` — key.
- `transcript_version: str` (renamed from `transcript_id_version`),
  `transcript_type: str`,
  `transcript_level: str` (renamed from `level`),
  `chromosome: str`, `start_position: int32`, `end_position: int32`,
  `cds_length: int64`, `num_coding_exons: int64`.
- `gene_quality_metrics: struct{exome_prop_bp_AN90, exome_mean_AS_MQ, exome_prop_segdup, exome_prop_LCR}`
  — see §10.
- `gene_flags: set<str>`, `constraint_flags: set<str>` — propagated from §10/§11.

Each constraint group is exposed as a named struct (group names from
`RELEASE_GROUP_NAMES`):

- `syn: struct` — synonymous constraint (fields below).
- `mis: struct` — missense constraint.
- `lof_hc_lc: struct` — pLoF (LOFTEE HC + LC).
- `lof: struct` — pLoF (LOFTEE HC only; renamed from internal `lof_hc`). This
  is the canonical LOEUF group.

**Common sub-fields** on every `{syn,mis,lof_hc_lc,lof}` struct
(`RELEASE_CG_SELECT` order):

- `mu: float64` — total mu mass for this group, without the coverage
  correction: the internal `mu_snp` sum, renamed per `RELEASE_CG_RENAME`. It
  replaces the internal `mu`, which includes the coverage correction.
- `possible: int64` — total possible SNVs (from `possible_variants`).
- `obs: int64` — observed SNV count at the global-adj index.
- `exp: float64` — expected SNV count at the global-adj index.
- `oe: float64` — `obs / exp`.
- `z_raw: float64` — raw z-statistic.
- `z_score: float64` — `z_raw / sd_raw_z[<group>]`.
- `oe_ci: struct{lower: float64, upper: float64}` — OE confidence interval.
  Every group uses the Gamma CI (`oe_ci_gamma` in §10; for `lof` this is the
  LOEUF). `lof` **additionally** carries rank fields on `oe_ci` (the only group
  in `RELEASE_GROUPS_WITH_RANK`):
    - `upper_rank: int64` — 0-based rank by `oe_ci.upper` ascending.
    - `upper_bin_percentile: int32` — percentile bin (0–99).
    - `upper_bin_decile: int32` — decile bin (0–9; LOEUF decile).
    - `upper_bin_sextile: int32` — sextile bin (0–5).
- `gen_anc_obs: struct{global, afr, amr, eas, nfe, sas}` — per-genetic-ancestry
  observed counts. Each field is `array<int64>` parallel to
  `downsamplings.<group>`.
- `gen_anc_exp: struct{global, afr, amr, eas, nfe, sas}` — per-genetic-ancestry
  expected counts. Each field is `array<float64>` parallel to
  `downsamplings.<group>`.

**Additional sub-fields on `lof_hc_lc` and `lof`** (`RELEASE_GROUPS_WITH_PLI`):
- `pLI: float64`, `pNull: float64`, `pRec: float64` — probability of being
  loss-of-function intolerant / null / recessive (Lek et al. 2016).

---

## 13. `release/mutation_rate.ht`

Public-release mutation rate. Produced by `--prepare-release-mutation-rate`
from §3 by selecting/renaming.

- **Path:** `gs://gnomad/v{version}/constraint/release/gnomad.v{version}.mutation_rate.ht`
- **Resource fn:** [`get_release_mutation_ht`](resource_utils.py)
- **Key:** `context, ref, alt, methylation_level`
- A TSV mirror is exported to `release/gnomad.v{version}.mutation_rate.tsv`.

### Globals
- `version: str`.
- `calculate_mu_params: struct{ac_cutoff, min_cov, max_cov, gerp_lower_cutoff,
  gerp_upper_cutoff, downsampling_level, most_severe_consequence}` —
  release-cleaned `calculate_mu_globals`.

### Rows
- `context, ref, alt, methylation_level` — key.
- `mu: float64` — scalar mu at the canonical downsampling index (renamed from
  `mu_snp`).
- `cpg: bool`, `transition: bool`, `mutation_type: str` — mutation-type tags.

---

# Regional mutation rate tables

Tables used to compute the regional mutation rate (rr): the outputs of
[import_data.py](../pipeline/import_data.py) (§14-15) and
[regional_rate_pipeline.py](../pipeline/regional_rate_pipeline.py) (§16-19). Only
the penalty tuning, model, and rr Tables (§17-19) have globals.

The genomic feature columns are named `{feature}_{scale}` for the 13 features in
`GENOMIC_FEATURES` (`dist2telo, dist2cent, GC_content, LCR, SINE, LINE,
recomb_male, recomb_female, met_sperm, Nucleosome, cDNM_maternal_05M,
cDNM_paternal_05M, CpG_island`) at the 4 flanking scales in
`GENOMIC_FEATURE_SCALES` (`1k, 10k, 100k, 1M`); see [constants.py](constants.py).

---

## 14. `imported_data/dnm.ht`

De novo mutation (DNM) sites from the Satterstrom 2026 ASD and DDD files
([`get_dnm_paths`](resource_utils.py)), produced by `--import-dnms` /
[`format_dnm_ht`](../utils/import_data.py). Each variant has one row, however many
files or individuals report it. Indels and alleles that don't match the reference
are kept; they are removed when the DNMs are joined to `annotated_context.ht`.

- **Path:** `gs://gnomad/v{version}/constraint/imported_data/gnomad.v{version}.dnm.ht`
- **Resource fn:** [`get_dnm_ht`](resource_utils.py)
- **Key:** `locus`, `alleles`

### Rows
- `locus: locus<GRCh38>`, `alleles: array<str>` — key; no other fields.

---

## 15. `imported_data/genomic_features13_1kb.ht`

The published Gnocchi 1kb regional genomic feature matrix
([`get_published_genomic_features_path`](resource_utils.py)), produced by
`--import-genomic-features` /
[`format_genomic_features`](../utils/import_data.py). It covers chr1-22 only.

- **Path:** `gs://gnomad/v{version}/constraint/imported_data/gnomad.v{version}.genomic_features13_1kb.ht`
- **Resource fn:** [`get_genomic_features_ht`](resource_utils.py)
- **Key:** `interval`

### Rows
- `interval: interval<locus<GRCh38>>` — the BED window in `element_id` in 1-based
  positions, `[start + 1, end]`. (The Gnocchi rr Tables used `[start, end)`, 1bp
  left of the BED window.)
- `element_id: str` — window ID from the matrix (`chr-start-end`, BED coordinates).
- `{feature}_{scale}: float64` — the 52 feature columns. Values written as `NA`,
  `nan`, `NaN`, or empty in the matrix are missing.

---

## 16. `regional_rate/rr_training.ht`

rr training set, produced by `--create-rr-training-set` /
[`create_rr_training_ht`](../utils/regional_rate.py). The DNM sites are the
synonymous sites of `annotated_context.ht` in `dnm.ht`; the background is a
`RR_BACKGROUND_FRACTION` (5%) sample, with seed `RR_BACKGROUND_SEED`, of the other
synonymous sites that are absent from gnomAD exomes (`filters.exomes` is missing),
one per locus. Sites in the span of an Ensembl MANE Select transcript with a
`gene_flags` entry in the v4.1.2 constraint metrics Table are removed
([`remove_qc_flagged_genes`](../utils/regional_rate.py)). Each DNM site is one row,
however many DNMs were reported there.

- **Path:** `gs://gnomad/v{version}/constraint/regional_rate/gnomad.v{version}.rr_training.ht`
- **Resource fn:** [`get_rr_training_ht`](resource_utils.py)
- **Key:** `locus`

### Rows
- `locus: locus<GRCh38>` — key.
- `alleles: array<str>`, `context: str` — the site's alleles and strand-collapsed
  context.
- `y: int32` — 1 for a DNM site, 0 for a background site.
- `{feature}_{scale}: float64` — the 52 feature columns of the site's 1kb window
  (missing outside a window).

---

## 17. `regional_rate/rr_penalty_tuning.ht`

Held-out log loss of the rr models for each ridge penalty tried, produced by
`--tune-rr-penalty` with `tune_regional_mutation_rate_penalty` from
`gnomad.utils.constraint.regional_mutation_rate`. For each penalty, context, and
fold, the context's model is fit on the training sites of the other fold's
chromosomes (odd or even) and scored on the fold's sites.

- **Path:** `gs://gnomad/v{version}/constraint/regional_rate/gnomad.v{version}.rr_penalty_tuning.ht`
- **Resource fn:** [`get_rr_penalty_tuning_ht`](resource_utils.py)
- **Key:** `penalty`, `context`, `fold`

### Globals
- `ridge_penalty: float64` — the penalty with the lowest held-out log loss summed
  over contexts and folds, or `RR_RIDGE_PENALTY` (0.3) if no held-out sites could be
  scored (a `--test` run covers chr20 only); `--train-rr-models` uses it for every
  context.

### Rows
- `penalty: float64`, `context: str`, `fold: str` — key (`fold` is the held-out
  chromosome set, `odd` or `even`).
- `n_test: int64` — held-out training sites with no missing feature.
- `log_loss: float64` — log loss summed over those sites.

---

## 18. `regional_rate/rr_models.ht`

Per-context ridge-penalized logistic rr models, produced by `--train-rr-models` with
`fit_regional_mutation_rate_models` from
`gnomad.utils.constraint.regional_mutation_rate`. The models are fit on z-scored
features, so a window's linear predictor is `intercept + sum(coefs * (x - mean) / sd)`.

- **Path:** `gs://gnomad/v{version}/constraint/regional_rate/gnomad.v{version}.rr_models.ht`
- **Resource fn:** [`get_rr_models_ht`](resource_utils.py)
- **Key:** `context`

### Globals
- `ridge_penalty: float64` — per-site ridge penalty the models were fit with, from
  §17.
- `non_converged_contexts: array<str>` — contexts whose logistic fit didn't
  converge (also logged as a warning).

### Rows
- `context: str` — strand-collapsed trinucleotide context (one of `RR_CONTEXTS`).
  Contexts without DNM or background training sites, or where no feature varies,
  have no row.
- `features: array<str>` — model features in coefficient order: the context's
  `RR_CONTEXT_FEATURES`, without any feature that has a single value in the
  context's training sites.
- `mean: array<float64>`, `sd: array<float64>` — training mean and standard
  deviation (ddof=0) of each feature.
- `coefs: array<float64>` — coefficient of each standardized feature.
- `intercept: float64` — intercept of the linear predictor.
- `ave: float64` — predicted probability for the average window (all features at
  their mean); rr divides by it. As in Gnocchi, rr is this ratio of probabilities
  rather than the odds ratio, which compresses the regional effect toward 1 by
  about `1 - ave` on the log scale and caps rr at `1 / ave`. `ave` is about the DNM
  share of the training sites, a few percent for non-CpG contexts but 20-50% for
  CpG contexts, so CpG rr is damped.
- `n_mutated: int32`, `n_background: int32` — DNM and background training sites
  with no missing feature.
- `converged: bool` — whether the logistic fit converged.

---

## 19. `regional_rate/syn_dnm_rr_1kb.ht`

Normalized rr per 1kb window and context, produced by `--compute-rr` with
`compute_regional_mutation_rate_ht` and `normalize_regional_mutation_rate_expr` from
`gnomad.utils.constraint.regional_mutation_rate`. Read by
`constraint_pipeline.py --aggregate-per-variant-expected --apply-regional-rr`.

- **Path:** `gs://gnomad/v{version}/constraint/regional_rate/gnomad.v{version}.syn_dnm_rr_1kb.ht`
- **Resource fn:** [`get_regional_rr_ht`](resource_utils.py)
- **Key:** `interval`

### Globals
- `rr_globals: struct`:
  - `ridge_penalty: float64` — ridge penalty the models were fit with (§17).
  - `normalization_factors: dict<str, float64>` — per-context expected-weighted
    mean of the unnormalized rr over the exome synonymous Ensembl MANE Select sites
    of `per_variant_expected.ht` that fall in an rr window; rr is divided by it.
  - `per_variant_apply_models_globals: struct` — the `apply_models_globals` of the
    per-variant expected Table rr was normalized on (see §7).

### Rows
- `interval: interval<locus<GRCh38>>` — 1kb window, as in
  `genomic_features13_1kb.ht`. Windows where no modeled context has an rr are
  removed.
- `adj_r: dict<str, float64>` — rr for each of the 32 `RR_CONTEXTS`:
  `sigmoid(linear predictor) / ave / normalization factor`. Before normalization,
  rr is set to 1.0 in windows missing one of the context's model features, so those
  windows get `1 / normalization factor`; contexts without a model are 1.0.

---

# Constraint-specific input tables

These are constraint-pipeline-owned input tables (not part of `gnomad_methods`)
that are read by `prepare_ht_for_constraint_calculations` to annotate
`annotated_context.ht`.

## `adj_r_per_context_methyl_genome_1kb_autosome.agg.ht`

Per-context regional-depletion correction. Aggregated to 1kb autosomal
intervals over the genome reference; values vary by trinucleotide context.

- **Path:** `gs://gnomad/v4.1/constraint/resources/annotations/ht/adj_r_per_context_methyl_genome_1kb_autosome.agg.ht`
- **Resource fn:** [`get_adj_r_ht`](resource_utils.py)
- **Key:** `interval`

### Globals
*(none)*

### Rows
- `interval: interval<locus<GRCh38>>` — 1kb autosomal interval.
- `adj_r: dict<str, float64>` — context-keyed correction. Look up with
  `adj_r_ht[ht.locus].adj_r[ht.context]`.

## `adj_r_syn_dnm_per_context_methyl_genome_1kb_autosome.agg.ht`

Identical schema to `adj_r`, but computed from a synonymous de-novo-mutation
constraint baseline rather than the general per-context baseline.

- **Path:** `gs://gnomad/v4.1/constraint/resources/annotations/ht/adj_r_syn_dnm_per_context_methyl_genome_1kb_autosome.agg.ht`
- **Resource fn:** [`get_syn_adj_r_ht`](resource_utils.py)
- **Key:** `interval`

### Globals
*(none)*

### Rows
- `interval: interval<locus<GRCh38>>`.
- `adj_r: dict<str, float64>` — context-keyed correction. Joined into
  `annotated_context.syn_adj_r` as `adj_r_ht[ht.locus].adj_r[ht.context]`.

---

## Updating this README

After a change to an output's schema, describe the output of a `--test` run and
update its section. In a Hail session with access to GCS:

```python
import hail as hl

from gnomad_constraint.resources.resource_utils import get_constraint_metrics_dataset

ht = get_constraint_metrics_dataset(
    "transcript_consequences", version="4.1.1", test=True
).ht()
ht.describe()
print(hl.eval(ht.globals))
```

The path conventions are sourced from [`get_constraint_data`](resource_utils.py)
and [`get_constraint_root`](resource_utils.py). When in doubt, instantiate a
resource via the corresponding `get_*` function and read its `.path`.
