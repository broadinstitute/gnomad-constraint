# gnomad-constraint

This repo contains the code that generates the gnomAD gene and transcript
constraint metrics: observed and expected variant counts, the observed/expected
(o/e) ratio and its confidence interval, pLI, and z-scores for synonymous,
missense, and predicted loss-of-function (pLoF) variants, plus the mutation rate
table they are built on. Expected counts are computed per variant and then
aggregated by transcript and constraint group.

The current version is v4.1.1 (GRCh38). Version 2.1.1 (GRCh37) is no longer
supported; commit 39928d1 is the last version of the pipeline that supports it.

An overview of the pipeline and the functions it uses is in
[gnomad_constraint/flowchart/constraint_pipeline_v4.pdf](./gnomad_constraint/flowchart/constraint_pipeline_v4.pdf).
Many of the functions are imported from
[gnomad_methods](https://github.com/broadinstitute/gnomad_methods).

## Repository layout

| Path | Contents |
|---|---|
| `gnomad_constraint/pipeline/` | The pipeline scripts, run in order: `prepare_data.py`, `expected_models.py`, `expected_counts.py`, and `constraint_metrics.py`; each step is a command-line flag |
| `gnomad_constraint/utils/constraint.py` | Functions used by the pipeline steps |
| `gnomad_constraint/utils/import_data.py` | Functions used by the import steps |
| `gnomad_constraint/resources/resource_utils.py` | Input and output resources (GCS paths) |
| `gnomad_constraint/resources/constants.py` | Pipeline and release constants |
| `gnomad_constraint/resources/README.md` | Schema reference for every pipeline output |
| `gnomad_constraint/plots/` | R scripts for constraint figures; see its [README](gnomad_constraint/plots/README.md) for setup |
| `gnomad_constraint/flowchart/` | Pipeline overview |
| `tests/` | Unit tests; see its [README](tests/README.md) |
| `release_readme/` | Field descriptions for the released constraint and mutation rate Tables |
| `requirements.txt` | What the pipeline needs to run: gnomad_methods, gnomad_qc, and Hail |
| `requirements-dev.in`, `requirements-dev.txt` | Development tools; the `.txt` is generated from the `.in` |
| `.claude/` | Shared Claude Code configuration (`settings.json`): a hook that clones [gnomad-skills](https://github.com/broadinstitute/gnomad-skills) into `~/.claude/gnomad-skills` (or `$GNOMAD_SKILLS_DIR`), links its skills into `~/.claude/skills/`, and pulls it at most once a day (set `GNOMAD_SKIP_SKILLS_HOOK=1` to turn it off), and `attribution` settings that put `Assisted by Claude Code.` in PR descriptions and nothing on commits. |

## Pipeline steps

The pipeline is four scripts in `gnomad_constraint/pipeline/`, run in this order.
Each runs one or more of its steps per invocation, each enabled by a flag:

| Script | Flag | Step |
|---|---|---|
| `prepare_data.py` | `--import-dnms`, `--import-genomic-features` | Import the de novo mutations and genomic features the regional mutation rate is modeled from; see [Regional mutation rate (rr)](#regional-mutation-rate-rr) |
| | `--prepare-context-ht` | Annotate the VEP context Table with methylation, coverage, allele number, frequency, filters, GERP, and genomic region |
| | `--calculate-gerp-cutoffs` | Optional: compute GERP percentile cutoffs. Nothing is saved, so it must run in the same command as `--preprocess-data`, which uses them in place of `--calculate-mutation-rate-gerp-{lower,upper}-cutoff` |
| | `--preprocess-data` | Add per-site observed and possible counts, mutation rate inputs, and model groups |
| | `--compute-gene-quality-metrics` | Per-transcript coverage, mapping quality, segmental duplication, and low-complexity region metrics |
| `expected_models.py` | `--calculate-mutation-rate` | Compute the mutation rate for each substitution and context |
| | `--create-training-set` | Count observed and possible synonymous variants for model training |
| | `--build-models` | Build the plateau and coverage models |
| `expected_counts.py` | `--apply-models-per-variant` | Apply the models to compute expected counts per variant |
| | `--apply-models-aggregated` | Alternative to `--apply-models-per-variant` and `--aggregate-per-variant-expected`: sum observed and possible counts first, then apply the models |
| `constraint_metrics.py` | `--aggregate-per-variant-expected` | Sum the per-variant counts by transcript and consequence |
| | `--aggregate-by-constraint-groups` | Sum by constraint group (synonymous, missense, pLoF). Add `--use-aggregated-expected` to read the `--apply-models-aggregated` output |
| | `--compute-constraint-metrics` | Compute o/e, confidence intervals, z-scores, flags, pLI, ranks, and bins. `--skip-pre-rank-metrics` reruns only the ranking |
| | `--prepare-release` | Format the constraint metrics as the release Table |
| | `--prepare-release-mutation-rate` | Format the mutation rate Table for release |
| | `--export-release-tsv`, `--export-release-downsampling-tsv` | Export the release Table and its per-ancestry downsampling counts to TSV |

`--version` selects the gnomAD version (default 4.1.1). Autosomes, chrX, and chrY
are all computed; each site is labeled with a `genomic_region` of
`autosome_or_par`, `chrx_nonpar`, or `chry_nonpar`. Run a script with `--help`
for its options.

Outputs are written under `gs://gnomad/v{version}/constraint/`.
[gnomad_constraint/resources/README.md](gnomad_constraint/resources/README.md)
has the full layout and the schema of every output. With `--test` (chr20, chrX,
and chrY only) or `--test-gene-list` (six genes), every output, including the
release files, goes to `gs://gnomad-tmp/gnomad_v{version}_testing/constraint/`
instead, so a test run never touches production outputs. `--directory-post-fix`
and `--path-post-fix` keep separate runs of the same version apart. The directory
post-fix applies to every output; the path post-fix only to the outputs from
`--create-training-set` onward, except the release mutation rate files.

### Regional mutation rate (rr)

The regional mutation rate (rr) is a factor for each 1kb window and trinucleotide
context, modeled as in Gnocchi from synonymous de novo mutations (DNMs) and 13
regional genomic features of each site's window at four flanking scales.

`prepare_data.py` imports the DNMs and genomic features rr is modeled from, writing
under `imported_data/` in the constraint output directory (with `--test`, in the
testing bucket, and with neither post-fix):

| Flag | Step |
|---|---|
| `--import-dnms` | Import the Satterstrom 2026 ASD and DDD DNMs as a Table of DNM sites |
| `--import-genomic-features` | Import the Gnocchi 1kb regional genomic feature matrix (chr1-22, so chrX and chrY get no rr) |

### Release constraint metrics Table

Keyed by `(gene, gene_id, transcript, canonical, mane_select)`, with one struct
per constraint group (`syn`, `mis`, `lof_hc_lc`, `lof`) holding the observed and
expected counts, o/e and its 90% gamma confidence interval (for `lof`, the LOEUF,
with rank and bin fields), raw and standardized z-scores, and, for `lof_hc_lc` and
`lof`, pLI, pNull, and pRec computed from that group's own counts.
[release_readme/gene_constraint_ht_README.md](release_readme/gene_constraint_ht_README.md)
describes every field.

## Setup

Python 3.10 or later and a Java runtime for Hail:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pip install pre-commit && pre-commit install
python -m pytest
```

The pipeline needs gnomad_methods at or after the per-base constraint changes
(gnomad_methods #824, not yet on PyPI), installed from GitHub at the commit
pinned in `requirements.txt`, and Hail 0.2.137 or later for the gamma
confidence interval.

To change a development dependency, edit `requirements-dev.in` and regenerate
the `.txt` with Python 3.11 (the version CI runs):

```bash
pip-compile requirements-dev.in
```

## Running the pipeline

Steps run on a Dataproc cluster. Package `gnomad_constraint`, `gnomad_qc`, and
`gnomad_methods` (at the commit pinned in `requirements.txt`) into a single zip:
if you pass hailctl more than one file, it repackages them and nests the
packages incorrectly.

```bash
rm -f /tmp/pyfiles.zip
(cd <gnomad-constraint checkout> && zip -r /tmp/pyfiles.zip gnomad_constraint \
  -x '*.pyc' '*__pycache__*' '*.RData' '*/renv/*' '*.DS_Store' '*.png' '*.pdf')
(cd <gnomad_qc checkout> && zip -r /tmp/pyfiles.zip gnomad_qc -x '*.pyc' '*__pycache__*')
(cd <gnomad_methods checkout at the pinned commit> && zip -r /tmp/pyfiles.zip gnomad -x '*.pyc' '*__pycache__*')

# From the gnomad-constraint checkout:
hailctl dataproc submit <cluster> gnomad_constraint/pipeline/constraint_metrics.py \
  --pyfiles /tmp/pyfiles.zip -- --compute-constraint-metrics --test
```

Run a step with `--test` and check its output before running it at full scale.

## Continuous integration

CI (`.github/workflows/ci.yml`) runs on pull requests and pushes to `main`. It
checks formatting (black, isort, autopep8), docstrings (pydocstyle), and pylint
errors on the Python code, runs the unit tests, and runs styler and lintr on the
R code. The same formatters run as pre-commit hooks.
