# gnomad-constraint

This repo contains the code that generates the gnomAD gene and transcript
constraint metrics: observed and expected variant counts, the observed/expected
(o/e) ratio and its confidence interval, pLI, and z-scores for synonymous,
missense, and predicted loss-of-function (pLoF) variants, plus the mutation rate
table they are built on.

An overview of the pipeline and the functions it uses is in
[gnomad_constraint/flowchart/constraint_pipeline_v4.pdf](./gnomad_constraint/flowchart/constraint_pipeline_v4.pdf).
Many of the functions are imported from
[gnomad_methods](https://github.com/broadinstitute/gnomad_methods).

## Repository layout

| Path | Contents |
|---|---|
| `gnomad_constraint/pipeline/constraint_pipeline.py` | The pipeline script; each step is a command-line flag |
| `gnomad_constraint/utils/constraint.py` | Functions used by the pipeline steps |
| `gnomad_constraint/resources/resource_utils.py` | Input and output resources (GCS paths) and pipeline constants |
| `gnomad_constraint/plots/` | R scripts for constraint figures; see its [README](gnomad_constraint/plots/README.md) for setup |
| `gnomad_constraint/flowchart/` | Pipeline overview |
| `release_readme/` | Field descriptions for the released constraint and mutation rate Tables |
| `requirements.txt` | What the pipeline needs to run: gnomad_methods, gnomad_qc, and Hail |
| `requirements-dev.in`, `requirements-dev.txt` | Development tools; the `.txt` is generated from the `.in` |
| `.claude/` | Shared Claude Code configuration (`settings.json`): a hook that clones [gnomad-skills](https://github.com/broadinstitute/gnomad-skills) into `~/.claude/gnomad-skills` (or `$GNOMAD_SKILLS_DIR`), links its skills into `~/.claude/skills/`, and pulls it at most once a day (set `GNOMAD_SKIP_SKILLS_HOOK=1` to turn it off), and `attribution` settings that put `Assisted by Claude Code.` in PR descriptions and nothing on commits. |

## Pipeline steps

`constraint_pipeline.py` runs one or more steps per invocation, each enabled by a
flag, in this order:

| Flag | Step |
|---|---|
| `--prepare-context-ht` | Split multiallelic sites in the VEP context Table and add methylation, coverage, and GERP annotations |
| `--preprocess-data` | Prepare the exome, genome, and context Tables for constraint calculations, one Table per genomic region |
| `--calculate-gerp-cutoffs` | Optional: compute GERP percentile cutoffs to use in place of `--gerp-lower-cutoff` and `--gerp-upper-cutoff`. Nothing is saved, so run it in the same command as `--calculate-mutation-rate`, the only step that uses them |
| `--calculate-mutation-rate` | Compute the baseline mutation rate for each substitution and context from genome downsamplings |
| `--create-training-set` | Count observed and possible variants by exome coverage at synonymous sites |
| `--build-models` | Build the plateau and coverage models |
| `--apply-models` | Apply the models to compute expected variant counts |
| `--compute-constraint-metrics` | Compute o/e, its confidence interval, pLI, and z-scores |
| `--export-tsv` | Export the constraint metrics to a TSV |

`--version` selects the gnomAD version (default 4.1). For v4, only autosomes are
computed: chrX and chrY are skipped, and the `autosome_par` Tables leave out the
pseudoautosomal regions until there is X/Y methylation data (see the TODO in
`--preprocess-data`). Run `constraint_pipeline.py --help` for every option.

Outputs are written under `gs://gnomad/v{version}/constraint/`. With `--test`
(chr20, chrX, and chrY only) or `--test-gene-list` (five genes), they go to
`gs://gnomad-tmp/gnomad_v{version}_testing/constraint/` instead, so a test run
never touches production outputs.

## Setup

Python 3.10 or later and a Java runtime for Hail:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pip install pre-commit && pre-commit install
```

gnomad_methods is installed from GitHub at the commit pinned in
`requirements.txt`, because the pipeline depends on that version of its API.

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
hailctl dataproc submit <cluster> gnomad_constraint/pipeline/constraint_pipeline.py \
  --pyfiles /tmp/pyfiles.zip -- --compute-constraint-metrics --test
```

Run a step with `--test` and check its output before running it at full scale.

## Continuous integration

CI (`.github/workflows/ci.yml`) runs on pull requests and pushes to `main`. It
checks formatting (black, isort, autopep8), docstrings (pydocstyle), and pylint
errors on the Python code, and styler and lintr on the R code. The same
formatters run as pre-commit hooks.
