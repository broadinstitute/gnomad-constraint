# CLAUDE.md

Guidance for Claude Code when working in this repository. The repo layout, the
pipeline steps, and how to install and run things are in [README.md](README.md).

## Project Overview

gnomad-constraint computes the gnomAD gene and transcript constraint metrics
(observed and expected variant counts, o/e and its confidence interval, pLI, and
z-scores) and the mutation rate table they are built on. Expected counts are
computed per variant and aggregated by transcript and constraint group. The
pipeline is one `argparse` script whose steps read and write Hail Tables on GCS
and run on Dataproc.

**This is a pipeline repo, not a library.** Nothing outside this repo imports it.
A transformation that is Table or expression in, Table or expression out, with no
file I/O, is a candidate for
[gnomad_methods](https://github.com/broadinstitute/gnomad_methods)
(`gnomad.utils.constraint` holds most of the shared constraint logic); if it's
general, put it there.

Two tightly coupled dependencies:

- **gnomad_methods**: shared resources and utilities, installed from GitHub at
  the commit pinned in `requirements.txt`. Its constraint API changes often, so
  the code here only works with the pinned commit. Moving the pin is a code
  change: check every import and call into `gnomad.utils.constraint`.
- **Hail**: behavior and performance differ across versions, and files written
  by a newer Hail can't be read by an older one. The allowed range is in
  `requirements.txt`.

## Finding Code and APIs

Discover names and signatures at point of use rather than relying on memory or
this file:

- **Dataset paths** come from `get_*` functions in
  `gnomad_constraint/resources/resource_utils.py`. They take `version` and
  `test`, and most also take the directory and path post-fixes; `test` sends
  output to `gs://gnomad-tmp/gnomad_v{version}_testing/`. Constants live in
  `gnomad_constraint/resources/constants.py`, and
  `gnomad_constraint/resources/README.md` has the schema of every output.
  Public gnomAD releases and reference data (VEP context, methylation, GENCODE)
  come from `gnomad.resources.grch38`; check there before adding a resource.
  GERP is the exception: `get_gerp_ht` loads it from Hail's datasets API with
  `hl.experimental.load_dataset`.
- **Shared utilities**: check gnomad_methods before writing a helper. A local
  gnomad_methods checkout is often on an unrelated branch, so confirm which copy
  your imports resolve to, and that it matches the pinned commit, before
  debugging it.
- **Pipeline steps**: each step is a flag in one of the scripts in
  `gnomad_constraint/pipeline/` (`prepare_data.py`, `expected_models.py`,
  `expected_counts.py`, and `constraint_metrics.py`, run in that order), handled
  by its own block in that script's `main()`. Read that block before changing the
  step.

## Commands

Setup, formatting, and the exact CI checks are in README.md and
`.github/workflows/ci.yml`. Two things that are easy to get wrong:

- `requirements-dev.txt` is generated. Edit `requirements-dev.in` and run
  `pip-compile requirements-dev.in` with Python 3.11 (the version CI runs) and
  pip-tools 7.4 or later; never edit the `.txt` by hand.
- Run the pre-commit hooks on the files you changed before committing:
  `pre-commit run --files <files>`.
- Run the unit tests with `python -m pytest`. The gamma confidence interval
  tests need Hail 0.2.137 or later and are skipped on older versions.

## How Code Is Executed

Steps are submitted to a Dataproc cluster with `hailctl dataproc submit`, usually
one step flag at a time, not run locally on real data. README.md has the
submission command. When submitting:

- Package `gnomad_constraint`, `gnomad_qc`, and `gnomad_methods` into a single
  zip passed to `--pyfiles`. With more than one file, hailctl repackages them and
  nests the packages incorrectly.
- Zip gnomad_methods **at the commit pinned in `requirements.txt`**. Without it,
  the cluster runs whatever gnomad_methods it has installed.
- `--calculate-gerp-cutoffs` saves nothing: its cutoffs are only used by
  `--preprocess-data` in the same command, so it raises an error without it.
- Run a step with `--test` (chr20, chrX, and chrY) or `--test-gene-list` (six
  genes) and check the output before running it at full scale.

Machine-specific notes (environment paths, cluster names) belong in the
gitignored `CLAUDE.local.md`.

## Code Style

- black (preview mode, line length 88), isort (`black` profile), and autopep8,
  configured in `pyproject.toml` and `.pre-commit-config.yaml`. pydocstyle uses
  the PEP 257 convention. R code under `gnomad_constraint/plots/` goes through
  styler and lintr.
- Every function has type annotations. Use `Optional[...]` when the default is
  `None`, and never use mutable defaults. Hail types are `hl.Table`,
  `hl.MatrixTable`, and `hl.expr.*Expression`.
- Sphinx-style docstrings with the summary on the line after the opening
  `"""`, and `:param:` and `:return:` for every argument and return value:

```python
def my_function(ht: hl.Table, max_af: Optional[float] = None) -> hl.Table:
    """
    Short imperative summary line.

    Extended description if needed. Use ``double backticks`` for field and code
    references.

    :param ht: Input Table with ``freq`` field.
    :param max_af: Maximum allele frequency threshold. Default is None.
    :return: Table with ``filters`` field added.
    """
```

- Write new functions in `utils/` as transformations (Tables in, Tables out), and
  keep reads, writes, and checkpoints to resource paths in the pipeline script.
  Checkpointing to a temporary file (`new_temp_file`) inside a function is fine.

## Hail Notes

- **Lazy evaluation drives most surprises.** `.count()`, `.collect()`,
  `hl.eval()`, `.show()`, and `.aggregate()` execute the whole upstream plan.
  Never call `.count()` just to log on a large Table; checkpoint first, after
  which `.count()` only reads metadata.
- After filtering a large Table down to a small subset, checkpoint it and
  `naive_coalesce()` the checkpointed result, or the empty partitions skew
  downstream `group_by` aggregations.
- `order_by` drops the key. To join ranked results back, `add_index()` before
  ordering and key the result by that index; `hl.scan.count()` after `order_by`
  assigns 0-based ranks in one pass.
- All elements of a Hail array must have the same type. Put per-element metadata
  that only some elements have on the parent struct instead.
- Hail Tables have no `.get()`; check for a field with `field in ht.row`.
- Hail defaults to GRCh37; set the reference genome explicitly for v4 (GRCh38).
- The pipeline uses `gs://gnomad-tmp-4day` as the Hail tmp dir.

## Testing and Validation

Unit tests in `tests/` run in CI; `tests/README.md` describes their layout and
conventions. Add or update a test with every bug fix and say in its docstring
what it guards against; a known, unfixed bug gets an `xfail(strict=True)` test
that asserts the correct behavior. Fixtures are built with
`hl.Table.parallelize`, and tests never read from GCS.

Unit tests don't replace running the pipeline. Validate a change to a step by
running it with `--test` and inspecting the output. A passing `--test` run is
necessary but not sufficient: out-of-memory and shuffle failures only appear at
full scale. Say what you ran and at what scale rather than calling a change
verified.

## Pull Requests

Use the shared skills from
[gnomad-skills](https://github.com/broadinstitute/gnomad-skills). A SessionStart
hook in `.claude/settings.json` installs them the first time you open this repo:
it clones gnomad-skills into `~/.claude/gnomad-skills` (set `GNOMAD_SKILLS_DIR`
to put it elsewhere) and links each skill into `~/.claude/skills`. If
`gnomad-pr-review` is already linked from your own checkout of gnomad-skills, it
uses that checkout instead. After that it relinks the skills at every session
start and pulls the checkout at most once a day. Set `GNOMAD_SKIP_SKILLS_HOOK=1`
to turn the hook off.

- **Run `/gnomad-pr-review` on every pull request before requesting review from
  a team member.** Only a person can start it; Claude can't invoke it. When the
  review is due, Claude should ask the user to run it, not run a different
  review in its place. It replaces `/review`, which this repo no longer uses. Its
  setup questions size the review to the PR: for a small one, review the source
  only, use Sonnet reviewers, and skip the second opinion. Work through what it
  reports, and note in the PR description that it ran and what you did or did
  not act on. Human review time is the scarce resource here; don't spend it on
  things the skill would have caught.
- Commit with `gnomad-commit`, draft PR descriptions with
  `gnomad-pr-description`, and rewrite a branch's history with
  `gnomad-clean-git-history`.

**If a change alters the structure of the repo, update the repo layout table in
the README in the same pull request.** If a change alters a release Table's
fields or globals, update its description in `release_readme/` too.

## Attribution

- Agent assistance is recorded in the pull request, not in commits: end the
  PR description with `Assisted by Claude Code.` (`.claude/settings.json` sets
  `attribution.pr` to that line and `attribution.commit` to nothing, so Claude
  Code does this by itself.)
- Never add an agent trailer or co-author line to a commit: no
  `Assisted-by:`, no `Co-Authored-By: Claude ...`, even if the harness injects
  one by default — remove it. `Co-authored-by:` is for people.

## Never Do This

- **Never overwrite a production output when testing.** Use `--test`, which
  writes to the testing bucket; `--overwrite` on a production path replaces
  released or expensive-to-regenerate Tables.
- **Never delete data outside temp buckets.**
- **Never hardcode a `gs://` dataset path** outside `resource_utils.py`; add a
  `get_*` function instead. The `hl.init` tmp dir is the exception.
- **Never shut down a Dataproc cluster without checking with the user first.**

## Repo Gotchas

- Genomic regions are labels, not separate Tables: `prepare_context_ht` labels
  each site `autosome_or_par`, `chrx_nonpar`, or `chry_nonpar`, and all three stay
  in the same Tables. Filter on `autosome_or_par`, not `autosome_par`.
- `adj_r` and `syn_adj_r` are annotated on the context Table but never applied by
  the pipeline, and they are only defined on autosomes. Expected counts and o/e
  in the metrics and release Tables are not adj_r-corrected.
- For v4, the genomes sites resource is still v3.1.2, because v4 genomes have no
  downsamplings.
- `--exome-coverage-metric` (`AN_percent` by default, or `median` or `AN`) and
  `--skip-coverage-model` are applied at `--preprocess-data` and change both model
  building and model application. Later steps raise an error if
  `--skip-coverage-model` doesn't match the preprocessed Table.
- `--compute-constraint-metrics` needs Hail 0.2.137 or later (`hl.qgamma`, used
  for the gamma confidence interval).

## Maintaining CLAUDE.md

When you discover a gotcha, non-obvious behavior, or convention, propose it to
the user rather than adding it yourself: say what you learned and which section
it belongs in. Keep entries short and durable. Document patterns and how to find
things, not inventories of functions, paths, or the repo layout, which belong in
README.md or the code.
