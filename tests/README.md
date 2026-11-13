# tests

Unit tests for `gnomad_constraint`. `pyproject.toml` sets `testpaths = ["tests"]`,
so only this directory is collected.

## Running

```sh
python -m pytest
python -m pytest -q tests/utils/test_constraint_metrics.py
```

The tests need a JVM for Hail and the packages in `requirements.txt`, including
gnomad_methods at the pinned commit; a different gnomad_methods can fail to import
or give different results. `conftest.py` starts Hail once for the whole
session, before tests are collected, on the local Spark backend even when
`~/.config/hail/config.ini` selects Hail Batch. It gives the Spark driver a 4 GB
heap (the 1 GB default runs out of memory), binds it to `127.0.0.1` so it starts
on machines whose hostname does not resolve, and prints the path of the session's
Hail log in the test header. Spark fails to start on Java 23 or later; CI uses
Java 11. CI runs the suite on every pull request; it takes
about 1.5 minutes locally.

The gamma confidence interval needs `hl.qgamma` (Hail >= 0.2.137). On older Hail
the tests that compute it are skipped.

## Layout

Test directories mirror the package. `gnomad_constraint/utils/constraint.py` is
large, so its tests are split by pipeline stage into
`tests/utils/test_constraint_<stage>.py`. Each tested function gets a
`Test<FunctionName>` class, with two exceptions: path helpers that take the
same arguments share one parametrized class, and a test that checks two
implementations against each other is named for the comparison.

| Module | Covers |
|---|---|
| `utils/test_constraint_preprocess.py` | `prepare_context_ht`, frequency filtering, mutation rate site selection, exome coverage metrics, exome observed/possible counts, build-model groups, `prepare_ht_for_constraint_calculations`, the training set, and the mutation rate |
| `utils/test_constraint_apply.py` | Applying the models with and without a coverage model, the per-variant and aggregated expected paths (checked against each other), and `aggregate_by_constraint_groups` |
| `utils/test_constraint_metrics.py` | o/e, both CIs, z-scores and flags, transcript selection, ranks and percentile bins, pLI per LoF group, `compute_constraint_metrics`, and the release Tables (values, globals, and the published TSV columns) |
| `utils/test_constraint_gene_quality.py` | Coverage, mapping quality, segdup/LCR, and flag thresholds in `compute_gene_quality_metrics` |
| `utils/test_import_data.py` | Formatting the imported de novo mutations and the genomic feature windows |
| `resources/test_resource_utils.py` | Parameter checks, resource version mapping, every pipeline path (production, test, and post-fixes), the step wiring of each script's resource collection, and `filter_for_test` |
| `pipeline/test_constraint_pipeline.py` | The checks in the pipeline scripts: `--skip-coverage-model` consistency, `--calculate-gerp-cutoffs` without `--preprocess-data`, the coverage model scale, the constraint group input warning, and the release downsamplings |

## How these tests are written

- **Fixtures are hand-written** with `hl.Table.parallelize`, not slices of real
  data, so tests run in seconds with no GCS access.
- **Cases include missing and boundary inputs**: sites absent from the exomes,
  filtered sites, sites with no coverage, groups with no expected variants.
- **Where two implementations compute the same thing, they are tested against
  each other**: the per-variant and aggregated expected paths must give the same
  sums, and the gamma CI is checked against `scipy.stats.gamma`.
- **Regression tests say what they guard against** in their docstring, with the
  issue number when there is one.
- **Known bugs are `xfail(strict=True)`, not pinned.** The test asserts the correct
  behavior and its `reason` describes the bug, so fixing it turns the test into an
  unexpected pass that has to be acknowledged. Behavior that is pinned as-is,
  because it looks intended but was not checked against a spec, says so in the
  docstring.
- **Expensive Hail runs are shared.** Fixtures that run a whole pipeline step are
  module- or class-scoped, or cached, so each runs once.
- **`data/` holds schemas, not data.** `constraint_metrics_row_type.txt` is the
  row type of the internal constraint metrics Table, and
  `release_tsv_columns.txt` is the header of the published v4.1.2 release TSV.
  The release test builds an empty Table with that schema and checks that
  `prepare_release_ht` and `flatten_release_ht` produce exactly those columns,
  and `TestComputeConstraintMetrics` checks the row type against what
  `compute_constraint_metrics` and `compute_constraint_percentile_bins` write.
  Update both files when the schema changes on purpose.

## Not covered yet

Anything that reads real inputs or outputs on GCS: `get_gencode_ht`,
`get_gencode_cds_ht`, `get_gerp_ht`, and the pipeline `main()`, whose argument
parser is only built under `if __name__ == "__main__"`.
