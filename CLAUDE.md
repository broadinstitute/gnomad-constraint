# gnomad-constraint Claude Reference

## Code Style

### Formatting

Code is formatted with **black** (preview mode, line length 88), **isort** (profile `"black"`), and **autopep8** (aggressive=1, ignoring E201/E202/E203/E731). Linting uses **pylint** and **pydocstyle** (PEP 257 convention, ignoring D100/D104). Config is in `pyproject.toml`.

```bash
# Manual formatting
black gnomad_constraint/
isort --profile black gnomad_constraint/
```

### Docstrings

Use **Sphinx-style** (`:param:`, `:return:`) docstrings following the gnomad_methods convention (see gnomad_mnv CLAUDE.md for full examples).

### Type Annotations

- **All functions** must have type annotations on parameters and return values.
- Use `typing.List`, `typing.Optional`, etc. for generic types.
- For Hail expression parameters, use `hl.expr.StructExpression`, `hl.expr.BooleanExpression`, etc.
- For Hail table/matrix types, use `hl.Table`, `hl.MatrixTable`.

## Key Constants (`constants.py`)

```python
VERSIONS = ["2.1.1", "4.0", "4.1", "4.1.1"]
CURRENT_VERSION = "4.1.1"
DATA_TYPES = ["context", "exomes", "genomes"]
MODEL_TYPES = ["plateau", "coverage"]
GENOMIC_REGIONS = ["autosome_par", "chrx_nonpar", "chry_nonpar"]
POPS = ("global", "afr", "amr", "eas", "nfe", "sas")
COVERAGE_CUTOFF = 40
CUSTOM_VEP_ANNOTATIONS = ["transcript_consequences", "worst_csq_by_gene"]
MU_GROUPING = ("context", "ref", "alt", "methylation_level")
```

## Known Gotchas

- **MU_GROUPING not exported**: `gnomad_constraint.resources.resource_utils` does NOT export `MU_GROUPING`. It must be defined locally as `("context", "ref", "alt", "methylation_level")`.
- **v4 drops chrY/chrX**: The v4 pipeline removes `chry_nonpar` and `chrx_nonpar` from regions early in `main()`. Don't assume all 3 genomic regions are present.
- **Coverage metric**: v4 can use `"exomes_AN_percent"` instead of `"exome_coverage"`. This affects model building and application.
- **Genomes v3.1.2 for v4**: Even in v4, the genomes sites resource uses v3.1.2 (downsamplings dropped in v4).
- **Hail `hl.init` tmp_dir**: Pipeline uses `gs://gnomad-tmp-4day` as the Hail temp directory. Ensure this bucket exists and is writable.
- **`ht.get()` doesn't exist on Hail Tables**: Use `field in ht.row` to check field existence.
- **constraint_metrics.ht pLoF lacks adj_r**: The pre-computed `constraint_metrics.ht` has `lof.exp` WITHOUT the regional depletion correction (`adj_r`). To get adj_r-corrected pLoF, compute it from the per-SNV table by filtering to LOFTEE HC + `possible_variants == 1` and aggregating with `expected_variants[0] * adj_r`.
- **Table version awareness**: Undated tables (e.g., `annotate_with_oe.ht`) and dated tables (e.g., `annotate_with_oe.12_23_25.ht`) may have different schemas and column names. Always verify which version you're using.

## Hail / Dataproc Best Practices

- **Never use `.count()` for logging on large tables**: `count()` forces a full table materialization. On a large per-SNV table this triggers a massive Spark job and can cause shuffle failures. Use it only when the result is actually needed for computation.
- **Use `naive_coalesce()` after aggressive filters**: When filtering a large table down to a small subset (e.g., LOFTEE HC LoF from all variants), most partitions become empty. This causes shuffle skew in downstream `group_by` aggregations. Call `naive_coalesce(200)` after the filter to rebalance.
- **Per-SNV table `calibrate_mu` struct**: The per-variant expected table (`gnomad.v4.1.per_variant_expected.coverage_corrected.with_downsamplings.ht`) stores transcript-level fields (`gene`, `transcript`, `canonical`, `modifier`, `observed_variants`, `expected_variants`, `possible_variants`) inside a `calibrate_mu` struct. Flatten it with `ht = ht.annotate(**ht.calibrate_mu)` before accessing those fields.
- **Log full row field lists for debugging**: When debugging schema issues, log `list(ht.row)` (all fields), not `list(ht.row)[:20]` (truncated). Important fields like `calibrate_mu` may be beyond the first 20.
- **`order_by` destroys the key — use `add_index` to rekey cheaply**: After `ht.order_by(expr)`, the table is unkeyed. To rejoin ranked results back to the original table, call `ht.add_index("_rank_idx")` before ordering, then `rank_ht.key_by("_rank_idx")` after. Rejoining via an integer index is an O(1) lookup vs a full key scan.
- **`hl.scan.count()` for rank assignment**: After `order_by`, annotate with `hl.scan.count()` to assign 0-based ascending ranks in a single pass: `ht = ht.order_by(ht.val).annotate(rank=hl.scan.count())`.
- **Checkpoint small select-then-order tables, not the full wide table**: When computing ranks for many `(group, field)` combinations, select only the columns needed (`ht.select("_rank_idx", _val=expr)`), order, rank, checkpoint, and join results back in one pass. Avoids checkpointing the full wide table once per iteration.
- **`.count()` after checkpoint is free**: `count()` on a checkpointed table reads already-materialized metadata rather than re-executing the query. Place `count()` after a checkpoint to avoid computing the table twice.
- **Python list comprehensions over Hail arrays for indexed access**: When you need to index a Hail array with a known Python integer (e.g., `ht.constraint_groups[i]`), use a Python list comprehension rather than `hl.enumerate` + lambda. This also allows Python-time dict lookups like `rank_hts[(i, key)]` inside the expression.
- **`hl.Table.parallelize` to reconstruct a small HT from collected data**: `hl.Table.parallelize(hl.eval(ht.my_array_global), schema=ht.my_array_global.dtype.element_type).key_by(...)` reconstructs a small Hail Table from a global array without re-running any jobs.
- **Hail array elements must share the same struct schema**: All elements of a Hail array field must have identical types. You cannot annotate only `array[0]` with extra fields while leaving `array[1+]` unchanged — Hail will reject the mixed schema. Instead, promote such metadata to the parent struct level (e.g., add a `{field}_rank` struct directly on the constraint group rather than inside `oe_info[0]`).

## Dataproc Submission

**Important**: hailctl repackages `--pyfiles` into a temp zip using `os.walk`, which nests packages incorrectly. Use the single-zip workaround (same as gnomad_mnv):

```bash
# Build a single zip with correct top-level package structure
cd <gnomad-constraint root> && \
  rm -f /tmp/pyfiles.zip && \
  zip -r /tmp/pyfiles.zip gnomad_constraint/ -x '*.pyc' '*__pycache__*' && \
  cd <gnomad_qc root> && \
  zip -r /tmp/pyfiles.zip gnomad_qc/ -x '*.pyc' '*__pycache__*' '*.DS_Store'

# Submit to cluster (single zip = used directly, not repackaged)
hailctl dataproc submit <CLUSTER> \
  gnomad_constraint/pipeline/constraint_pipeline.py \
  --pyfiles /tmp/pyfiles.zip \
  -- --compute-constraint-metrics --test --overwrite
```

## gnomad_methods / gnomad_qc API

See the gnomad_mnv CLAUDE.md for shared API reference (`public_release`, `TableResource`, `get_gnomad_v4_vds`, etc.). Key constraint-specific imports:

```python
from gnomad.utils.constraint import build_models, compute_pli, oe_confidence_interval
from gnomad.resources.grch38.gnomad import public_release, DOWNSAMPLINGS, all_sites_an
from gnomad_qc.resource_utils import PipelineResourceCollection, PipelineStepResourceCollection
```
