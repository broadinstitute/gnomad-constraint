#!/bin/bash
# Constraint metrics, release Table and release TSV for the downstream runs, on a
# cluster with Hail 0.2.137+ (the metrics step needs hl.qgamma). For each version in
# rr_downstream given as an argument, waits for its constraint groups and runs only the
# steps whose outputs aren't written yet; nothing is rerun or overwritten. "none" also
# runs the aggregation and constraint groups without rr (--skip-regional-rr), as the
# baseline for the LOEUF checks.
#   run_metrics.sh ridge_qcexcl bonferroni_full ...
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg3-137}
MAX_WAIT_HOURS=${MAX_WAIT_HOURS:-24}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_metrics_$$.zip
OUT=gs://gnomad-tmp/jg/rr_downstream
C=gnomad_v4.1.1/constraint

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}
done_ht() { gsutil -q stat "$1/_SUCCESS"; }

while [ "$(gcloud dataproc clusters describe "$CLUSTER" --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 60; done
echo "=== $CLUSTER is running"
source "$F/cluster/snapshot_code.sh"

run_version() {
  local v=$1 root=$OUT/$1
  local groups=$root/$C/apply_models/transcript_consequences/gnomad.v4.1.1.constraint_group.ht
  local metrics=$root/$C/metrics/transcript_consequences/gnomad.v4.1.1.constraint_metrics.ht
  local release=$root/$C/release/gnomad.v4.1.1.constraint_metrics.ht
  local tsv=$root/$C/release/gnomad.v4.1.1.constraint_metrics.tsv.bgz
  if [ "$v" = none ] && ! done_ht "$groups"; then
    local agg=(--aggregate-by-constraint-groups)
    done_ht "$root/$C/apply_models/transcript_consequences/gnomad.v4.1.1.per_variant_expected.aggregated.ht" ||
      agg=(--aggregate-per-variant-expected --skip-regional-rr "${agg[@]}")
    submit "groups_$v" "$P/constraint_pipeline.py" "${agg[@]}" --output-root "$root" || return 1
  fi
  local waited=0
  until done_ht "$groups"; do
    [ $waited -ge $((MAX_WAIT_HOURS * 60)) ] && { echo "=== $v: no constraint groups after ${MAX_WAIT_HOURS} h"; return 1; }
    sleep 300; waited=$((waited + 5))
  done
  local steps=()
  done_ht "$metrics" || steps+=(--compute-constraint-metrics)
  done_ht "$release" || steps+=(--prepare-release)
  gsutil -q stat "$tsv" || steps+=(--export-release-tsv)
  if [ ${#steps[@]} -eq 0 ]; then echo "=== $v already done"; return 0; fi
  submit "metrics_$v" "$P/constraint_pipeline.py" "${steps[@]}" --output-root "$root"
}

pids=()
for v in "$@"; do run_version "$v" & pids+=($!); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "=== $(date '+%H:%M:%S') finished (rc=$rc)"
for v in "$@"; do gsutil -q stat "$OUT/$v/$C/release/gnomad.v4.1.1.constraint_metrics.tsv.bgz" && echo "    $v: release TSV written"; done
exit $rc
