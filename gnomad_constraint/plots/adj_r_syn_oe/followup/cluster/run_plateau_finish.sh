#!/bin/bash
# The steps after run_plateau.sh for the plateau versions. Only the constraint metrics
# need Hail 0.2.137+ (hl.qgamma); everything else runs on a 0.2.134 cluster. Each step
# waits for its inputs and is skipped if its output exists; nothing is overwritten.
#   CLUSTER=jg1     run_plateau_finish.sh groupoe          (all versions' group obs/exp)
#   CLUSTER=jg4-137 run_plateau_finish.sh metrics          (needs 0.2.137+)
#   CLUSTER=jg2     run_plateau_finish.sh release slim     (release, TSV, slim copies)
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg1}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_plateau_finish_$$.zip
OUT=gs://gnomad-tmp/jg/rr_downstream
C=gnomad_v4.1.1/constraint
VERSIONS=${VERSIONS:-"plateau_saturating plateau_logistic plateau_by_substitution"}

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag on $CLUSTER"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}
done_ht() { gsutil -q stat "$1/_SUCCESS"; }
wait_for() { until "$@"; do sleep 300; done; }
groups() { echo "$OUT/$1/$C/apply_models/transcript_consequences/gnomad.v4.1.1.constraint_group.ht"; }
metrics() { echo "$OUT/$1/$C/metrics/transcript_consequences/gnomad.v4.1.1.constraint_metrics.ht"; }
release() { echo "$OUT/$1/$C/release/gnomad.v4.1.1.constraint_metrics.ht"; }
tsv() { echo "$OUT/$1/$C/release/gnomad.v4.1.1.constraint_metrics.tsv.bgz"; }

while [ "$(gcloud dataproc clusters describe "$CLUSTER" --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 60; done
echo "=== $CLUSTER is running"
source "$F/cluster/snapshot_code.sh"

for step in "$@"; do
  case $step in
    groupoe)
      for v in $VERSIONS; do wait_for done_ht "$(groups "$v")"; done
      gsutil -q stat "$OUT/group_oe/plateau_by_substitution.group_oe.tsv.bgz" ||
        submit export_group_oe_plateau "$S/export_group_oe.py" --versions $VERSIONS || exit 1
      ;;
    metrics)
      for v in $VERSIONS; do
        wait_for done_ht "$(groups "$v")"
        done_ht "$(metrics "$v")" ||
          submit "metrics_$v" "$P/constraint_pipeline.py" --compute-constraint-metrics --output-root "$OUT/$v" || exit 1
      done
      ;;
    release)
      for v in $VERSIONS; do
        wait_for done_ht "$(metrics "$v")"
        steps=()
        done_ht "$(release "$v")" || steps+=(--prepare-release)
        gsutil -q stat "$(tsv "$v")" || steps+=(--export-release-tsv)
        [ ${#steps[@]} -eq 0 ] ||
          submit "release_$v" "$P/constraint_pipeline.py" "${steps[@]}" --output-root "$OUT/$v" || exit 1
      done
      ;;
    slim)
      for v in $VERSIONS; do wait_for gsutil -q stat "$(tsv "$v")"; done
      submit export_release_slim_plateau "$S/export_release_slim.py" --versions $VERSIONS || exit 1
      ;;
  esac
done
echo "=== $(date '+%H:%M:%S') all done on $CLUSTER"
