#!/bin/bash
# Finish the constraint metrics of the downstream versions from the checkpoints their
# failed --compute-constraint-metrics runs left in gs://gnomad-tmp-4day (matched to
# versions by resume_constraint_metrics.py --map), then write the release Table and
# TSV. One version at a time on jg2; nothing is rerun or overwritten.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg2}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_resume_run.zip
OUT=gs://gnomad-tmp/jg/rr_downstream
C=gnomad_v4.1.1/constraint
T=gs://gnomad-tmp-4day/constraint_metrics

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

source "$F/cluster/snapshot_code.sh"

# version:checkpoint:stage (a list, not a heredoc, so the submit commands can't read it).
RUNS=(
  ridge_qcexcl:z_scores-pSIb1MqjzsTGOdKrOcZvBh.ht:z_scores
  ridge_full:oe_ci_z-pGTSfqVtiNQFKuD0qrqMjD.ht:oe_ci_z
  none:oe_ci_z-La7lV7j5SMT4rsOrLS0EzH.ht:oe_ci_z
  bonferroni_full:oe_ci_z-bDX3AxW8he0o6EFLNSeRBv.ht:oe_ci_z
  bonferroni_qcexcl:z_scores-83JJfeDptBjNvCZuYEM5kC.ht:z_scores
  roulette_ridge:z_scores-mK2ygh3xjDNSnHQwUGtocK.ht:z_scores
)
# Or only the runs given as arguments, in the same form.
[ $# -gt 0 ] && RUNS=("$@")
for run in "${RUNS[@]}"; do
  IFS=: read -r v ck stage <<< "$run"
  metrics=$OUT/$v/$C/metrics/transcript_consequences/gnomad.v4.1.1.constraint_metrics.ht
  gsutil -q stat "$metrics/_SUCCESS" ||
    submit "resume_metrics_$v" "$S/resume_constraint_metrics.py" --resume --checkpoint "$T.$ck" --stage "$stage" --output-root "$OUT/$v" < /dev/null || continue
  steps=()
  gsutil -q stat "$OUT/$v/$C/release/gnomad.v4.1.1.constraint_metrics.ht/_SUCCESS" || steps+=(--prepare-release)
  gsutil -q stat "$OUT/$v/$C/release/gnomad.v4.1.1.constraint_metrics.tsv.bgz" || steps+=(--export-release-tsv)
  if [ ${#steps[@]} -gt 0 ]; then
    submit "release_$v" "$P/constraint_pipeline.py" "${steps[@]}" --output-root "$OUT/$v" < /dev/null
  fi
done
echo "=== $(date '+%H:%M:%S') all done"
