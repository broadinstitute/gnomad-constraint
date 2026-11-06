#!/bin/bash
# Plateau-model test, without rr: build per-variant multipliers that swap the stored
# plateau lines for the saturating curve, the logistic curve and one line per
# substitution type (build_plateau_multipliers.py), then run the constraint pipeline's
# aggregation and constraint groups with each (--skip-regional-rr
# --variant-rate-ht-path) under rr_downstream/plateau_<version>. Metrics and the release
# need Hail 0.2.137+, so they run with run_plateau_finish.sh on a newer cluster. Each step is
# skipped if its output exists; nothing is overwritten.
#   CLUSTER=jg1 run_plateau.sh build saturating
#   CLUSTER=jg2 run_plateau.sh logistic
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg1}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_plateau_$$.zip
M=gs://gnomad-tmp/jg/rr_plateau
OUT=gs://gnomad-tmp/jg/rr_downstream
C=gnomad_v4.1.1/constraint

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag on $CLUSTER"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}
done_ht() { gsutil -q stat "$1/_SUCCESS"; }

while [ "$(gcloud dataproc clusters describe "$CLUSTER" --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 60; done
echo "=== $CLUSTER is running"
source "$F/cluster/snapshot_code.sh"

for step in "$@"; do
  if [ "$step" = build ]; then
    gsutil -q stat "$M/plateau_params.json" || gsutil -q cp "$F/plateau_params.json" "$M/plateau_params.json" || exit 1
    done_ht "$M/gnomad.v4.1.1.plateau_by_substitution_multiplier.ht" ||
      submit plateau_build "$S/build_plateau_multipliers.py" --params "$M/plateau_params.json" --out "$M" || exit 1
    continue
  fi
  rates=$M/gnomad.v4.1.1.plateau_${step}_multiplier.ht
  until done_ht "$rates"; do sleep 300; done
  root=$OUT/plateau_$step
  done_ht "$root/$C/apply_models/transcript_consequences/gnomad.v4.1.1.constraint_group.ht" ||
    submit "constraint_plateau_$step" "$P/constraint_pipeline.py" --aggregate-per-variant-expected \
      --skip-regional-rr --variant-rate-ht-path "$rates" --aggregate-by-constraint-groups \
      --output-root "$root" || exit 1
done
echo "=== $(date '+%H:%M:%S') all done on $CLUSTER"
