#!/bin/bash
# Export the variant cells for the single-model try (export_model_cells.py): the
# synonymous cells to fit on, then the missense and pLoF cells to apply the fits to;
# a chr22 test first each time, then all autosomes, copied into data/one_model/.
# Nothing is overwritten.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg1}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_one_model_$$.zip
O=gs://gnomad-tmp/jg/one_model
mkdir -p "$F/data/one_model"

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag on $CLUSTER"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

while [ "$(gcloud dataproc clusters describe "$CLUSTER" --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 60; done
echo "=== $CLUSTER is running"
source "$F/cluster/snapshot_code.sh"
gsutil -q stat "$O/syn_model_cells.chr22.window_mq.tsv.bgz" ||
  submit model_cells_test "$S/export_model_cells.py" --test --out "$O/syn_model_cells.chr22.tsv.bgz" || exit 1
for f in syn_model_cells.chr22.tsv.bgz syn_model_cells.chr22.window_mq.tsv.bgz; do [ -e "$F/data/one_model/$f" ] || gsutil -q cp "$O/$f" "$F/data/one_model/"; done
gsutil -q stat "$O/syn_model_cells.window_mq.tsv.bgz" ||
  submit model_cells "$S/export_model_cells.py" --out "$O/syn_model_cells.tsv.bgz" || exit 1
for f in syn_model_cells.tsv.bgz syn_model_cells.window_mq.tsv.bgz; do [ -e "$F/data/one_model/$f" ] || gsutil -q cp "$O/$f" "$F/data/one_model/"; done
# The missense and pLoF cells, to apply the fitted models to (one_model_fit.py
# --apply-cells). The window MQ summary is the same as for the synonymous cells.
gsutil -q stat "$O/mislof_model_cells.chr22.tsv.bgz" ||
  submit model_cells_mislof_test "$S/export_model_cells.py" --test --csqs mis lof --skip-window-mq \
    --out "$O/mislof_model_cells.chr22.tsv.bgz" || exit 1
[ -e "$F/data/one_model/mislof_model_cells.chr22.tsv.bgz" ] || gsutil -q cp "$O/mislof_model_cells.chr22.tsv.bgz" "$F/data/one_model/"
gsutil -q stat "$O/mislof_model_cells.tsv.bgz" ||
  submit model_cells_mislof "$S/export_model_cells.py" --csqs mis lof --skip-window-mq \
    --out "$O/mislof_model_cells.tsv.bgz" || exit 1
[ -e "$F/data/one_model/mislof_model_cells.tsv.bgz" ] || gsutil -q cp "$O/mislof_model_cells.tsv.bgz" "$F/data/one_model/"
echo "=== $(date '+%H:%M:%S') all done"
