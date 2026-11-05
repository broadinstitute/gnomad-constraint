#!/bin/bash
# Put the per-variant-rate versions (Roulette, home-built) on the cell-level scoring:
# add their cell factors to Ruchit's exported cells (export_rate_cells.py) and score
# them like the full run on the driver (run_full_run_scores.py --rr-names), then copy
# rr_rates.* into followup/. Waits for the cluster; nothing is overwritten.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg1}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_rate_cells.zip
CELLS=gs://gnomad-tmp/jg/rr_downstream/syn_cells_rr_rates.tsv.bgz
OUT=gs://gnomad-tmp/jg/rr_downstream/full_run_scores/out
NAMES="roulette_only roulette_ridge homebuilt_only homebuilt_ridge homebuilt_regional homebuilt_sequence"

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

while [ "$(gcloud dataproc clusters describe "$CLUSTER" --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 60; done
echo "=== $CLUSTER is running"
source "$F/cluster/snapshot_code.sh"
M=$(mktemp -d /tmp/rr_ratemods.XXXX)
cp "$F/full_run_scores.py" "$F/remap_rr_to_shet.py" "$F/rr_shrinkage_experiments.py" "$M/"
(cd "$M" && zip -q mods.zip full_run_scores.py remap_rr_to_shet.py rr_shrinkage_experiments.py)

gsutil -q stat "$CELLS" ||
  submit rate_cells "$S/export_rate_cells.py" --pyfiles "$ZIP" -- --out "$CELLS" || exit 1
# The driver scoring needs scikit-learn and statsmodels on a fresh cluster.
submit rate_install "$S/install_driver_packages.py" || exit 1
gsutil -q stat "$OUT/rr_rates.metrics.tsv" ||
  submit rate_scores "$S/run_full_run_scores.py" --pyfiles "$M/mods.zip" -- --cells "$CELLS" --out "$OUT" --prefix rr_rates --rr-names $NAMES || exit 1
for f in $(gsutil ls "$OUT/rr_rates.*"); do
  [ -e "$F/$(basename "$f")" ] || gsutil -q cp "$f" "$F/"
done
tr '\r' '\n' < "$D/logs/rate_cells.log" | grep "Base cells without" || true
echo "=== $(date '+%H:%M:%S') all done"
