#!/bin/bash
# Score the background runs (run_backgrounds.sh) on jg2: once every arm's rr is written,
# add each arm's ridge rr (default penalty and 10 per DNM; all chromosomes and held out)
# to Ruchit's exported synonymous cells, score them like the full run on the driver
# (run_full_run_scores.py --rr-names) and copy rr_bg_pipeline.* into followup/.
# Nothing is overwritten.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg2}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
BG=gs://gnomad-tmp/jg/rr_bg
RD=gnomad_v4.1.1/constraint/regional_rate
CELLS=gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export/syn_cells.tsv.bgz
OUT_CELLS=$BG/syn_cells_rr_backgrounds.tsv.bgz
OUT=$BG/scores
ARMS="absent_locus all_locus absent_variant all_variant ctx_proportional ctx_even"
CUTS="ridge ridge_dnm_10"

rr_path() { echo "$BG/$1/$RD/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_$2.calibrated${3:+.out_of_fold}.ht"; }
submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

# Snapshot the scripts so later edits can't change this run.
S=$(mktemp -d /tmp/rr_bgs.XXXX)
cp "$F/cluster/annotate_cells_rr.py" "$F/cluster/run_full_run_scores.py" "$F/full_run_scores.py" \
  "$F/remap_rr_to_shet.py" "$F/rr_shrinkage_experiments.py" "$S/"
(cd "$S" && zip -q mods.zip full_run_scores.py remap_rr_to_shet.py rr_shrinkage_experiments.py)

TABLES=(); NAMES=()
for arm in $ARMS; do
  for cut in $CUTS; do
    until gsutil -q stat "$(rr_path "$arm" "$cut" oof)/_SUCCESS"; do sleep 300; done
    TABLES+=(--rr-tables "${arm}_$cut" "$(rr_path "$arm" "$cut")" --rr-tables "${arm}_${cut}_oof" "$(rr_path "$arm" "$cut" oof)")
    NAMES+=("${arm}_$cut")
  done
done
echo "=== every arm's rr is written"

gsutil -q stat "$OUT_CELLS" ||
  submit annotate_cells_bg "$S/annotate_cells_rr.py" -- --cells "$CELLS" "${TABLES[@]}" --out "$OUT_CELLS" || exit 1
gsutil -q stat "$OUT/rr_bg_pipeline.metrics.tsv" ||
  submit scores_bg "$S/run_full_run_scores.py" --pyfiles "$S/mods.zip" -- --cells "$OUT_CELLS" --out "$OUT" --prefix rr_bg_pipeline --rr-names "${NAMES[@]}" || exit 1
for f in $(gsutil ls "$OUT/rr_bg_pipeline.*"); do
  [ -e "$F/$(basename "$f")" ] || gsutil -q cp "$f" "$F/"
done
echo "=== $(date '+%H:%M:%S') all done"
