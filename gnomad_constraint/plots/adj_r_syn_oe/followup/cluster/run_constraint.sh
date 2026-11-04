#!/bin/bash
# Score the full-run rr tables downstream: annotate the exported synonymous cells with
# every rr, then run the constraint pipeline's post-calibration steps once per rr
# (aggregation with rr, constraint groups, metrics, release TSV). Every output goes to
# a new root under gs://gnomad-tmp/jg/; nothing is run with --overwrite.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_constraint.zip
source "$F/cluster/snapshot_code.sh"
A=gs://gnomad-tmp/jg/rr_full
B=gs://gnomad-tmp/jg/rr_full_qcexcl
OUT=gs://gnomad-tmp/jg/rr_downstream
LEG=gs://gnomad/v4.1/constraint/resources/annotations/ht
CELLS=gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export/syn_cells.tsv.bgz

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit jg1 "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

rr_path() {  # rr_path ROOT CUTOFF [oof]
  gsutil ls -d "$1/*/constraint/regional_rate/*rr_1kb.fs_$2.*ht" 2>/dev/null | sed 's#/$##' | sort -u |
    if [ "${3:-}" = oof ]; then grep -i 'oof\|out_of_fold'; else grep -v -i 'oof\|out_of_fold'; fi | head -1
}

TABLES=(--rr-tables adj_r "$LEG/adj_r_per_context_methyl_genome_1kb_autosome.agg.ht"
        --rr-tables syn_adj_r "$LEG/adj_r_syn_dnm_per_context_methyl_genome_1kb_autosome.agg.ht")
for run in full qcexcl; do
  root=$A; [ $run = qcexcl ] && root=$B
  for cut in bonferroni 1e-2 ridge; do
    TABLES+=(--rr-tables "${cut}_${run}" "$(rr_path $root $cut)" --rr-tables "${cut}_${run}_oof" "$(rr_path $root $cut oof)")
  done
done
printf '%s\n' "${TABLES[@]}" | paste - - - | sed 's/^/    /'
gsutil -q stat "$OUT/syn_cells_rr_full_run.tsv.bgz" && echo "=== annotate_cells already done" ||
submit annotate_cells $S/annotate_cells_rr.py --cells "$CELLS" "${TABLES[@]}" --out "$OUT/syn_cells_rr_full_run.tsv.bgz" || exit 1

for v in ridge_qcexcl bonferroni_full ridge_full bonferroni_qcexcl; do
  cut=${v%_*}; run=${v##*_}; root=$A; [ $run = qcexcl ] && root=$B
  if gsutil ls "$OUT/$v/**" 2>/dev/null | grep -q "\.tsv"; then echo "=== constraint_$v already done"; continue; fi
  submit "constraint_$v" "$P/constraint_pipeline.py" --aggregate-per-variant-expected --regional-rr-ht-path "$(rr_path $root $cut)" \
    --aggregate-by-constraint-groups --compute-constraint-metrics --prepare-release --export-release-tsv --output-root "$OUT/$v" || exit 1
done
echo "=== $(date '+%H:%M:%S') all done"
gsutil ls "$OUT/*/*/constraint/**" 2>/dev/null | grep "\.tsv"
