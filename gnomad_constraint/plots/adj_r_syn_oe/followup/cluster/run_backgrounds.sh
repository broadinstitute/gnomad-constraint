#!/bin/bash
# Proper pipeline runs of the rr training background designs on jg2. Every arm reuses
# the DNM and feature Tables from the full run (copied into its own root) and writes
# only to its own new root under gs://gnomad-tmp/jg/rr_bg/. Nothing uses --overwrite.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=jg2
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_bg.zip
P=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/pipeline
A=gs://gnomad-tmp/jg/rr_full
OUT=gs://gnomad-tmp/jg/rr_bg
CUT="bonferroni ridge ridge_dnm_3 ridge_dnm_10 ridge_dnm_30"

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit $CLUSTER "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}
ht_of() { gsutil ls -d "$A/**$1*.ht" 2>/dev/null | sed 's#/$##' | grep -v '\.ht/' | head -1; }
wait_for() {  # wait_for NAME_PATTERN
  local ht=""
  while :; do ht=$(ht_of "$1"); [ -n "$ht" ] && gsutil -q stat "$ht/_SUCCESS" && break; sleep 60; done
  echo "$ht"
}
copy_inputs() {  # copy_inputs ROOT HT...
  local root=$1; shift
  for ht in "$@"; do gsutil -q -m cp -r "$ht" "$(dirname "${ht/$A/$root}")/" || return 1; done
}
train_and_rr() {  # train_and_rr NAME ROOT
  submit "$1_train" "$P/regional_rate_pipeline.py" --train-rr-models --feature-pval-cutoffs $CUT --cv-chrom-split --output-root "$2" &&
  submit "$1_rr" "$P/regional_rate_pipeline.py" --compute-rr --apply-rr-safeguards --out-of-fold --feature-pval-cutoffs $CUT --output-root "$2"
}
arm() {  # arm NAME OPTIONS...
  local name=$1; shift; local root=$OUT/$name
  copy_inputs "$root" "$DNM" "$FEAT" || return 1
  submit "${name}_training" "$P/regional_rate_pipeline.py" --create-rr-training-set --output-root "$root" "$@" && train_and_rr "$name" "$root"
}

while [ "$(gcloud dataproc clusters describe $CLUSTER --format='value(status.state)' 2>/dev/null)" != RUNNING ]; do sleep 30; done
echo "=== $CLUSTER is running"
rm -f "$ZIP"
(cd /Users/jgoodric/PycharmProjects/gnomad-constraint && zip -rq "$ZIP" gnomad_constraint/ -x 'gnomad_constraint/experimental/*' 'gnomad_constraint/plots/*' 'gnomad_constraint/flowchart/*' 'gnomad_constraint/renv/*' '*.pyc' '*__pycache__*' '*.RData' '*.DS_Store' '*.png' '*.pdf' '*.ipynb' '*.Rproj.user*')
(cd /Users/jgoodric/PycharmProjects/gnomad_qc && zip -rq "$ZIP" gnomad_qc/ -x '*.pyc' '*__pycache__*' '*.DS_Store' '*.ipynb')
(cd /Users/jgoodric/PycharmProjects/gnomad_methods && zip -rq "$ZIP" gnomad/ -x '*.pyc' '*__pycache__*' '*.DS_Store')
submit bg_install "$F/cluster/install_driver_packages.py" || exit 1
submit bg_test "$F/cluster/test_background_options.py" || exit 1
tr '\r' '\n' < "$D/logs/bg_test.log" | grep -E "^variants|background rows"

echo "=== waiting for the full run's DNM and feature Tables"
DNM=$(wait_for syn_dnm); FEAT=$(wait_for genomic_features)
echo "    $DNM"; echo "    $FEAT"
for f in $(gsutil ls "$A/**syn_dnm*.tsv*" 2>/dev/null); do for n in all_locus absent_variant all_variant ctx_proportional ctx_even absent_locus; do gsutil -q cp "$f" "${f/$A/$OUT/$n}"; done; done

( arm all_locus --background-sites all &&
  arm absent_variant --background-unit variant &&
  arm ctx_proportional --background-context-sampling proportional ) & P1=$!
( arm all_variant --background-sites all --background-unit variant &&
  arm ctx_even --background-context-sampling even &&
  { TRAIN=$(wait_for rr_training); echo "    $TRAIN"; copy_inputs "$OUT/absent_locus" "$DNM" "$FEAT" "$TRAIN" && train_and_rr absent_locus "$OUT/absent_locus"; } ) & P2=$!
wait $P1; R1=$?; wait $P2; R2=$?
echo "=== $(date '+%H:%M:%S') finished: chain 1 rc=$R1, chain 2 rc=$R2"
gsutil ls -d "$OUT/**rr_1kb*.ht" 2>/dev/null | grep -v '\.ht/' | sort -u
