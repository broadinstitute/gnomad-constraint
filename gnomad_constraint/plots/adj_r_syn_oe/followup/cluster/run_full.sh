#!/bin/bash
# Full rr pipeline run: baseline (A) and with QC-flagged genes left out of training (B).
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
C=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup/data/cluster
ZIP=$C/pyfiles.zip
P=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/pipeline
A=gs://gnomad-tmp/jg/rr_full
B=gs://gnomad-tmp/jg/rr_full_qcexcl
CUT="bonferroni 1e-2 ridge"

run() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit jg1 "$P/$script" --pyfiles "$ZIP" -- "$@" > "$C/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$C/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

chain() {
  local root=$1 name=$2; shift 2
  run "${name}_3_training" regional_rate_pipeline.py --create-rr-training-set --output-root "$root" "$@" &&
  run "${name}_4_train" regional_rate_pipeline.py --train-rr-models --feature-pval-cutoffs $CUT --cv-chrom-split --output-root "$root" &&
  run "${name}_5_rr" regional_rate_pipeline.py --compute-rr --apply-rr-safeguards --out-of-fold --feature-pval-cutoffs $CUT --output-root "$root"
}

run A_1_dnms prepare_syn_dnms.py --prepare-syn-dnms --output-root "$A" || exit 1
run A_2_features regional_rate_pipeline.py --import-genomic-features --output-root "$A" || exit 1

echo "=== copying the DNM and feature Tables to $B"
for ht in $(gsutil ls -d "$A/**.ht" 2>/dev/null | grep -E 'syn_dnm|genomic_features' | sed 's#/$##' | sort -u); do
  case "$ht" in *.ht/*) continue;; esac
  dst="${ht/$A/$B}"
  echo "    $ht -> $dst"
  gsutil -q -m cp -r "$ht" "$(dirname "$dst")/" || exit 1
done
for f in $(gsutil ls "$A/**syn_dnm*.tsv*" 2>/dev/null); do gsutil -q cp "$f" "${f/$A/$B}"; done

chain "$A" A & PA=$!
chain "$B" B --exclude-qc-flagged-genes & PB=$!
wait $PA; RA=$?; wait $PB; RB=$?
echo "=== $(date '+%H:%M:%S') finished: baseline rc=$RA, QC-excluded rc=$RB"
gsutil ls -d "$A/**rr_1kb*.ht" "$B/**rr_1kb*.ht" 2>/dev/null | grep -v '\.ht/' | sort -u
