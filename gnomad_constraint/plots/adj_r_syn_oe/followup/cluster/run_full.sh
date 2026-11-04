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

# gsutil's ** only works as a whole path segment, so list the regional_rate folder.
hts() { gsutil ls -d "$1/*/constraint/regional_rate/*.ht" 2>/dev/null | sed 's#/$##'; }
done_ht() { hts "$1" | grep -q "$2" && gsutil -q stat "$(hts "$1" | grep "$2" | head -1)/_SUCCESS"; }

# The utilities come from the zip, so rebuild it from the current code.
rm -f "$ZIP"
(cd /Users/jgoodric/PycharmProjects/gnomad-constraint && zip -rq "$ZIP" gnomad_constraint/ -x 'gnomad_constraint/experimental/*' 'gnomad_constraint/plots/*' 'gnomad_constraint/flowchart/*' 'gnomad_constraint/renv/*' '*.pyc' '*__pycache__*' '*.RData' '*.DS_Store' '*.png' '*.pdf' '*.ipynb' '*.Rproj.user*')
(cd /Users/jgoodric/PycharmProjects/gnomad_qc && zip -rq "$ZIP" gnomad_qc/ -x '*.pyc' '*__pycache__*' '*.DS_Store' '*.ipynb')
(cd /Users/jgoodric/PycharmProjects/gnomad_methods && zip -rq "$ZIP" gnomad/ -x '*.pyc' '*__pycache__*' '*.DS_Store')

# Steps 1 and 2 are skipped when an earlier attempt already wrote their Tables.
done_ht "$A" syn_dnm.ht || run A_1_dnms prepare_syn_dnms.py --prepare-syn-dnms --output-root "$A" || exit 1
done_ht "$A" genomic_features || run A_2_features regional_rate_pipeline.py --import-genomic-features --output-root "$A" || exit 1

echo "=== copying the DNM and feature Tables to $B"
for ht in $(hts "$A" | grep -E 'syn_dnm|genomic_features'); do
  dst="${ht/$A/$B}"
  if gsutil -q stat "$dst/_SUCCESS"; then echo "    $dst already there"; continue; fi
  echo "    $ht -> $dst"
  gsutil -q -m cp -r "$ht" "$(dirname "$dst")/" || exit 1
done
for f in $(gsutil ls "$A/*/constraint/regional_rate/*syn_dnm*.tsv*" 2>/dev/null); do gsutil -q cp -n "$f" "${f/$A/$B}"; done

chain "$A" A & PA=$!
chain "$B" B --exclude-qc-flagged-genes & PB=$!
wait $PA; RA=$?; wait $PB; RB=$?
echo "=== $(date '+%H:%M:%S') finished: baseline rc=$RA, QC-excluded rc=$RB"
{ hts "$A"; hts "$B"; } | grep rr_1kb
