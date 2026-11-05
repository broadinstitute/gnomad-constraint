#!/bin/bash
# Home-built Roulette-style rates on jg1: build per-variant rate multipliers from
# gnomAD v4.1 genome rare variants (cluster/build_homebuilt_rates.py), fit the ridge rr
# on top of them (fixed offset), and run the constraint pipeline downstream with both.
# Reuses the full run's DNM and feature Tables (copied into its own root). Every output
# goes to a new root; nothing uses --overwrite.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg1}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_homebuilt.zip
A=gs://gnomad-tmp/jg/rr_full
H=gs://gnomad-tmp/jg/rr_homebuilt
OUT=gs://gnomad-tmp/jg/rr_downstream
CUT="ridge ridge_dnm_10"
RD=gnomad_v4.1.1/constraint/regional_rate
RATES=$H/rates/homebuilt_variant_rates.ht

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

source "$F/cluster/snapshot_code.sh"

# A quick end-to-end check on chr22, then the autosomes.
gsutil -q stat "$H/rates_test/homebuilt_variant_rates.ht/_SUCCESS" ||
  submit homebuilt_rates_test "$S/build_homebuilt_rates.py" --count-reference --count-variants --build-rates --contigs chr22 --out "$H/rates_test" || exit 1
gsutil -q stat "$RATES/_SUCCESS" ||
  submit homebuilt_rates "$S/build_homebuilt_rates.py" --count-reference --count-variants --build-rates --out "$H/rates" || exit 1

echo "=== copying the DNM and feature Tables to $H"
for ht in syn_dnm genomic_features13_1kb; do
  gsutil -q stat "$H/$RD/gnomad.v4.1.1.$ht.ht/_SUCCESS" || gsutil -q -m cp -r "$A/$RD/gnomad.v4.1.1.$ht.ht" "$H/$RD/" || exit 1
done
gsutil -q cp -n "$A/$RD/gnomad.v4.1.1.syn_dnm.tsv.bgz" "$H/$RD/"

RATE_ARGS=(--use-variant-rates --variant-rate-ht-path "$RATES")
submit homebuilt_training "$P/regional_rate_pipeline.py" --create-rr-training-set "${RATE_ARGS[@]}" --output-root "$H" &&
submit homebuilt_train "$P/regional_rate_pipeline.py" --train-rr-models --feature-pval-cutoffs $CUT --cv-chrom-split --output-root "$H" &&
submit homebuilt_rr "$P/regional_rate_pipeline.py" --compute-rr --apply-rr-safeguards --out-of-fold --feature-pval-cutoffs $CUT "${RATE_ARGS[@]}" --output-root "$H" || exit 1

# The clusters run Hail 0.2.134, so stop at the constraint groups (metrics need 0.2.137+).
submit constraint_homebuilt_ridge "$P/constraint_pipeline.py" --aggregate-per-variant-expected \
  --regional-rr-ht-path "$H/$RD/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_ridge.calibrated.ht" --variant-rate-ht-path "$RATES" \
  --aggregate-by-constraint-groups --output-root "$OUT/homebuilt_ridge" || exit 1
echo "=== $(date '+%H:%M:%S') all done"
