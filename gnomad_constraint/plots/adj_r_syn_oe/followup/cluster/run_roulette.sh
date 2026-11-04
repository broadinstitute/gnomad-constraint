#!/bin/bash
# Roulette run on jg2: import Roulette's per-variant rates over coding sites, fit the
# ridge rr on top of them (log multiplier as a fixed offset), and run the constraint
# pipeline downstream with both. Reuses the full run's DNM and feature Tables (copied
# into its own root). Every output goes to a new root; nothing uses --overwrite.
set -uo pipefail
export PATH="/Users/jgoodric/google-cloud-sdk/bin:/usr/bin:/bin:$PATH"
export CLOUDSDK_PYTHON=/Users/jgoodric/miniconda3/bin/python CLOUDSDK_CORE_PROJECT=broad-mpg-gnomad CLOUDSDK_DATAPROC_REGION=us-central1
HAILCTL=/Users/jgoodric/miniconda3/envs/hail-134/bin/hailctl
CLUSTER=${CLUSTER:-jg2}
F=/Users/jgoodric/PycharmProjects/gnomad-constraint/gnomad_constraint/plots/adj_r_syn_oe/followup
D=$F/data/cluster
ZIP=$D/pyfiles_roulette.zip
A=gs://gnomad-tmp/jg/rr_full
R=gs://gnomad-tmp/jg/rr_roulette
T=gs://gnomad-tmp/jg/rr_roulette_test
OUT=gs://gnomad-tmp/jg/rr_downstream
VCF=gs://gnomad-tmp/jg/roulette/v5.2
CUT="ridge ridge_dnm_10"
RD=gnomad_v4.1.1/constraint/regional_rate

submit() {
  local tag=$1 script=$2; shift 2
  echo "=== $(date '+%H:%M:%S') start $tag"
  if ! "$HAILCTL" dataproc submit "$CLUSTER" "$script" --pyfiles "$ZIP" -- "$@" > "$D/logs/$tag.log" 2>&1; then
    echo "=== $(date '+%H:%M:%S') FAILED $tag"; tr '\r' '\n' < "$D/logs/$tag.log" | grep -v '^\[Stage' | tail -30; return 1
  fi
  echo "=== $(date '+%H:%M:%S') done $tag"
}

source "$F/cluster/snapshot_code.sh"

VCFS=()
for c in $(seq 1 22); do VCFS+=("$VCF/${c}_rate_v5.2_TFBS_correction_all.vcf.bgz"); done

# A quick end-to-end check on chr22 before the full import.
submit roulette_import_test "$P/regional_rate_pipeline.py" --import-variant-rates --roulette-vcf-paths "${VCFS[21]}" --output-root "$T" || exit 1
submit roulette_import "$P/regional_rate_pipeline.py" --import-variant-rates --roulette-vcf-paths "${VCFS[@]}" --output-root "$R" || exit 1

echo "=== copying the DNM and feature Tables to $R"
for ht in syn_dnm genomic_features13_1kb; do
  src="$A/$RD/gnomad.v4.1.1.$ht.ht"; dst="$R/$RD/gnomad.v4.1.1.$ht.ht"
  gsutil -q stat "$dst/_SUCCESS" || gsutil -q -m cp -r "$src" "$R/$RD/" || exit 1
done
gsutil -q cp -n "$A/$RD/gnomad.v4.1.1.syn_dnm.tsv.bgz" "$R/$RD/"

submit roulette_training "$P/regional_rate_pipeline.py" --create-rr-training-set --use-variant-rates --output-root "$R" &&
submit roulette_train "$P/regional_rate_pipeline.py" --train-rr-models --feature-pval-cutoffs $CUT --cv-chrom-split --output-root "$R" &&
submit roulette_rr "$P/regional_rate_pipeline.py" --compute-rr --apply-rr-safeguards --out-of-fold --feature-pval-cutoffs $CUT --use-variant-rates --output-root "$R" || exit 1

RATES="$R/$RD/gnomad.v4.1.1.roulette_variant_rates.ht"
submit constraint_roulette_ridge "$P/constraint_pipeline.py" --aggregate-per-variant-expected \
  --regional-rr-ht-path "$R/$RD/gnomad.v4.1.1.syn_dnm_rr_1kb.fs_ridge.calibrated.ht" --variant-rate-ht-path "$RATES" \
  --aggregate-by-constraint-groups --compute-constraint-metrics --prepare-release --export-release-tsv --output-root "$OUT/roulette_ridge" || exit 1
echo "=== $(date '+%H:%M:%S') all done"
gsutil ls "$OUT/roulette_ridge/**" 2>/dev/null | grep "\.tsv"
