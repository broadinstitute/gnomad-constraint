"""
Test script to verify PAE and pLDDT filtering parameters work with test data.

This script loads the test datasets and runs determine_regions_with_min_oe_upper
with different parameter combinations to verify they work correctly.
"""

import sys

import hail as hl

from gnomad_constraint.experimental.proemis3d.utils import (
    determine_regions_with_min_oe_upper,
)

# Initialize Hail
hl.init(tmp_dir="gs://gnomad-tmp-4day")

# Load test datasets
# Use GCS path for cloud environments, or local path for local testing

if "--local" in sys.argv:
    test_data_dir = "test_data"
else:
    # Default to GCS path
    test_data_dir = "gs://gnomad-tmp-4day/proemis3d_test_data"

af2_ht = hl.read_table(f"{test_data_dir}/af2_test.ht")
oe_codon_ht = hl.read_table(f"{test_data_dir}/oe_codon_test.ht")
pae_ht = hl.read_table(f"{test_data_dir}/pae_test.ht")
plddt_ht = hl.read_table(f"{test_data_dir}/plddt_test.ht")

print("Loaded test datasets:")
print(f"  af2_ht: {af2_ht.count()} rows")
print(f"  oe_codon_ht: {oe_codon_ht.count()} rows")
print(f"  pae_ht: {pae_ht.count()} rows")
print(f"  plddt_ht: {plddt_ht.count()} rows")


def test_determine_regions_with_min_oe_upper(
    max_pae=None,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method=None,
):
    print(
        f"\n"
        "=============================================================================\n"
        "=== Test:                                                            ===\n"
        "===    max_pae={max_pae:>5}                                          ===\n"
        "===    min_plddt={min_plddt:>5}                                          ===\n"
        "===    plddt_cutoff_method={plddt_cutoff_method:>50}===\n"
        "===    pae_cutoff_method={pae_cutoff_method:>52}===\n"
        "=============================================================================\n"
    )
    result = determine_regions_with_min_oe_upper(
        af2_ht=af2_ht,
        oe_codon_ht=oe_codon_ht,
        pae_ht=pae_ht if max_pae is not None else None,
        plddt_ht=plddt_ht if min_plddt is not None else None,
        max_pae=max_pae,
        min_plddt=min_plddt,
        plddt_cutoff_method=plddt_cutoff_method,
        pae_cutoff_method=pae_cutoff_method,
        debug=True,
    )
    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

    return result


# Test: No filtering
test_determine_regions_with_min_oe_upper(
    max_pae=None,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method=None,
)

# Test: PAE truncate method
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method="truncate_on_pairwise_pae_with_center",
)

# Test: PAE filter method
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method="filter_on_pairwise_pae_with_center",
)

# Test: PAE filter_on_pairwise_pae_in_region
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
)

# Test: pLDDT truncate method
test_determine_regions_with_min_oe_upper(
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="truncate_at_first_low_plddt",
    pae_cutoff_method=None,
)

# Test: pLDDT remove method
test_determine_regions_with_min_oe_upper(
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="remove_low_plddt_residues",
    pae_cutoff_method=None,
)

# Test: pLDDT exclude method
test_determine_regions_with_min_oe_upper(
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    pae_cutoff_method=None,
)

# Test: Combined PAE (filter_on_pairwise_pae_in_region) and pLDDT
# (remove_low_plddt_residues)
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=70.0,
    plddt_cutoff_method="remove_low_plddt_residues",
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
)

# Test: Combined PAE (truncate_on_pairwise_pae_with_center) and pLDDT
# (exclude_low_plddt_from_stats)
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    pae_cutoff_method="truncate_on_pairwise_pae_with_center",
)

# Test: Combined PAE (filter_on_pairwise_pae_with_center) and pLDDT
# (exclude_low_plddt_from_stats)
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    pae_cutoff_method="filter_on_pairwise_pae_with_center",
)

# Test: Combined PAE (filter_on_pairwise_pae_in_region) and pLDDT
# (exclude_low_plddt_from_stats)
test_determine_regions_with_min_oe_upper(
    max_pae=15.0,
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
)


print("\nAll tests completed!")

hl.stop()
