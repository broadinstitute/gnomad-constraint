"""
Test script to verify PAE and pLDDT filtering parameters work with test data.

This script loads the test datasets and runs determine_regions_with_min_oe_upper
with different parameter combinations to verify they work correctly.
"""

import hail as hl
from gnomad_constraint.experimental.proemis3d.utils import determine_regions_with_min_oe_upper

# Initialize Hail
hl.init(tmp_dir="gs://gnomad-tmp-4day")

# Load test datasets
# Use GCS path for cloud environments, or local path for local testing
import sys
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
"""
# Test: No filtering
print("\n=== Test: No PAE or pLDDT filtering ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=None,
    plddt_ht=None,
    max_pae=None,
    min_plddt=None,
    plddt_cutoff_method=None,
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: PAE truncate method
print("\n=== Test: PAE truncate_on_pairwise_pae_with_center (max_pae=15) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=None,
    max_pae=15.0,
    pae_cutoff_method="truncate_on_pairwise_pae_with_center",
    min_plddt=None,
    plddt_cutoff_method=None,
    debug=True,
)

print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: PAE filter method
print("\n=== Test 3: PAE filter_on_pairwise_pae_with_center (max_pae=15) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=None,
    max_pae=15.0,
    pae_cutoff_method="filter_on_pairwise_pae_with_center",
    min_plddt=None,
    plddt_cutoff_method=None,
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

print("\n=== Test: PAE filter_on_pairwise_pae_in_region (max_pae=15) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=None,
    max_pae=15.0,
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
    min_plddt=None,
    plddt_cutoff_method=None,
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: pLDDT truncate method
print("\n=== Test: pLDDT truncate_at_first_low_plddt (min_plddt=70) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=None,
    plddt_ht=plddt_ht,
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="truncate_at_first_low_plddt",
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: pLDDT remove method
print("\n=== Test: pLDDT remove_low_plddt_residues (min_plddt=70) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=None,
    plddt_ht=plddt_ht,
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="remove_low_plddt_residues",
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: pLDDT exclude method
print("\n=== Test: pLDDT exclude_low_plddt_from_stats (min_plddt=70) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=None,
    plddt_ht=plddt_ht,
    max_pae=None,
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    debug=True,
)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: Combined PAE (filter_on_pairwise_pae_in_region) and pLDDT (remove_low_plddt_residues)
print("\n=== Test: Combined PAE and pLDDT filtering ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=plddt_ht,
    max_pae=15.0,
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
    min_plddt=70.0,
    plddt_cutoff_method="remove_low_plddt_residues",
    debug=True,
)

print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Test: Combined PAE (truncate_on_pairwise_pae_with_center) and pLDDT (exclude_low_plddt_from_stats)
print("\n=== Test: Combined PAE (truncate_on_pairwise_pae_with_center) and pLDDT (exclude_low_plddt_from_stats) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=plddt_ht,
    max_pae=15.0,
    pae_cutoff_method="truncate_on_pairwise_pae_with_center",
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    debug=True,
)

print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")


# Test: Combined PAE (filter_on_pairwise_pae_with_center) and pLDDT (exclude_low_plddt_from_stats)
print("\n=== Test: Combined PAE (filter_on_pairwise_pae_with_center) and pLDDT (exclude_low_plddt_from_stats) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=plddt_ht,
    max_pae=15.0,
    pae_cutoff_method="filter_on_pairwise_pae_with_center",
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    debug=True,
)

print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
"""
# Test: Combined PAE (filter_on_pairwise_pae_in_region) and pLDDT (exclude_low_plddt_from_stats)
print("\n=== Test: Combined PAE (filter_on_pairwise_pae_in_region) and pLDDT (exclude_low_plddt_from_stats) ===")
result = determine_regions_with_min_oe_upper(
    af2_ht=af2_ht,
    oe_codon_ht=oe_codon_ht,
    pae_ht=pae_ht,
    plddt_ht=plddt_ht,
    max_pae=15.0,
    pae_cutoff_method="filter_on_pairwise_pae_in_region",
    min_plddt=70.0,
    plddt_cutoff_method="exclude_low_plddt_from_stats",
    debug=True,
)


print("\nAll tests completed!")

hl.stop()

