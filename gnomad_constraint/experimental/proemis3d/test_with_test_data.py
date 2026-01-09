"""
Test script to verify PAE and pLDDT filtering parameters work with test data.

This script can create test datasets and/or run tests with determine_regions_with_min_oe_upper
and optionally run_forward using different parameter combinations to verify they work correctly.

Usage:
    # Create test datasets
    python test_with_test_data.py --create-test-data
    python test_with_test_data.py --create-test-data --local  # Create locally

    # Run all tests
    python test_with_test_data.py

    # Run all tests including run_forward
    python test_with_test_data.py --run-forward

    # Run a specific test by name
    python test_with_test_data.py --test no_filtering
    python test_with_test_data.py --test pae_truncate

    # Run a specific test with run_forward
    python test_with_test_data.py --test no_filtering --run-forward

    # Run a group of tests
    python test_with_test_data.py --group plddt_methods
    python test_with_test_data.py --group pae_methods
    python test_with_test_data.py --group combined

    # Use local test data
    python test_with_test_data.py --local

    # List available tests
    python test_with_test_data.py --list-tests

    # Run with a specific uniprot/transcript from production data
    python test_with_test_data.py --uniprot-id P12345 --transcript-id ENST00000123456
    python test_with_test_data.py --uniprot-id P12345 --transcript-id ENST00000123456 --run-forward
"""

import argparse
import os
import sys

import hail as hl

from gnomad_constraint.experimental.proemis3d.resources import (
    get_af2_dist_ht,
    get_af2_pae_ht,
    get_af2_plddt_ht,
    get_gencode_pos_ht,
    get_obs_exp_ht,
)
from gnomad_constraint.experimental.proemis3d.utils import (
    determine_regions_with_min_oe_upper,
    generate_codon_oe_table,
    run_forward,
)

# ANSI formatting codes.
RESET = "\033[0m"
BOLD = "\033[1m"

# Test data creation constants.
NUM_RESIDUES = 15
TEST_UNIPROT_ID = "P12345"
TEST_TRANSCRIPT_ID = "ENST00000123456"

# Define all available tests.
ALL_TESTS = {
    "no_filtering": {
        "max_pae": None,
        "min_plddt": None,
        "plddt_cutoff_method": None,
        "pae_cutoff_method": None,
    },
    "pae_truncate": {
        "max_pae": 15.0,
        "min_plddt": None,
        "plddt_cutoff_method": None,
        "pae_cutoff_method": "truncate_on_pairwise_pae_with_center",
    },
    "pae_filter_center": {
        "max_pae": 15.0,
        "min_plddt": None,
        "plddt_cutoff_method": None,
        "pae_cutoff_method": "filter_on_pairwise_pae_with_center",
    },
    "pae_filter_region": {
        "max_pae": 15.0,
        "min_plddt": None,
        "plddt_cutoff_method": None,
        "pae_cutoff_method": "filter_on_pairwise_pae_in_region",
    },
    "plddt_truncate": {
        "max_pae": None,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "truncate_at_first_low_plddt",
        "pae_cutoff_method": None,
    },
    "plddt_remove": {
        "max_pae": None,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "remove_low_plddt_residues",
        "pae_cutoff_method": None,
    },
    "plddt_exclude": {
        "max_pae": None,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "exclude_low_plddt_from_stats",
        "pae_cutoff_method": None,
    },
    "combined_pae_region_plddt_remove": {
        "max_pae": 15.0,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "remove_low_plddt_residues",
        "pae_cutoff_method": "filter_on_pairwise_pae_in_region",
    },
    "combined_pae_truncate_plddt_exclude": {
        "max_pae": 15.0,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "exclude_low_plddt_from_stats",
        "pae_cutoff_method": "truncate_on_pairwise_pae_with_center",
    },
    "combined_pae_filter_center_plddt_exclude": {
        "max_pae": 15.0,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "exclude_low_plddt_from_stats",
        "pae_cutoff_method": "filter_on_pairwise_pae_with_center",
    },
    "combined_pae_filter_region_plddt_exclude": {
        "max_pae": 15.0,
        "min_plddt": 70.0,
        "plddt_cutoff_method": "exclude_low_plddt_from_stats",
        "pae_cutoff_method": "filter_on_pairwise_pae_in_region",
    },
}

# Define test groups.
TEST_GROUPS = {
    "all": list(ALL_TESTS.keys()),
    "plddt_methods": [
        "plddt_truncate",
        "plddt_remove",
        "plddt_exclude",
    ],
    "pae_methods": [
        "pae_truncate",
        "pae_filter_center",
        "pae_filter_region",
    ],
    "combined": [
        "combined_pae_region_plddt_remove",
        "combined_pae_truncate_plddt_exclude",
        "combined_pae_filter_center_plddt_exclude",
        "combined_pae_filter_region_plddt_exclude",
    ],
    "basic": [
        "no_filtering",
        "pae_truncate",
        "plddt_truncate",
    ],
}


def test_determine_regions_with_min_oe_upper(
    af2_ht,
    oe_codon_ht,
    pae_ht,
    plddt_ht,
    max_pae=None,
    min_plddt=None,
    plddt_cutoff_method=None,
    pae_cutoff_method=None,
):
    """Run a single test with the given parameters."""
    # Format values, handling None.
    max_pae_str = str(max_pae) if max_pae is not None else "None"
    min_plddt_str = str(min_plddt) if min_plddt is not None else "None"
    plddt_method_str = (
        str(plddt_cutoff_method) if plddt_cutoff_method is not None else "None"
    )
    pae_method_str = str(pae_cutoff_method) if pae_cutoff_method is not None else "None"

    print(
        "\n"
        + f"{BOLD}={RESET}" * 77
        + "\n"
        + f"{BOLD}=== Test:{RESET}"
        + " " * 65
        + f"{BOLD}==={RESET}\n"
        + f"{BOLD}===    max_pae: {max_pae_str:<58}==={RESET}\n"
        + f"{BOLD}===    min_plddt: {min_plddt_str:<56}==={RESET}\n"
        + f"{BOLD}===    plddt_cutoff_method: {plddt_method_str:<46}==={RESET}\n"
        + f"{BOLD}===    pae_cutoff_method: {pae_method_str:<48}==={RESET}\n"
        + f"{BOLD}={RESET}" * 77
        + "\n"
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
    print("\n\n\n")

    return result


def test_run_forward(
    regions_ht,
    min_exp_mis=16,
    oe_upper_method="gamma",
    model_comparison_method="aic",
    lrt_alpha=0.001,
    lrt_df_added=1,
    bonferroni_per_round=True,
    aic_weight_thresh=0.80,
):
    """Run run_forward test with the given parameters."""
    print(
        "\n"
        + f"{BOLD}={RESET}" * 77
        + "\n"
        + f"{BOLD}=== Testing run_forward with parameters:{RESET}"
        + " " * 34
        + f"{BOLD}==={RESET}\n"
        + f"{BOLD}===    min_exp_mis: {min_exp_mis:<54}==={RESET}\n"
        + f"{BOLD}===    oe_upper_method: {oe_upper_method:<50}==={RESET}\n"
        + f"{BOLD}===    model_comparison_method: {model_comparison_method:<42}==={RESET}\n"
        + f"{BOLD}===    lrt_alpha: {lrt_alpha:<56}==={RESET}\n"
        + f"{BOLD}===    lrt_df_added: {lrt_df_added:<53}==={RESET}\n"
        + f"{BOLD}===    bonferroni_per_round: {bonferroni_per_round:<45}==={RESET}\n"
        + f"{BOLD}===    aic_weight_thresh: {aic_weight_thresh:<48}==={RESET}\n"
        + f"{BOLD}={RESET}" * 77
        + "\n"
    )

    # Prepare the table for run_forward
    # run_forward expects:
    # - oe: array of structs with obs and exp (without residue_index)
    # - min_oe_upper: array of structs with region information
    # determine_regions_with_min_oe_upper returns oe with residue_index, so we need to remove it
    # The table should already be keyed by uniprot_id and transcript_id from group_by
    regions_ht_prepared = regions_ht.annotate(
        oe=regions_ht.oe.map(lambda x: x.select("obs", "exp"))
    )

    result = run_forward(
        regions_ht_prepared,
        min_exp_mis=min_exp_mis,
        oe_upper_method=oe_upper_method,
        model_comparison_method=model_comparison_method,
        lrt_alpha=lrt_alpha,
        lrt_df_added=lrt_df_added,
        bonferroni_per_round=bonferroni_per_round,
        aic_weight_thresh=aic_weight_thresh,
        debug=True,
    )
    print("\n\n\n")

    return result


def create_symmetric_dist_matrix(num_residues: int) -> dict:
    """Create a symmetric distance matrix where dist[i][j] = dist[j][i].

    Self-distances (i==j) are always 0.0.
    """
    dist_matrix = {}
    for i in range(num_residues):
        dist_matrix[i] = {}
        for j in range(num_residues):
            if i == j:
                dist_matrix[i][j] = 0.0  # Self-distance is always 0
            elif j < i:
                # Use symmetric value from lower triangle (already computed)
                dist_matrix[i][j] = dist_matrix[j][i]
            else:
                # Create realistic distances (5-50 Angstroms)
                # Closer residues (in sequence) have smaller distances
                seq_dist = abs(i - j)
                # Add some 3D structure variation
                struct_variation = ((i + j) % 5) * 1.2
                dist = 5.0 + seq_dist * 2.0 + struct_variation
                # Ensure distance is never 0 for non-self pairs
                dist = max(1.0, dist)
                dist_matrix[i][j] = float(dist)
    return dist_matrix


def create_dist_mat(num_residues: int, center_idx: int, dist_matrix: dict) -> list:
    """Create a distance matrix for a residue from the symmetric matrix.

    Includes self-distance (0) and all other residues, in residue index order (not sorted).
    The distances are kept in their original order by residue_index.
    """
    dist_mat = []
    # Add all residues including self, in residue index order
    for i in range(num_residues):
        dist = dist_matrix[center_idx][i]
        # Verify self-distance is 0.0
        if i == center_idx and dist != 0.0:
            raise ValueError(
                f"Self-distance for residue {center_idx} is {dist}, should be 0.0"
            )
        # Verify non-self distances are > 0
        if i != center_idx and dist == 0.0:
            raise ValueError(
                f"Non-self distance from {center_idx} to {i} is 0.0, should be > 0"
            )
        dist_mat.append({"dist": float(dist), "residue_index": int(i)})
    # Do NOT sort - keep in residue index order
    return dist_mat


def create_pae_matrix(num_residues: int, dist_matrix: dict) -> dict:
    """Create a PAE matrix where PAE[i][j] is the error at residue i when aligned on residue j.

    PAE should be:
    - Low (0-5) for residues close in sequence (likely same domain)
    - Medium (5-15) for residues at moderate sequence distance
    - High (>15) for residues far in sequence (likely different domains)
    - Correlated with 3D distance but also sequence distance
    """
    pae_matrix = {}
    for i in range(num_residues):
        pae_matrix[i] = {}
        for j in range(num_residues):
            if i == j:
                pae_matrix[i][j] = 0.0  # Self-PAE is 0
            else:
                seq_dist = abs(i - j)
                struct_dist = dist_matrix[i][j]

                # Base PAE increases with sequence distance (different domains)
                # Residues far apart in sequence are more likely in different domains
                base_pae = seq_dist * 0.8

                # Add 3D structure influence (but less than sequence distance)
                # Residues far in 3D space also have higher PAE
                struct_influence = (
                    struct_dist - 5.0
                ) * 0.1  # Scale 3D distance contribution

                # Add some noise/variation
                noise = ((i + j) % 7) * 0.3

                pae = base_pae + struct_influence + noise

                # Create test pattern: residues 0-4 have high PAE with residues 10-14
                # (simulating different domains)
                if (i < 5 and j >= 10) or (i >= 10 and j < 5):
                    pae = 20.0 + seq_dist * 0.5  # High PAE for cross-domain pairs

                # Ensure minimum PAE for non-self pairs
                pae = max(1.0, pae)

                pae_matrix[i][j] = float(pae)

    return pae_matrix


def create_test_data(local: bool = False):
    """Create small test datasets for testing PAE and pLDDT filtering parameters.

    This function creates minimal test datasets that can be used to verify:
    - PAE cutoff methods (truncate_on_pairwise_pae_with_center, filter_on_pairwise_pae_with_center, filter_on_pairwise_pae_in_region)
    - pLDDT cutoff methods (truncate_at_first_low_plddt, remove_low_plddt_residues, exclude_low_plddt_from_stats)

    :param local: If True, write to local directory. If False, write to GCS.
    """
    # Initialize Hail
    hl.init(tmp_dir="gs://gnomad-tmp-4day")

    # Create symmetric distance matrix first
    symmetric_dist_matrix = create_symmetric_dist_matrix(NUM_RESIDUES)

    # Create AF2 HT (AlphaFold2 data)
    # Key: uniprot_id, aa_index
    # Fields: enst, dist_mat
    af2_dicts = []
    for aa_idx in range(NUM_RESIDUES):
        dist_mat_list = create_dist_mat(NUM_RESIDUES, aa_idx, symmetric_dist_matrix)
        af2_dicts.append(
            {
                "uniprot_id": TEST_UNIPROT_ID,
                "aa_index": aa_idx,
                "enst": TEST_TRANSCRIPT_ID,
                "dist_mat": dist_mat_list,
            }
        )

    # Use hl.Table.parallelize to create from Python data
    af2_ht = hl.Table.parallelize(
        [
            hl.struct(
                uniprot_id=d["uniprot_id"],
                aa_index=d["aa_index"],
                enst=d["enst"],
                dist_mat=hl.array(
                    [
                        hl.struct(
                            dist=hl.float32(x["dist"]),
                            residue_index=hl.int32(x["residue_index"]),
                        )
                        for x in d["dist_mat"]
                    ]
                ),
            )
            for d in af2_dicts
        ]
    )
    af2_ht = af2_ht.key_by("uniprot_id", "aa_index")

    # Create OE codon HT (Observed/Expected data)
    # Key: uniprot_id
    # Fields: oe_by_transcript (array of structs with enst and oe)
    #   where oe is an array of structs with obs and exp
    oe_array = []
    for aa_idx in range(NUM_RESIDUES):
        # Create realistic OE data
        # Some residues have low OE (constrained), some have high OE (tolerant)
        obs = 2 if aa_idx < 5 else 8  # First 5 residues are constrained
        exp = 10.0 + (aa_idx % 3) * 2.0
        oe_array.append(
            {
                "obs": int(obs),
                "exp": float(exp),
            }
        )

    # Structure: oe_by_transcript is array of {enst: str, oe: array<{obs:
    # int64, exp: float64}>}
    oe_codon_ht = hl.Table.parallelize(
        [
            hl.struct(
                uniprot_id=TEST_UNIPROT_ID,
                oe_by_transcript=hl.array(
                    [
                        hl.struct(
                            enst=TEST_TRANSCRIPT_ID,
                            oe=hl.array(
                                [
                                    hl.struct(
                                        obs=hl.int64(oe_item["obs"]),
                                        exp=hl.float64(oe_item["exp"]),
                                    )
                                    for oe_item in oe_array
                                ]
                            ),
                        )
                    ]
                ),
            )
        ]
    )
    oe_codon_ht = oe_codon_ht.key_by("uniprot_id")

    # Create PAE matrix
    pae_matrix = create_pae_matrix(NUM_RESIDUES, symmetric_dist_matrix)

    # Convert to the format needed for the HT
    pae_dicts = []
    for aa_idx in range(NUM_RESIDUES):
        pae_array = []
        for other_idx in range(NUM_RESIDUES):
            pae_array.append(pae_matrix[aa_idx][other_idx])
        pae_dicts.append(
            {
                "uniprot_id": TEST_UNIPROT_ID,
                "aa_index": aa_idx,
                "pae": pae_array,
            }
        )

    pae_ht = hl.Table.parallelize(
        [
            hl.struct(
                uniprot_id=d["uniprot_id"],
                aa_index=d["aa_index"],
                pae=hl.array([hl.float32(x) for x in d["pae"]]),
            )
            for d in pae_dicts
        ]
    )
    pae_ht = pae_ht.key_by("uniprot_id", "aa_index")

    # Create pLDDT HT (per-residue confidence data)
    # Key: uniprot_id, aa_index
    # Fields: plddt (pLDDT score)
    # For testing: some residues will have low pLDDT (<70) to test filtering
    plddt_dicts = []
    for aa_idx in range(NUM_RESIDUES):
        # Create pLDDT values: some low (<70) to test filtering
        # Residues 8-12 have low pLDDT
        if 8 <= aa_idx <= 12:
            plddt = 50.0 + (aa_idx % 3) * 5.0  # Low pLDDT (50-65)
        else:
            plddt = 80.0 + (aa_idx % 5) * 3.0  # High pLDDT (80-95)
        plddt_dicts.append(
            {
                "uniprot_id": TEST_UNIPROT_ID,
                "aa_index": aa_idx,
                "plddt": float(plddt),
            }
        )

    plddt_ht = hl.Table.parallelize(
        [
            hl.struct(
                uniprot_id=d["uniprot_id"],
                aa_index=d["aa_index"],
                plddt=hl.float32(d["plddt"]),
            )
            for d in plddt_dicts
        ]
    )
    plddt_ht = plddt_ht.key_by("uniprot_id", "aa_index")

    # Write test datasets
    # Use GCS path for cloud environments, or local path for local testing
    if local:
        output_dir = "test_data"
        os.makedirs(output_dir, exist_ok=True)
    else:
        # Default to GCS path
        output_dir = "gs://gnomad-tmp-4day/proemis3d_test_data"

    af2_ht.write(f"{output_dir}/af2_test.ht", overwrite=True)
    oe_codon_ht.write(f"{output_dir}/oe_codon_test.ht", overwrite=True)
    pae_ht.write(f"{output_dir}/pae_test.ht", overwrite=True)
    plddt_ht.write(f"{output_dir}/plddt_test.ht", overwrite=True)

    print(f"Test datasets created in {output_dir}/")
    print(f"  - af2_test.ht: {af2_ht.count()} rows")
    print(f"  - oe_codon_test.ht: {oe_codon_ht.count()} rows")
    print(f"  - pae_test.ht: {pae_ht.count()} rows")
    print(f"  - plddt_test.ht: {plddt_ht.count()} rows")

    # Print summary
    print("\nTest data summary:")
    print(f"  UniProt ID: {TEST_UNIPROT_ID}")
    print(f"  Transcript ID: {TEST_TRANSCRIPT_ID}")
    print(f"  Number of residues: {NUM_RESIDUES}")
    print("\nPAE test pattern:")
    print("  - Residues 0-4 have high PAE (>15) with residues 10-14")
    print("  - Other pairs have low PAE (<15)")
    print("\npLDDT test pattern:")
    print("  - Residues 8-12 have low pLDDT (50-65)")
    print("  - Other residues have high pLDDT (80-95)")

    hl.stop()


def list_tests():
    """Print all available tests and groups."""
    print("\n=== Available Tests ===")
    for test_name, test_params in ALL_TESTS.items():
        print(f"\n  {test_name}:")
        for key, value in test_params.items():
            print(f"    {key}: {value}")

    print("\n=== Available Test Groups ===")
    for group_name, test_names in TEST_GROUPS.items():
        print(f"\n  {group_name}:")
        for test_name in test_names:
            print(f"    - {test_name}")


def run_tests(
    test_names, af2_ht, oe_codon_ht, pae_ht, plddt_ht, run_forward_tests=False
):
    """Run the specified tests."""
    for test_name in test_names:
        if test_name not in ALL_TESTS:
            print(f"Warning: Unknown test '{test_name}', skipping...")
            continue

        test_params = ALL_TESTS[test_name]
        # Run determine_regions_with_min_oe_upper
        regions_ht = test_determine_regions_with_min_oe_upper(
            af2_ht=af2_ht,
            oe_codon_ht=oe_codon_ht,
            pae_ht=pae_ht,
            plddt_ht=plddt_ht,
            **test_params,
        )

        # Optionally run run_forward on the result
        if run_forward_tests and regions_ht is not None:
            # Test with different model comparison methods
            for model_method in ["aic", "lrt", "aic_weight"]:
                test_run_forward(
                    regions_ht,
                    min_exp_mis=16,
                    oe_upper_method="gamma",
                    model_comparison_method=model_method,
                    lrt_alpha=0.001,
                    lrt_df_added=1,
                    bonferroni_per_round=True,
                    aic_weight_thresh=0.80,
                )


def main(args):
    """Run tests with test data."""
    # Handle create-test-data flag
    if args.create_test_data:
        create_test_data(local=args.local)
        return

    # Initialize Hail
    hl.init(tmp_dir="gs://gnomad-tmp-4day")

    # Check if running with a specific uniprot/transcript
    if args.uniprot_id and args.transcript_id:
        # Load production data and filter to specified uniprot/transcript
        print(
            f"Loading production data for UniProt ID: {args.uniprot_id}, Transcript ID: {args.transcript_id}"
        )

        # Load production tables
        af2_dist_ht = get_af2_dist_ht().ht()
        obs_exp_ht = get_obs_exp_ht().ht()
        gencode_pos_ht = get_gencode_pos_ht().ht()
        pae_ht = get_af2_pae_ht().ht()
        plddt_ht = get_af2_plddt_ht().ht()

        # Filter to specified uniprot_id and transcript_id
        af2_ht = af2_dist_ht.filter(
            af2_dist_ht.uniprot_id == args.uniprot_id
        ).checkpoint(
            f"gs://gnomad-tmp-4day/proemis3d_test_data/af2_test.uniprot_id_{args.uniprot_id}.ht",
            _read_if_exists=True,
        )
        pae_ht = pae_ht.filter(pae_ht.uniprot_id == args.uniprot_id).checkpoint(
            f"gs://gnomad-tmp-4day/proemis3d_test_data/pae_test.uniprot_id_{args.uniprot_id}.ht",
            _read_if_exists=True,
        )
        plddt_ht = plddt_ht.filter(plddt_ht.uniprot_id == args.uniprot_id).checkpoint(
            f"gs://gnomad-tmp-4day/proemis3d_test_data/plddt_test.uniprot_id_{args.uniprot_id}.ht",
            _read_if_exists=True,
        )

        # Filter gencode_pos_ht to specified uniprot_id and transcript_id
        # This is the key filtering - generate_codon_oe_table will only process
        # positions in this table
        gencode_pos_ht = (
            gencode_pos_ht.filter(
                (gencode_pos_ht.uniprot_id == args.uniprot_id)
                & (gencode_pos_ht.enst == args.transcript_id)
            )
            .checkpoint(
                f"gs://gnomad-tmp-4day/proemis3d_test_data/gencode_pos_test.uniprot_id_{args.uniprot_id}.transcript_id_{args.transcript_id}.ht",
                _read_if_exists=True,
            )
            .naive_coalesce(1)
        )

        # Filter obs_exp_ht to specified transcript_id.
        obs_exp_ht = (
            obs_exp_ht.filter(obs_exp_ht.transcript == args.transcript_id)
            .checkpoint(
                f"gs://gnomad-tmp-4day/proemis3d_test_data/obs_exp_test.transcript_id_{args.transcript_id}.ht",
                _read_if_exists=True,
            )
            .naive_coalesce(1)
        )

        # obs_exp_ht is keyed by (locus, alleles) and needs to be aggregated by (locus, transcript)
        # before being passed to generate_codon_oe_table
        print("Aggregating obs_exp_ht by (locus, transcript)...")
        obs_exp_ht = (
            obs_exp_ht.group_by("locus", "transcript")
            .aggregate(
                obs=hl.agg.sum(obs_exp_ht.calibrate_mu.observed_variants[0]),
                exp=hl.agg.sum(obs_exp_ht.expected_variants[0]),
            )
            .checkpoint(
                f"gs://gnomad-tmp-4day/proemis3d_test_data/obs_exp_test.transcript_id_{args.transcript_id}.agg.ht",
                _read_if_exists=True,
            )
            .naive_coalesce(1)
        )

        # Generate oe_codon_ht from obs_exp_ht and gencode_pos_ht
        print("Generating oe_codon_ht from obs_exp_ht and gencode_pos_ht...")
        oe_codon_ht = generate_codon_oe_table(obs_exp_ht, gencode_pos_ht)

        # Filter oe_codon_ht to specified uniprot_id and transcript_id
        # oe_by_transcript is an array of structs with 'enst' field
        oe_codon_ht = oe_codon_ht.filter(oe_codon_ht.uniprot_id == args.uniprot_id)
        oe_codon_ht = oe_codon_ht.annotate(
            oe_by_transcript=oe_codon_ht.oe_by_transcript.filter(
                lambda x: x.enst == args.transcript_id
            )
        )
        oe_codon_ht = (
            oe_codon_ht.filter(hl.len(oe_codon_ht.oe_by_transcript) > 0)
            .checkpoint(
                f"gs://gnomad-tmp-4day/proemis3d_test_data/oe_codon_test.uniprot_id_{args.uniprot_id}.transcript_id_{args.transcript_id}.ht",
                _read_if_exists=True,
            )
            .naive_coalesce(1)
        )

        print("Loaded production data (filtered):")
        print(f"  af2_ht: {af2_ht.count()} rows")
        print(f"  oe_codon_ht: {oe_codon_ht.count()} rows")
        print(f"  pae_ht: {pae_ht.count()} rows")
        print(f"  plddt_ht: {plddt_ht.count()} rows")

        if af2_ht.count() == 0:
            print(f"Error: No data found for UniProt ID {args.uniprot_id}")
            sys.exit(1)
        if oe_codon_ht.count() == 0:
            print(
                f"Error: No OE data found for UniProt ID {args.uniprot_id} and Transcript ID {args.transcript_id}"
            )
            sys.exit(1)
    else:
        # Load test datasets
        # Use GCS path for cloud environments, or local path for local testing
        if args.local:
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

    # Handle list-tests flag
    if args.list_tests:
        list_tests()
        return

    # Determine which tests to run
    # When using production data (uniprot_id/transcript_id), default to
    # no_filtering test
    if args.uniprot_id and args.transcript_id:
        if args.test:
            # Run a specific test
            if args.test not in ALL_TESTS:
                print(f"Error: Unknown test '{args.test}'")
                print("Use --list-tests to see available tests")
                sys.exit(1)
            tests_to_run = [args.test]
        elif args.group:
            print(
                "Warning: --group is not supported when using --uniprot-id/--transcript-id"
            )
            print("Running 'no_filtering' test instead")
            tests_to_run = ["no_filtering"]
        else:
            # Default to no_filtering when using production data
            tests_to_run = ["no_filtering"]
    else:
        # Using test data - normal test selection logic
        if args.test:
            # Run a specific test
            if args.test not in ALL_TESTS:
                print(f"Error: Unknown test '{args.test}'")
                print("Use --list-tests to see available tests")
                sys.exit(1)
            tests_to_run = [args.test]
        elif args.group:
            # Run a group of tests
            if args.group not in TEST_GROUPS:
                print(f"Error: Unknown test group '{args.group}'")
                print("Use --list-tests to see available groups")
                sys.exit(1)
            tests_to_run = TEST_GROUPS[args.group]
        else:
            # Run all tests by default
            tests_to_run = TEST_GROUPS["all"]

    # Run the selected tests
    run_tests(
        tests_to_run,
        af2_ht,
        oe_codon_ht,
        pae_ht,
        plddt_ht,
        run_forward_tests=args.run_forward,
    )

    print("\nAll tests completed!")

    hl.stop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test PAE and pLDDT filtering parameters with test data"
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="Use local test data instead of GCS",
    )
    parser.add_argument(
        "--create-test-data",
        action="store_true",
        help="Create test datasets instead of running tests",
    )
    parser.add_argument(
        "--test",
        type=str,
        help="Run a specific test by name (e.g., 'no_filtering', 'pae_truncate')",
    )
    parser.add_argument(
        "--group",
        type=str,
        help="Run a group of tests (e.g., 'plddt_methods', 'pae_methods', 'combined')",
    )
    parser.add_argument(
        "--list-tests",
        action="store_true",
        help="List all available test names and groups",
    )
    parser.add_argument(
        "--run-forward",
        action="store_true",
        help="Also run run_forward tests after determine_regions_with_min_oe_upper",
    )
    parser.add_argument(
        "--uniprot-id",
        type=str,
        help="UniProt ID to use instead of test data (requires --transcript-id)",
    )
    parser.add_argument(
        "--transcript-id",
        type=str,
        help="Transcript ID (ENST) to use instead of test data (requires --uniprot-id)",
    )

    args = parser.parse_args()

    # Validate that both uniprot_id and transcript_id are provided together
    if (args.uniprot_id is None) != (args.transcript_id is None):
        parser.error("--uniprot-id and --transcript-id must be provided together")

    main(args)
