"""
Create small test datasets for testing PAE and pLDDT filtering parameters.

This script creates minimal test datasets that can be used to verify:
- PAE cutoff methods (truncate_on_pairwise_pae_with_center, filter_on_pairwise_pae_with_center, filter_on_pairwise_pae_in_region)
- pLDDT cutoff methods (truncate_at_first_low_plddt, remove_low_plddt_residues, exclude_low_plddt_from_stats)
"""

import hail as hl
import os

# Initialize Hail
hl.init(tmp_dir="gs://gnomad-tmp-4day")

# Test parameters
NUM_RESIDUES = 15  # Small protein with 15 residues
TEST_UNIPROT_ID = "P12345"
TEST_TRANSCRIPT_ID = "ENST00000123456"

# Create distance matrix data
# For each residue, create distances to all other residues
# Structure: dist_mat is an array of structs with {dist: float32, residue_index: int32}
# First, create a symmetric distance matrix
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
            raise ValueError(f"Self-distance for residue {center_idx} is {dist}, should be 0.0")
        # Verify non-self distances are > 0
        if i != center_idx and dist == 0.0:
            raise ValueError(f"Non-self distance from {center_idx} to {i} is 0.0, should be > 0")
        dist_mat.append({"dist": float(dist), "residue_index": int(i)})
    # Do NOT sort - keep in residue index order
    return dist_mat


# Create symmetric distance matrix first
symmetric_dist_matrix = create_symmetric_dist_matrix(NUM_RESIDUES)

# Create AF2 HT (AlphaFold2 data)
# Key: uniprot_id, aa_index
# Fields: enst, dist_mat
af2_dicts = []
for aa_idx in range(NUM_RESIDUES):
    dist_mat_list = create_dist_mat(NUM_RESIDUES, aa_idx, symmetric_dist_matrix)
    af2_dicts.append({
        "uniprot_id": TEST_UNIPROT_ID,
        "aa_index": aa_idx,
        "enst": TEST_TRANSCRIPT_ID,
        "dist_mat": dist_mat_list,
    })

# Use hl.Table.parallelize to create from Python data
af2_ht = hl.Table.parallelize([
    hl.struct(
        uniprot_id=d["uniprot_id"],
        aa_index=d["aa_index"],
        enst=d["enst"],
        dist_mat=hl.array([
            hl.struct(dist=hl.float32(x["dist"]), residue_index=hl.int32(x["residue_index"]))
            for x in d["dist_mat"]
        ])
    )
    for d in af2_dicts
])
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
    oe_array.append({
        "obs": int(obs),
        "exp": float(exp),
    })

# Structure: oe_by_transcript is array of {enst: str, oe: array<{obs: int64, exp: float64}>}
oe_codon_ht = hl.Table.parallelize([
    hl.struct(
        uniprot_id=TEST_UNIPROT_ID,
        oe_by_transcript=hl.array([
            hl.struct(
                enst=TEST_TRANSCRIPT_ID,
                oe=hl.array([
                    hl.struct(
                        obs=hl.int64(oe_item["obs"]),
                        exp=hl.float64(oe_item["exp"]),
                    )
                    for oe_item in oe_array
                ])
            )
        ])
    )
])
oe_codon_ht = oe_codon_ht.key_by("uniprot_id")

# Create PAE HT (Predicted Aligned Error data)
# Key: uniprot_id, aa_index
# Fields: pae (array of PAE values, one per residue)
# PAE[x][y] = expected position error at residue x if aligned on residue y
# PAE should be:
# - Low for residues close in sequence (same domain)
# - Higher for residues far in sequence (different domains)
# - Correlated with 3D distance (but not perfectly)
# - Generally symmetric-ish (PAE[x][y] ≈ PAE[y][x] in many cases)
# For testing: some residues will have high PAE (>15) to test filtering

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
                struct_influence = (struct_dist - 5.0) * 0.1  # Scale 3D distance contribution
                
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

# Create PAE matrix
pae_matrix = create_pae_matrix(NUM_RESIDUES, symmetric_dist_matrix)

# Convert to the format needed for the HT
pae_dicts = []
for aa_idx in range(NUM_RESIDUES):
    pae_array = []
    for other_idx in range(NUM_RESIDUES):
        pae_array.append(pae_matrix[aa_idx][other_idx])
    pae_dicts.append({
        "uniprot_id": TEST_UNIPROT_ID,
        "aa_index": aa_idx,
        "pae": pae_array,
    })

pae_ht = hl.Table.parallelize([
    hl.struct(
        uniprot_id=d["uniprot_id"],
        aa_index=d["aa_index"],
        pae=hl.array([hl.float32(x) for x in d["pae"]])
    )
    for d in pae_dicts
])
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
    plddt_dicts.append({
        "uniprot_id": TEST_UNIPROT_ID,
        "aa_index": aa_idx,
        "plddt": float(plddt),
    })

plddt_ht = hl.Table.parallelize([
    hl.struct(
        uniprot_id=d["uniprot_id"],
        aa_index=d["aa_index"],
        plddt=hl.float32(d["plddt"])
    )
    for d in plddt_dicts
])
plddt_ht = plddt_ht.key_by("uniprot_id", "aa_index")

# Write test datasets
# Use GCS path for cloud environments, or local path for local testing
import sys
if "--local" in sys.argv:
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
