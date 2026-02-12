# Proemis3D Filtering Options and Model Comparison Methods

This document provides a comprehensive guide to the filtering options and forward algorithm model comparison methods available in Proemis3D.

---

## Table of Contents

1. [Overview](#overview)
2. [AlphaFold2 Confidence Filtering](#alphafold2-confidence-filtering)
   - [pLDDT Filtering Methods](#plddt-filtering-methods)
   - [PAE Filtering Methods](#pae-filtering-methods)
   - [Filtering Order](#filtering-order)
3. [Forward Algorithm Model Comparison](#forward-algorithm-model-comparison)
   - [AIC Method](#aic-method)
   - [AIC Weight Method](#aic-weight-method)
   - [LRT Method](#lrt-method)
4. [Quick Reference Tables](#quick-reference-tables)
5. [Example: Complete Walkthrough](#example-complete-walkthrough)
6. [Command Line Usage](#command-line-usage)

---

## Overview

Proemis3D uses AlphaFold2 structural predictions to identify protein regions that are intolerant to missense variation. The pipeline applies two types of filtering based on AlphaFold2 confidence metrics:

1. **pLDDT (predicted Local Distance Difference Test)**: Per-residue confidence score (0-100)
2. **PAE (Predicted Aligned Error)**: Pairwise confidence in relative residue positions (in Ångströms)

After filtering, the forward algorithm uses model comparison to decide whether to add each candidate region.

---

## AlphaFold2 Confidence Filtering

### Understanding AlphaFold2 Confidence Scores

#### pLDDT Score
- **Range**: 0-100 (higher = more confident)
- **Interpretation**:
  - > 90: Very high confidence (typically for well-structured cores)
  - 70-90: High confidence (typical for most structured regions)
  - 50-70: Low confidence (often disordered regions or surface loops)
  - < 50: Very low confidence (likely disordered)
- **Typical cutoff**: 70

#### PAE Score
- **Range**: 0-31.75 Å (lower = more confident in relative position)
- **Interpretation**:
  - < 5 Å: Very high confidence in relative position (same domain)
  - 5-15 Å: Moderate confidence (possibly different domains)
  - > 15 Å: Low confidence (domains may move independently)
- **Typical cutoff**: 15 Å

---

### pLDDT Filtering Methods

There are **two** pLDDT filtering methods, each with different behavior:

#### 1. `remove_low_plddt_residues` (Hard Removal)

**Behavior**: Removes only residues with low pLDDT; continues evaluating others.

**Use case**: When low-confidence residues are scattered but nearby residues may still be valid.

```
Residues by distance from center:
Position:   [C] [1] [2] [3] [4] [5] [6]
pLDDT:       85  90  45  80  30  75  82
                    ✗       ✗
                 Removed  Removed

Result:     [C] [1] [3] [5] [6]
            Region = [C, 1, 3, 5, 6]
```

**Effect on null model**: Removed residues are NOT included in null model.

---

#### 2. `exclude_low_plddt_from_stats` (Soft Exclusion)

**Behavior**: Keeps all residues in the region for assignment purposes, but excludes low-confidence residues from statistical calculations (Obs, Exp, OE, NLL).

**Use case**: When you want regions to be spatially contiguous but not trust the statistics of uncertain residues.

```
Residues by distance from center:
Position:   [C] [1] [2] [3] [4] [5] [6]
pLDDT:       85  90  45  80  30  75  82
                   (ex)    (ex)
               Excluded from stats

Result:     Region = [C, 1, 2, 3, 4, 5, 6]  (all included)
            Stats calculated from: [C, 1, 3, 5, 6]
```

**Effect on null model**: All residues remain in null model, but low-pLDDT ones excluded from stats.

---

### PAE Filtering Methods

There are **five** PAE filtering methods. The first three remove residues; the last two keep residues but exclude them from stats (analogous to pLDDT `remove` vs `exclude`).

#### 1. `truncate_on_pairwise_pae_with_center` (Sequential Truncation)

**Behavior**: Stops at the first residue whose PAE to the center exceeds the cutoff.

**Use case**: When residues are sorted by distance, and you want a contiguous region with good PAE to center.

```
Center residue C, residues sorted by distance:
Position:   [C] [1] [2] [3] [4] [5]
PAE to C:    0   5  10  20   8   4
                        ↑
               First high PAE to center

Result:     [C] [1] [2] — TRUNCATED
            Region = [C, 1, 2]
```

**Effect on null model**: Since PAE is a matrix and truncation is center-specific, a residue excluded from candidate regions for one center may still be included for another center. The null model includes all residues that appear in **any** candidate region from **any** center. A residue is only excluded if it comes after the truncation point for **all** possible centers (which is rare in practice).

---

#### 2. `filter_on_pairwise_pae_with_center` (Individual Filtering)

**Behavior**: Removes only residues with high PAE to center; keeps others regardless of position.

**Use case**: When you want to include distant residues that still have good PAE to center.

```
Center residue C, residues sorted by distance:
Position:   [C] [1] [2] [3] [4] [5]
PAE to C:    0   5  10  20   8   4
                        ✗
                     Removed

Result:     [C] [1] [2] [4] [5]
            Region = [C, 1, 2, 4, 5]
```

**Effect on null model**: Since PAE is a matrix, each center has different PAE values to each residue. A residue excluded from candidate regions for one center (PAE(center, residue) > max) may still be included for another center. The null model includes all residues that appear in **any** candidate region from **any** center. A residue is only excluded if it has high PAE to **all** possible centers (rare in practice).

---

#### 3. `filter_on_pairwise_pae_in_region` (Region-Wide Filtering)

**Behavior**: For each candidate residue, check PAE to ALL residues already in the region. Only add if PAE is below cutoff for all.

**Use case**: Most stringent—ensures all residues in the region have good pairwise PAE to each other.

**Full PAE Matrix** (max_pae = 15, values > 15 marked with *):
```
      C    1    2    3    4    5
C     0    5   10   12    6    4
1     5    0    8   18*   9    7
2    10    8    0   14   13   11
3    12   18*  14    0   16*  17*
4     6    9   13   16*    0    8
5     4    7   11   17*    8    0
```

**Building region starting from center C**:

```
Step 1: Region = [C]
        Check: None (center is always included)

Step 2: Add residue 1?
        Check: PAE(1,C) = 5  ✓ (≤ 15)
        → Region = [C, 1]

Step 3: Add residue 2?
        Check: PAE(2,C) = 10  ✓ (≤ 15)
               PAE(2,1) = 8   ✓ (≤ 15)
        → Region = [C, 1, 2]

Step 4: Add residue 3?
        Check: PAE(3,C) = 12  ✓ (≤ 15)
               PAE(3,1) = 18* ✗ (> 15) → EXCLUDED
        (No need to check PAE(3,2) since PAE(3,1) already fails)

Step 5: Add residue 4?
        Check: PAE(4,C) = 6   ✓ (≤ 15)
               PAE(4,1) = 9   ✓ (≤ 15)
               PAE(4,2) = 13  ✓ (≤ 15)
        → Region = [C, 1, 2, 4]

Step 6: Add residue 5?
        Check: PAE(5,C) = 4   ✓ (≤ 15)
               PAE(5,1) = 7   ✓ (≤ 15)
               PAE(5,2) = 11  ✓ (≤ 15)
               PAE(5,4) = 8   ✓ (≤ 15)
        → Region = [C, 1, 2, 4, 5]

Final Result: Region = [C, 1, 2, 4, 5]
```

**Key insight**: Residue 3 is excluded because PAE(3,1) = 18 > 15, even though PAE(3,C) = 12 ≤ 15. This ensures all pairwise PAE values within the region are below the threshold.

**Effect on null model**: Since candidate regions are computed for multiple center residues, a residue excluded from candidate regions for one center may still be included for another center. A residue is only excluded if it cannot be added to any region from any center (i.e., it has high PAE to all other residues, making it impossible to include in any candidate region). In practice, most residues can be added to at least one region from at least one center, so the null model typically includes all residues that passed pLDDT filtering.

---

#### 4. `exclude_on_pairwise_pae_with_center` (Exclude from stats: center–residue PAE)

**Behavior**: Like `filter_on_pairwise_pae_with_center`, but residues with PAE(center, neighbor) > max_pae are **not removed**. They stay in the region and are marked **exclude_from_stats** (obs/exp set to NA, excluded from cumulative OE, NLL, etc.).

**Use case**: When you want to keep spatially contiguous assignment (all residues in region) but not use high-PAE-to-center residues in statistical calculations—analogous to pLDDT `exclude_low_plddt_from_stats`.

```
Center residue C, residues sorted by distance:
Position:   [C] [1] [2] [3] [4] [5]
PAE to C:    0   5  10  20   8   4
                        (ex)
                     Excluded from stats only

Result:     Region = [C, 1, 2, 3, 4, 5]  (all included)
            Stats calculated from: [C, 1, 2, 4, 5]
```

**Effect on null model**: All residues remain in the region (and in valid_residues). The null model includes all residues; only obs/exp for high-PAE residues are NA in stats.

---

#### 5. `exclude_on_pairwise_pae_in_region` (Exclude from stats: region-wide PAE)

**Behavior**: Like `filter_on_pairwise_pae_in_region`, but residues that would fail the pairwise-PAE-in-region check are **not removed**. They stay in the region and are marked **exclude_from_stats** (obs/exp set to NA, excluded from cumulative OE, NLL, etc.).

**Use case**: When you want to keep all residues in the region for assignment but exclude from stats those with high pairwise PAE to others in the region—analogous to pLDDT `exclude_low_plddt_from_stats` for PAE.

```
Same PAE matrix as filter_on_pairwise_pae_in_region:
Residue 3 would be excluded (PAE(3,1) > 15). With exclude_on_pairwise_pae_in_region:
Residue 3 stays in region but is marked exclude_from_stats.

Result:     Region = [C, 1, 2, 3, 4, 5]  (all included)
            Stats calculated from: [C, 1, 2, 4, 5]  (3 excluded from stats)
```

**Effect on null model**: All residues remain in the region (and in valid_residues). The null model includes all residues; only obs/exp for residues with high pairwise PAE in region are NA in stats.

---

### Filtering Order

**Important**: pLDDT filtering is always applied BEFORE PAE filtering.

```
Input residues
      ↓
┌─────────────────────┐
│   pLDDT Filtering   │  ← Applied first
└─────────────────────┘
      ↓
┌─────────────────────┐
│    PAE Filtering    │  ← Applied second
└─────────────────────┘
      ↓
Final candidate region
```

This order makes sense because:
1. pLDDT is a per-residue confidence—if a residue has low confidence, its PAE values are also unreliable
2. Filtering low-pLDDT residues first reduces the work for PAE filtering

---

### Null Region Construction

The **null region** represents the set of residues that would remain if a candidate region is selected. It is computed for each candidate during `determine_regions_with_min_oe_upper` and stored in the `null_region` field of each `min_oe_upper` entry.

**General rule**: Null region = All valid residues - Selected region residues

**Valid residues** are determined by the filtering methods applied. The null region includes all residues that:
1. Passed pLDDT filtering (if applied)
2. Passed PAE filtering (if applied)
3. Are not in the selected candidate region

#### pLDDT Filtering Effects on Null Region

| Method | Valid Residues Definition |
|--------|--------------------------|
| `remove_low_plddt_residues` | All residues with pLDDT ≥ min_plddt (excluding removed ones) |
| `exclude_low_plddt_from_stats` | All residues (low pLDDT included but excluded from stats) |
| No pLDDT filtering | All residues (no pLDDT-based exclusion) |

#### PAE Filtering Effects on Null Region

The null region is constructed from `valid_residues`, which is the **union of all residues present in any candidate region**. Therefore, residues that are excluded from ALL candidate regions by PAE filtering will also be excluded from the null model.

| Method | Effect on Null Region |
|--------|----------------------|
| `truncate_on_pairwise_pae_with_center` | **Residues are only excluded if they come after truncation point for ALL centers**<br><br>Since PAE is a matrix, each center has its own truncation point. A residue excluded from candidate regions for center A may still be included for center B. The null model includes all residues that appear in **any** candidate region from **any** center. A residue is only excluded if it comes after the truncation point for **all** possible centers (rare in practice).<br><br>**Example**: Residue 5 might be excluded when center 0 is used (PAE(0,5) > max), but included when center 3 is used (PAE(3,5) ≤ max). Residue 5 would still be in the null model. |
| `filter_on_pairwise_pae_with_center` | **Residues are only excluded if they have PAE(center, residue) > max_pae for ALL centers**<br><br>Since PAE is a matrix, each center has different PAE values to each residue. A residue excluded from candidate regions for center A (PAE(A, residue) > max) may still be included for center B (PAE(B, residue) ≤ max). The null model includes all residues that appear in **any** candidate region from **any** center. A residue is only excluded if it has high PAE to **all** possible centers (rare in practice).<br><br>**Example**: Residue 3 might be excluded when center 0 is used (PAE(0,3) = 20 > 15), but included when center 4 is used (PAE(4,3) = 10 ≤ 15). Residue 3 would still be in the null model. |
| `filter_on_pairwise_pae_in_region` | **Residues that cannot be added to any region from any center are excluded from null model**<br><br>Since candidate regions are computed for multiple center residues, a residue excluded from candidate regions for one center may still be included for another center. A residue is only excluded if it cannot be added to any region from any center (i.e., it has high PAE to all other residues, making it impossible to include in any candidate region). In practice, most residues can be added to at least one region from at least one center, so the null model typically includes all residues that passed pLDDT filtering.<br><br>**Example**: Residue 3 might be excluded from candidate regions for center 0 (PAE(3,0) > 15 or PAE(3,X) > 15 for some X in region), but included in candidate regions for center 4 (all pairwise PAE ≤ 15). Residue 3 would still be in the null model. |
| `exclude_on_pairwise_pae_with_center` | **All residues remain in region and in null.** Residues with PAE(center, residue) > max_pae are kept but excluded from stats only (obs/exp NA). Valid residues = all residues (no removal). Null region = valid_residues − selected region. |
| `exclude_on_pairwise_pae_in_region` | **All residues remain in region and in null.** Residues with high pairwise PAE in region are kept but excluded from stats only (obs/exp NA). Valid residues = all residues (no removal). Null region = valid_residues − selected region. |
| No PAE filtering | All residues that passed pLDDT filtering (if any) are included in null model. |

**Key insight**: The null model includes only residues that appear in at least one candidate region. Since candidate regions are computed for multiple center residues, and PAE is a matrix (center-specific), a residue excluded from candidate regions for one center may still be included for another center. Therefore, residues are only excluded from the null model if they are excluded from candidate regions for **all** possible centers.

**Example with `truncate_on_pairwise_pae_with_center` (max_pae = 15)**:
```
PAE matrix (showing PAE from each center to each residue):
         Residue:  0    1    2    3    4    5    6    7
PAE from 0:        0    5   10   20    8    4   12   18
PAE from 1:        5    0    8   18   10    6   14   20
PAE from 4:        8   10   13    6    0    5    9   11

Center 0, sorted by distance: [0, 1, 2, 4, 5, 3, 6, 7]
  PAE to 0:                     0   5  10   8   4  20  12  18
                                ↑
                    First PAE > 15 at residue 3
  Candidate region (center 0): [0, 1, 2]  (stops at residue 3)

Center 4, sorted by distance: [4, 5, 0, 1, 2, 6, 7, 3]
  PAE to 4:                     0   5   8  10  13   9  11   6
  Candidate region (center 4): [4, 5, 0, 1, 2, 6, 7, 3]  (all pass)

Valid residues (union of all candidate regions): [0, 1, 2, 3, 4, 5, 6, 7]
  (Residue 3 excluded from center 0's region, but included in center 4's region)

If selected region = [0, 1, 2] (from center 0):
Null region: [3, 4, 5, 6, 7]  (valid - selected)
```

**Example with `filter_on_pairwise_pae_with_center` (max_pae = 15)**:
```
Center residue 0:
Residue:      0    1    2    3    4    5    6    7
PAE to 0:      0    5   10   20    8    4   12   18
                              ✗
                    Filtered out (PAE > 15)

Candidate regions include: [0, 1, 2, 4, 5, 6, 7]  (residue 3 excluded)
Valid residues:            [0, 1, 2, 4, 5, 6, 7]
If selected region = [0, 1, 2]:
Null region:                [4, 5, 6, 7]  (valid - selected)
```

**Example with `filter_on_pairwise_pae_in_region` (max_pae = 15)**:
```
Using PAE matrix from earlier example:
- Residue 3 has PAE(3,1) = 18 > 15, so it's excluded from regions containing residue 1
- But residue 3 might be included in other regions (e.g., with different centers)
- Most residues can be added to at least one region

Valid residues: Typically all residues that passed pLDDT filtering
                 (unless a residue has high PAE to ALL other residues)
```

The null region is precomputed and stored in `min_oe_upper` to avoid recomputation during the forward algorithm.

---

## Forward Algorithm Model Comparison

The forward algorithm iteratively adds regions to explain missense constraint. At each round, it compares:
- **Current model**: The null + already selected regions
- **Candidate model**: Current model + best candidate region

### Model Evaluation Metrics

All methods use the **Negative Log-Likelihood (NLL)**:

```math
NLL = Σ [λᵢ - obsᵢ · log(λᵢ) + log(obsᵢ!)]
```

Where `λᵢ = γ̂ · expᵢ` is the Poisson mean for residue i.

---

### AIC Method

**Parameter**: `model_comparison_method="aic"`

**Decision rule**: Accept candidate if `AIC_candidate < AIC_current`

**Formula**:
```math
AIC = 2k + 2·NLL
```
Where k = number of parameters (1 per region + 1 for null)

**Characteristics**:
- Simple and interpretable
- Balances fit (NLL) with complexity (k)
- No explicit significance threshold
- Tends to be permissive (adds regions easily)

**Example**:
```
Current model: 2 regions, NLL = 100
Candidate model: 3 regions, NLL = 95

AIC_current   = 2(3) + 2(100) = 206
AIC_candidate = 2(4) + 2(95)  = 198

198 < 206 → Accept candidate
```

---

### AIC Weight Method

**Parameters**:
- `model_comparison_method="aic_weight"`
- `aic_weight_thresh=0.80` (default)

**Decision rule**: Accept candidate if Akaike weight ≥ threshold

**Formula**:
```math
w_{candidate} = 1 / (1 + exp(0.5 · (AIC_{candidate} - AIC_{current})))
```

**Characteristics**:
- Transforms AIC difference to probability scale (0-1)
- Threshold controls stringency
- `w = 0.5` means models are equally likely
- Higher threshold = more conservative

**Weight interpretation**:
| Δ AIC | Weight |
|-------|--------|
| -10   | 0.993  |
| -5    | 0.924  |
| -2    | 0.731  |
| 0     | 0.500  |
| +2    | 0.269  |
| +5    | 0.076  |

**Example**:
```
AIC_current = 206, AIC_candidate = 198
Δ AIC = 198 - 206 = -8

w = 1 / (1 + exp(0.5 × (-8))) = 1 / (1 + 0.018) = 0.982

0.982 ≥ 0.80 → Accept candidate
```

**Choosing the threshold**
The threshold is the minimum Akaike weight required to add a region. Higher values require stronger evidence (larger AIC improvement) and produce fewer, more conservative regions.

| Threshold | Effect | Roughly requires |
|-----------|--------|-------------------|
| **0.5** | Permissive | Candidate at least as good as current (ΔAIC ≤ 0) |
| **0.8** (default) | Balanced | Candidate clearly better (ΔAIC ≲ −2.5 to −3) |
| **0.9** | Conservative | Candidate much better (ΔAIC ≲ −4 to −5) |

So: **0.5** = add when the candidate is not worse; **0.8** = add only when it is clearly better; **0.9** = add only when it is much better.

**How AIC weight differs from plain AIC**
Plain AIC uses a binary rule: accept if `AIC_candidate < AIC_current`. AIC weight uses the same AIC values but converts the difference to a number on (0, 1), then you choose a cutoff. So:

- **AIC**: One fixed rule (strict improvement only).
- **AIC weight**: Same information (AIC difference), but you control stringency via the threshold. You can make it as permissive as AIC (see below) or stricter (0.8, 0.9) so that only larger improvements add a region.

**Equivalence to plain AIC**
A threshold of **0.5** is equivalent to plain AIC in outcome: the weight is ≥ 0.5 exactly when ΔAIC ≤ 0, i.e. when the candidate model has AIC less than or equal to the current model. So `aic_weight` with `aic_weight_thresh=0.5` accepts whenever the candidate is at least as good as the current model by AIC. The only nuance is tie-breaking: at ΔAIC = 0, plain AIC rejects (strict `<`), while `aic_weight` with 0.5 accepts (≥ 0.5). So with threshold 0.5 you may add one extra region in the rare case of an exact AIC tie; otherwise the results match AIC. If you want identical behavior including ties, use `model_comparison_method="aic"`; if you want the same permissiveness with a tunable threshold, use `aic_weight` with `aic_weight_thresh=0.5`.

---

### LRT Method

**Parameters**:
- `model_comparison_method="lrt"`
- `lrt_alpha=0.001` (default)
- `lrt_df_added=1` (default)
- `bonferroni_per_round=True` (default)

**Decision rule**: Accept candidate if p-value ≤ α (possibly Bonferroni-adjusted)

**Formula**:
```math
D = 2 · (NLL_{current} - NLL_{candidate})
p = P(χ²_df ≥ D)
```

**Characteristics**:
- Most conservative method
- Formal hypothesis testing framework
- Bonferroni correction controls family-wise error rate
- α = 0.001 is stringent; typical values: 0.05, 0.01, 0.001

**Bonferroni adjustment**:
```math
α_{adjusted} = α / m_{candidates}
```
Where m_candidates is the number of candidate regions evaluated this round.

**Example**:
```
NLL_current = 100, NLL_candidate = 95
m_candidates = 50

D = 2 × (100 - 95) = 10
p = P(χ²₁ ≥ 10) ≈ 0.0016

α_adjusted = 0.001 / 50 = 0.00002

0.0016 > 0.00002 → Reject candidate (not significant after Bonferroni)
```

---

## Quick Reference Tables

### pLDDT Methods Summary

| Method | Low pLDDT residues | Effect on region | Effect on null |
|--------|-------------------|------------------|----------------|
| `remove_low_plddt_residues` | Remove only those | Gaps in region | Removed residues excluded |
| `exclude_low_plddt_from_stats` | Keep but exclude from stats | Spatially contiguous | In null but excluded from stats |

### PAE Methods Summary

| Method | What it checks | Effect on region | Effect on null |
|--------|---------------|------------------|----------------|
| `truncate_on_pairwise_pae_with_center` | PAE(residue, center) | Sequential region (stops at first high PAE) | Residues only excluded if after truncation for ALL centers |
| `filter_on_pairwise_pae_with_center` | PAE(residue, center) | Non-contiguous OK (removes high PAE residues) | High PAE to center residues excluded (if excluded from ALL centers) |
| `filter_on_pairwise_pae_in_region` | PAE(residue, all in region) | Most stringent (all pairwise PAE ≤ max) | Only residues with high PAE to ALL others excluded |
| `exclude_on_pairwise_pae_with_center` | PAE(residue, center) | Keep all; high PAE excluded from stats only | All residues in null (high PAE excluded from stats only) |
| `exclude_on_pairwise_pae_in_region` | PAE(residue, all in region) | Keep all; high pairwise PAE excluded from stats only | All residues in null (high PAE excluded from stats only) |

### Model Comparison Summary

| Method | Decision criterion | Stringency | Use when |
|--------|-------------------|------------|----------|
| `aic` | AIC_cand < AIC_curr | Permissive | Exploratory analysis |
| `aic_weight` | w ≥ threshold (default 0.8) | Tunable: 0.5 ≈ AIC, 0.8 balanced, 0.9 conservative | Balanced approach; same info as AIC with tunable cutoff |
| `lrt` | p ≤ α (Bonferroni) | Conservative | Statistical rigor required |

---

## Example: Complete Walkthrough

**Note**: This example is simplified to show the process for a single center. In practice, candidate regions are computed for **all** centers (residues 1-10), and the null model includes the **union of all residues** that appear in **any** candidate region from **any** center.

### Setup: Small Protein (10 residues)

```
Residue:      1    2    3    4    5    6    7    8    9   10
Observed:     0    0    0    0    1    1    1    1    1    1
Expected:   0.9  1.0  1.1  0.8  1.2  1.0  0.9  1.1  1.0  0.8
pLDDT:       85   90   45   80   70   95   60   85   90   75
```

### Step 1: Apply pLDDT Filtering (min_plddt = 70)

**Method: `remove_low_plddt_residues`**
```
[5] [4] [6] [3] [7] [2] [8] [1] [9] [10]
 70  80  95  45  60  90  85  85  90  75
              ✗   ✗
         Remove these

Candidate region for center 5: [5, 4, 6, 2, 8, 1, 9, 10]
```

**Method: `exclude_low_plddt_from_stats`**
```
Candidate region for center 5: [5, 4, 6, 3, 7, 2, 8, 1, 9, 10]
Stats calculated from:         [5, 4, 6,    2, 8, 1, 9, 10]  (3, 7 excluded)
```

### Step 2: Apply PAE Filtering (max_pae = 15)

**Note**: This example shows PAE filtering for center 5 only. In practice, candidate regions are computed for all centers (1-10). A residue excluded from one center's candidate region may still be included in another center's candidate region. For example, residue 3 is excluded from center 5's region (PAE(3,7)=16 > 15), but could be included in center 1's region (PAE(3,1)=8 ≤ 15, PAE(3,2)=6 ≤ 15, etc.).

Given PAE matrix (symmetric, showing values > 15 with *):
```
     1    2    3    4    5    6    7    8    9   10
1    0    5    8   12  18*  20*  22*  25*  20*  18*
2    5    0    6   10  16*  18*  20*  22*  18*  16*
3    8    6    0    4   12   14   16*  18*  15   13
4   12   10    4    0    6    8   12   14   10    8
5   18*  16*  12    6    0    3    6    8    5    4
6   20*  18*  14    8    3    0    4    6    4    5
7   22*  20*  16*  12    6    4    0    3    6    8
8   25*  22*  18*  14    8    6    3    0    5    7
9   20*  18*  15   10    5    4    6    5    0    3
10  18*  16*  13    8    4    5    8    7    3    0
```

**After `filter_on_pairwise_pae_in_region` from center 5:**
```
Starting region: [5]
Add 6? PAE(6,5)=3  ✓  → [5, 6]
Add 4? PAE(4,5)=6  ✓
       PAE(4,6)=8  ✓  → [5, 6, 4]
Add 7? PAE(7,5)=6  ✓
       PAE(7,6)=4  ✓
       PAE(7,4)=12 ✓  → [5, 6, 4, 7]
Add 10? PAE(10,5)=4 ✓
        PAE(10,6)=5 ✓
        PAE(10,4)=8 ✓
        PAE(10,7)=8 ✓ → [5, 6, 4, 7, 10]
Add 9? PAE(9,5)=5  ✓
       PAE(9,6)=4  ✓
       PAE(9,4)=10 ✓
       PAE(9,7)=6  ✓
       PAE(9,10)=3 ✓ → [5, 6, 4, 7, 10, 9]
Add 8? PAE(8,5)=8  ✓
       PAE(8,6)=6  ✓
       PAE(8,4)=14 ✓
       PAE(8,7)=3  ✓
       PAE(8,10)=7 ✓
       PAE(8,9)=5  ✓ → [5, 6, 4, 7, 10, 9, 8]
Add 3? PAE(3,5)=12 ✓
       PAE(3,6)=14 ✓
       PAE(3,4)=4  ✓
       PAE(3,7)=16* ✗ → EXCLUDED

Final region: [4, 5, 6, 7, 8, 9, 10]
```

### Step 3: Forward Algorithm Model Comparison

**Note**: In this simplified example, we only show candidate regions from center 5. In practice, candidate regions are computed for all centers, and the null model includes the union of all residues from all candidate regions. For example, residues 1, 2, or 3 might be included in candidate regions from other centers (e.g., center 1, 2, or 3), even though they're excluded from center 5's region.

**Null model** (all valid residues from center 5's candidate regions in this example):
```
Region: [4, 5, 6, 7, 8, 9, 10]
Obs = 0+1+1+1+1+1+1 = 6
Exp = 0.8+1.2+1.0+0.9+1.1+1.0+0.8 = 6.8
γ̂ = 6/6.8 = 0.882
NLL_null ≈ 8.12
```

**Best candidate region** [4, 5, 6] (center at 5):
```
Region obs/exp: Obs = 2, Exp = 3.0, γ̂_R = 0.667
Remaining:      Obs = 4, Exp = 3.8, γ̂_Rc = 1.053
NLL_candidate ≈ 7.89
```

**AIC comparison**:
```
AIC_null = 2(1) + 2(8.12) = 18.24
AIC_cand = 2(2) + 2(7.89) = 19.78

19.78 > 18.24 → Reject (AIC increased)
```

**AIC Weight**:
```
w = 1 / (1 + exp(0.5 × (19.78 - 18.24)))
  = 1 / (1 + exp(0.77))
  = 0.316

0.316 < 0.80 → Reject
```

**LRT**:
```
D = 2 × (8.12 - 7.89) = 0.46
p = P(χ²₁ ≥ 0.46) ≈ 0.50

0.50 > 0.001 → Reject
```

All three methods reject the candidate region—the null model is sufficient.

---

## Command Line Usage

### Basic Usage

```bash
# Run forward algorithm with AIC
python proemis_3d.py --run-forward --model-comparison-method aic

# Run with LRT and strict alpha
python proemis_3d.py --run-forward \
    --model-comparison-method lrt \
    --lrt-alpha 0.001 \
    --bonferroni-per-round

# Run with AIC weight
python proemis_3d.py --run-forward \
    --model-comparison-method aic_weight \
    --aic-weight-thresh 0.80
```

### Filtering Options

```bash
# pLDDT filtering only
python proemis_3d.py --run-forward \
    --plddt-cutoff 70 \
    --plddt-cutoff-method remove_low_plddt_residues

# PAE filtering only
python proemis_3d.py --run-forward \
    --pae-cutoff 15 \
    --pae-cutoff-method filter_on_pairwise_pae_in_region

# Combined filtering (recommended)
python proemis_3d.py --run-forward \
    --plddt-cutoff 70 \
    --plddt-cutoff-method exclude_low_plddt_from_stats \
    --pae-cutoff 15 \
    --pae-cutoff-method filter_on_pairwise_pae_in_region
```

### Testing with Debug Output

```bash
# Run test with debug output to file
python test_with_test_data.py \
    --test plddt_truncate \
    --run-forward-tests \
    --output-file debug_output.txt
```

---

## Recommended Configurations

### Conservative (Statistical Rigor)
```bash
--model-comparison-method lrt
--lrt-alpha 0.001
--bonferroni-per-round
--plddt-cutoff 70
--plddt-cutoff-method remove_low_plddt_residues
--pae-cutoff 15
--pae-cutoff-method filter_on_pairwise_pae_in_region
```

### Balanced (Default)
```bash
--model-comparison-method aic
--plddt-cutoff 70
--plddt-cutoff-method exclude_low_plddt_from_stats
--pae-cutoff 15
--pae-cutoff-method truncate_on_pairwise_pae_with_center
```

### Exploratory (Permissive)
```bash
--model-comparison-method aic_weight
--aic-weight-thresh 0.60
--pae-cutoff 20
--pae-cutoff-method filter_on_pairwise_pae_with_center
```

### PAE exclude-from-stats (keep residues, exclude from stats only)
```bash
# Center–residue PAE: high-PAE residues kept but excluded from stats
--pae-cutoff 15
--pae-cutoff-method exclude_on_pairwise_pae_with_center
```
```bash
# Region-wide PAE: high pairwise PAE residues kept but excluded from stats
--pae-cutoff 15
--pae-cutoff-method exclude_on_pairwise_pae_in_region
```

---

## See Also

- [README.md](README.md) - Main Proemis3D documentation
- [README_FORWARD_ALGORITHM.md](README_FORWARD_ALGORITHM.md) - Mathematical details of forward algorithm
- [README_GAMMA_UDF.md](README_GAMMA_UDF.md) - Gamma distribution for OE upper bounds
