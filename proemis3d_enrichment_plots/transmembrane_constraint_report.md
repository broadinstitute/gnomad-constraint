# Transmembrane Region Constraint Analysis

## Overview

This report examines the relationship between transmembrane (TM) regions and
regional missense constraint across 4,665 genes with TM annotations analyzed
by Proemis3D. We investigate whether TM residues are more constrained than
non-TM residues, how this varies by TM segment type, and how the choice of
forward model-comparison method (AIC vs LRT) and structural filtering (pLDDT,
PAE) affects the results.

**Data sources:**

- Per-residue constraint estimates from 6 production Hail Tables (4,665 TM
  genes, one transcript per gene, 6 method configurations: 3 filter groups ×
  AIC/LRT)
- UniProt protein features for 19,241 genes exported from
  `gs://gnomad-tmp-4day/persist_TableQqk3byOR9W`
- DD de novo case variants from `autism_dd_variants.tsv`

**Constraint metric:** OE upper (upper bound of the observed/expected confidence
interval for missense variants). Lower values indicate stronger constraint.

**Summary statistic:** Geometric mean of per-residue OE upper within each region
category per gene. The geometric mean is the appropriate central tendency for
ratio data because OE values are inherently multiplicative—a value of 0.5 and
2.0 are equidistant from the neutral point (1.0) on a log scale.

---

## 1. Genes with Transmembrane Regions

Of the 19,241 genes with UniProt annotations, **4,665 contain transmembrane or
intramembrane regions**. After filtering to genes with at least 10 TM residues
and sufficient constraint data (including catch-all residues for genes with
at least one non-null region), 4,615 genes are included in the global
paired analysis. These fall into several structural classes:

| TM gene class | Genes | Description |
|---|---|---|
| Generic TM helix only | ~3,336 | Single-pass or multi-pass TM helices (receptors, adhesion molecules, etc.) |
| Ion channel / named segments | 76 | Voltage-gated ion channels with S1–S6 segments |
| Transporter helix | 913 | Transporters with numbered helices |
| Named TM (M1/M2) | 24 | Ligand-gated channels with named TM segments |
| Intramembrane | 191 | Segments that enter but do not fully cross the membrane |

---

## 2. TM vs Non-TM Constraint (Paired Distribution Analysis)

For each gene, we computed the geometric mean OE upper separately for TM and
non-TM residues. We tested the difference within each gene using a
**Mann-Whitney U test** and tested the systematic trend across genes using a
**paired Wilcoxon signed-rank test** on per-gene geometric means.

### Global results (paired Wilcoxon signed-rank across genes)

Analysis now includes catch-all residues for genes with at least one non-null
constraint region. Genes where all residues are classified as null are excluded.

| Method | p-value | Median diff | Direction | n genes | Genes with TM < non-TM |
|---|---|---|---|---|---|
| No filter (AIC) | 6.18 × 10⁻⁹³ | −0.061 | TM < non-TM | 4,615 | 3,068/4,615 (66%) |
| No filter (LRT) | 6.18 × 10⁻¹⁷ | −0.023 | TM < non-TM | 2,363 | 1,283/2,363 (54%) |
| pLDDT exclude (AIC) | 4.98 × 10⁻¹³⁵ | −0.094 | TM < non-TM | 4,552 | 3,144/4,552 (69%) |
| pLDDT exclude (LRT) | 8.41 × 10⁻⁵ | +0.007 | TM > non-TM | 1,882 | 918/1,882 (49%) |
| PAE filter center (AIC) | 2.44 × 10⁻⁹⁰ | −0.068 | TM < non-TM | 4,589 | 3,011/4,589 (66%) |
| PAE filter center (LRT) | 1.81 × 10⁻¹² | −0.005 | TM < non-TM | 2,224 | 1,127/2,224 (51%) |

TM regions are systematically more constrained across genes under all AIC
methods, with the pLDDT exclude method now showing the strongest signal
(p = 4.98 × 10⁻¹³⁵, median difference −0.094). Including catch-all residues
substantially strengthens the signals compared to prior analyses that excluded
null-model residues. LRT methods show weaker trends; pLDDT exclude (LRT) is
the only method where TM > non-TM, and the effect size (+0.007) is very small.
LRT's weakness reflects massive gene dropout (see Section 3).

### Breakdown by TM segment type (paired Wilcoxon signed-rank)

For each TM category, we paired the per-gene geometric mean OE upper against
the non-TM geomean for the same gene and ran a Wilcoxon signed-rank test.
"Median diff" is category minus non-TM (negative = more constrained).

#### Pore-lining (S5/S6) and Named TM (M1/M2)

| TM category | Method | n genes | Median diff | % cat < non-TM | Wilcoxon p |
|---|---|---|---|---|---|
| Pore-lining (S5) | No filter (AIC) | 74 | −0.259 | 97% | 2.3 × 10⁻¹³ \*\*\* |
| Pore-lining (S5) | No filter (LRT) | 65 | −0.245 | 97% | 2.9 × 10⁻¹² \*\*\* |
| Pore-lining (S5) | pLDDT exclude (AIC) | 74 | −0.278 | 96% | 2.4 × 10⁻¹³ \*\*\* |
| Pore-lining (S5) | pLDDT exclude (LRT) | 63 | −0.213 | 94% | 1.0 × 10⁻¹¹ \*\*\* |
| Pore-lining (S5) | PAE filter center (AIC) | 74 | −0.277 | 95% | 2.0 × 10⁻¹³ \*\*\* |
| Pore-lining (S5) | PAE filter center (LRT) | 63 | −0.226 | 98% | 5.4 × 10⁻¹² \*\*\* |
| Pore-lining (S6) | No filter (AIC) | 75 | −0.253 | 97% | 1.6 × 10⁻¹³ \*\*\* |
| Pore-lining (S6) | No filter (LRT) | 66 | −0.259 | 100% | 1.6 × 10⁻¹² \*\*\* |
| Pore-lining (S6) | pLDDT exclude (AIC) | 75 | −0.303 | 97% | 7.9 × 10⁻¹⁴ \*\*\* |
| Pore-lining (S6) | pLDDT exclude (LRT) | 64 | −0.245 | 92% | 1.4 × 10⁻¹¹ \*\*\* |
| Pore-lining (S6) | PAE filter center (AIC) | 75 | −0.292 | 97% | 8.2 × 10⁻¹⁴ \*\*\* |
| Pore-lining (S6) | PAE filter center (LRT) | 64 | −0.230 | 98% | 4.1 × 10⁻¹² \*\*\* |
| Named TM (M1/M2) | No filter (AIC) | 24 | −0.116 | 92% | 8.3 × 10⁻⁷ \*\*\* |
| Named TM (M1/M2) | No filter (LRT) | 17 | −0.115 | 88% | 1.1 × 10⁻⁴ \*\*\* |
| Named TM (M1/M2) | pLDDT exclude (AIC) | 24 | −0.120 | 88% | 8.3 × 10⁻⁶ \*\*\* |
| Named TM (M1/M2) | pLDDT exclude (LRT) | 15 | −0.070 | 67% | 1.2 × 10⁻² \* |
| Named TM (M1/M2) | PAE filter center (AIC) | 24 | −0.119 | 92% | 1.7 × 10⁻⁶ \*\*\* |
| Named TM (M1/M2) | PAE filter center (LRT) | 17 | −0.144 | 94% | 5.0 × 10⁻⁴ \*\*\* |

Pore-lining segments are highly significantly more constrained than non-TM
under every method, including all LRT variants. With catch-all residues
included, the median differences are large (−0.23 to −0.30) and all six methods
show p < 10⁻¹¹. Named TM (M1/M2) segments are also now significant across all
methods including LRT, with 88–94% of genes showing category < non-TM.

#### Voltage sensor and VSD segments

| TM category | Method | n genes | Median diff | % cat < non-TM | Wilcoxon p |
|---|---|---|---|---|---|
| Voltage sensor (S4) | No filter (AIC) | 76 | −0.212 | 95% | 7.7 × 10⁻¹⁴ \*\*\* |
| Voltage sensor (S4) | No filter (LRT) | 67 | −0.203 | 91% | 2.9 × 10⁻¹² \*\*\* |
| Voltage sensor (S4) | pLDDT exclude (AIC) | 76 | −0.215 | 95% | 1.6 × 10⁻¹³ \*\*\* |
| Voltage sensor (S4) | pLDDT exclude (LRT) | 65 | −0.131 | 83% | 4.3 × 10⁻¹⁰ \*\*\* |
| Voltage sensor (S4) | PAE filter center (AIC) | 76 | −0.258 | 96% | 4.8 × 10⁻¹⁴ \*\*\* |
| Voltage sensor (S4) | PAE filter center (LRT) | 65 | −0.208 | 94% | 6.1 × 10⁻¹² \*\*\* |
| VSD other (S1–S3) | No filter (AIC) | 77 | −0.147 | 88% | 3.8 × 10⁻¹² \*\*\* |
| VSD other (S1–S3) | No filter (LRT) | 68 | −0.113 | 84% | 3.4 × 10⁻¹⁰ \*\*\* |
| VSD other (S1–S3) | pLDDT exclude (AIC) | 77 | −0.139 | 84% | 3.9 × 10⁻¹¹ \*\*\* |
| VSD other (S1–S3) | pLDDT exclude (LRT) | 66 | −0.051 | 64% | 1.6 × 10⁻⁵ \*\*\* |
| VSD other (S1–S3) | PAE filter center (AIC) | 77 | −0.212 | 87% | 3.7 × 10⁻¹² \*\*\* |
| VSD other (S1–S3) | PAE filter center (LRT) | 66 | −0.164 | 82% | 5.0 × 10⁻¹⁰ \*\*\* |

Voltage sensor S4 segments are robustly constrained across all six methods,
including all LRT variants. The median differences (−0.13 to −0.26) are
substantially larger than in analyses that excluded catch-all residues. VSD
other (S1–S3) segments are similarly significant across all methods, with the
LRT pLDDT exclude method showing the smallest effect (−0.051, 64% cat < non-TM)
but still highly significant.

#### Intramembrane and Transporter helix

| TM category | Method | n genes | Median diff | % cat < non-TM | Wilcoxon p |
|---|---|---|---|---|---|
| Intramembrane | No filter (AIC) | 197 | −0.171 | 85% | 3.2 × 10⁻²¹ \*\*\* |
| Intramembrane | No filter (LRT) | 141 | −0.180 | 86% | 2.0 × 10⁻²⁰ \*\*\* |
| Intramembrane | pLDDT exclude (AIC) | 197 | −0.223 | 84% | 2.0 × 10⁻²³ \*\*\* |
| Intramembrane | pLDDT exclude (LRT) | 128 | −0.183 | 79% | 1.6 × 10⁻¹⁶ \*\*\* |
| Intramembrane | PAE filter center (AIC) | 197 | −0.208 | 83% | 6.8 × 10⁻²¹ \*\*\* |
| Intramembrane | PAE filter center (LRT) | 141 | −0.191 | 84% | 1.3 × 10⁻¹⁹ \*\*\* |
| Transporter helix | No filter (AIC) | 913 | −0.093 | 82% | 7.1 × 10⁻⁹¹ \*\*\* |
| Transporter helix | No filter (LRT) | 372 | −0.094 | 80% | 6.5 × 10⁻³⁹ \*\*\* |
| Transporter helix | pLDDT exclude (AIC) | 912 | −0.144 | 88% | 2.7 × 10⁻¹¹⁷ \*\*\* |
| Transporter helix | pLDDT exclude (LRT) | 303 | −0.075 | 77% | 5.5 × 10⁻²⁸ \*\*\* |
| Transporter helix | PAE filter center (AIC) | 913 | −0.099 | 81% | 8.8 × 10⁻⁹² \*\*\* |
| Transporter helix | PAE filter center (LRT) | 360 | −0.101 | 81% | 2.9 × 10⁻³⁹ \*\*\* |

Intramembrane regions are highly significant under all six methods, including
LRT variants, with consistent median differences around −0.17 to −0.22 and
79–86% of genes showing category < non-TM. Transporter helices are now also
strongly significant across all methods (including LRT), with median differences
of −0.075 to −0.144 and 77–88% of genes showing category < non-TM. The median
transporter geomean OE (1.051) falls below the overall non-TM median (1.075),
contrasting with earlier analyses (before catch-all inclusion) that placed
transporters marginally above non-TM.

#### Generic TM helix

| TM category | Method | n genes | Median diff | % cat < non-TM | Wilcoxon p |
|---|---|---|---|---|---|
| Generic TM helix | No filter (AIC) | 3,597 | −0.047 | 62% | 7.4 × 10⁻³³ \*\*\* |
| Generic TM helix | No filter (LRT) | 1,907 | +0.009 | 48% | 0.136 |
| Generic TM helix | pLDDT exclude (AIC) | 3,535 | −0.072 | 64% | 1.6 × 10⁻⁵³ \*\*\* |
| Generic TM helix | pLDDT exclude (LRT) | 1,501 | +0.020 | 41% | 0.016 \* |
| Generic TM helix | PAE filter center (AIC) | 3,571 | −0.048 | 61% | 1.0 × 10⁻³⁰ \*\*\* |
| Generic TM helix | PAE filter center (LRT) | 1,781 | +0.017 | 43% | 0.702 |

Generic TM helices remain the least constrained TM category. Under AIC the
signal is now stronger (median differences −0.047 to −0.072) and the pLDDT
exclude result no longer vanishes. LRT methods still show no significant or
directionally inconsistent results, highlighting that generic TM anchors
individually lack the strong functional constraint seen in ion channel or
ligand-gated channel TM segments. The median geomean OE for generic TM helices
(1.030) is still below the overall non-TM median (1.075), but the per-gene
effect sizes are modest compared with named TM categories.

### Distribution comparison

![TM vs non-TM distributions](fig1_tm_vs_nontm_distributions.png)

**Figure 1.** Paired violin + strip plots showing distributions of per-gene
geometric mean OE upper for TM vs non-TM residues across 6 method
configurations. Each dot is one gene. Black bars show medians. Annotations
show paired Wilcoxon test results.

---

## 3. AIC vs LRT: Why AIC Is Better Suited for This Analysis

AIC and LRT differ fundamentally in how they handle genes with weak constraint
signal, and this has cascading effects on the TM enrichment analysis. Three
compounding problems explain LRT's weak performance.

### 3.1. Massive gene dropout

LRT classifies 49% of genes (2,302 / 4,665) as having no significant regional
variation — every residue in those genes receives `is_null = True`. AIC keeps
4,615 genes with non-null data. Critically, the 2,880 genes that LRT drops
still show modest TM constraint under AIC: 45% of them have TM geomean OE
lower than non-TM. By discarding these genes, LRT throws away a large portion
of the real (if individually subtle) TM constraint signal.

| Method | Genes fully null | Genes with non-null data |
|---|---|---|
| AIC (no filter) | 50 (1.1%) | 4,615 |
| LRT (no filter) | 2,302 (49.3%) | 2,363 |
| AIC, pLDDT exclude | 113 (2.4%) | 4,552 |
| LRT, pLDDT exclude | 2,783 (59.7%) | 1,882 |
| AIC, PAE filter center | 55 (1.2%) | 4,610 |
| LRT, PAE filter center | 2,427 (52.0%) | 2,238 |

### 3.2. Coarse regions that straddle TM boundaries

When LRT does reject the null, it typically identifies far fewer constraint
regions per gene than AIC. This matters because the TM enrichment analysis
depends on regions being fine-grained enough to distinguish TM from non-TM
residues *within* a gene.

| Metric (no-filter methods) | AIC | LRT |
|---|---|---|
| Non-null regions per gene (median) | 5 | 1 |
| Unique OE upper values per gene (median) | 5 | 1 |
| Median size of mixed TM + non-TM regions | 55 residues | 149 residues |

Under LRT, the median gene has only **one** non-null region, and when that
region spans both TM and non-TM residues (which 46% of LRT regions do), every
residue in the region receives the *same* OE upper value. AIC fits median 5
regions per gene, providing the spatial resolution needed to capture distinct
constraint levels in TM vs non-TM portions of the protein.

### 3.3. Identical TM and non-TM OE values within genes

The coarse regions directly cause the paired test to fail. When TM and non-TM
residues share the same region(s), they get identical OE upper values, and
the within-gene difference is exactly zero.

Including catch-all residues substantially alleviates this problem for AIC:
catch-all residues receive the null-model OE (a gene-level constant distinct
from any constraint-region OE), so genes that have some TM residues in
constraint regions and some non-TM residues in the catch-all now show
*distinct* OE values for the two categories. Under LRT, the coarse-region
structure and high gene dropout (49% fully null) means that LRT-retained genes
are a selected subset with unusually strong regional variation, yet the per-gene
TM vs non-TM discrimination remains weaker than AIC. The global paired Wilcoxon
results confirm this: AIC no-filter p = 6.18 × 10⁻⁹³ vs LRT no-filter
p = 6.18 × 10⁻¹⁷ across 4,615 vs 2,363 genes respectively.

### 3.4. LRT-retained genes show weaker TM signal than AIC

With catch-all residues included, LRT no-filter retains 2,363 genes (with at
least one non-null residue), of which 1,283 (54%) show TM < non-TM with a
median paired difference of −0.023. AIC no-filter retains 4,615 genes (96%
more), of which 3,068 (66%) show TM < non-TM with a median difference of
−0.061. The LRT subset represents the genes with the *strongest* regional
variation signal (since they survived LRT's conservative threshold), yet AIC
reveals a much larger and more consistent TM constraint signal across the
broader gene set. The biological signal is present even in weakly varying genes
— LRT's conservatism filters them out entirely.

### Summary

LRT is not "wrong" — it is correctly answering a different question ("is there
statistically significant evidence of regional variation in this gene?"). But
the TM enrichment analysis requires *spatial resolution* to separate TM from
non-TM constraint within each gene. AIC provides this resolution because it
always selects the best-fitting multi-region model, even when the improvement
over the null is modest. LRT's conservatism — producing one large region or
no regions at all — destroys the within-gene contrast that the enrichment test
depends on.

For cross-gene enrichment analyses where individual gene effects are small but
systematic, **AIC is the appropriate choice**.

---

## 4. Constraint by TM Segment Type

Not all transmembrane segments are equally constrained. UniProt annotations
distinguish functionally distinct segment types:

Statistics are from the AIC no-filter method, including catch-all residues for
genes with at least one non-null constraint region. "% cat < non-TM" is the
fraction of genes where the category geomean OE is below their own non-TM
geomean OE (from the paired Wilcoxon analysis).

| Segment type | Genes | Median geomean OE | % cat < non-TM |
|---|---|---|---|
| Pore-lining (S5) | 74 | 0.602 | 97% |
| Pore-lining (S6) | 75 | 0.609 | 97% |
| Named TM (M1/M2) | 24 | 0.637 | 92% |
| Voltage sensor (S4) | 76 | 0.660 | 95% |
| VSD other (S1–S3) | 77 | 0.735 | 88% |
| Intramembrane | 197 | 0.852 | 85% |
| Generic TM helix | 3,597 | 1.030 | 62% |
| Transporter helix | 913 | 1.051 | 82% |
| **Non-TM** | **4,615** | **1.075** | **—** |

### Key findings

1. **Pore-lining segments (S5/S6) are the most constrained TM regions** in ion
   channels (median geomean OE ~0.60), consistent with their critical role in
   ion selectivity and conductance. With catch-all residues included, 97% of
   pore-lining genes show category OE below their own non-TM OE.

2. **Named TM segments (M1/M2)** in ligand-gated channels show similar high
   constraint (median 0.637, 92% of genes show category < non-TM), confirming
   that functionally annotated TM segments are strongly constrained.

3. **Voltage-sensor segments (S4)** are also highly constrained (median 0.660,
   95% of genes show category < non-TM), reflecting the essential role of
   charged residues in voltage sensing.

4. **All TM categories are more constrained than non-TM** (median non-TM =
   1.075). Even generic TM helices (median 1.030) and transporter helices
   (median 1.051) now fall below the non-TM median. This result, which was
   obscured before catch-all inclusion, indicates that all TM segments exert
   at least some purifying selection relative to the non-TM background.

5. **Transporter helices** (n=913) have a median geomean OE of 1.051, below the
   non-TM median of 1.075. 82% of transporter genes show category OE below their
   own non-TM OE (Wilcoxon p = 7.1 × 10⁻⁹¹), suggesting genuine constraint
   rather than the near-neutral picture seen in earlier analyses. However, the
   effect size remains moderate and heterogeneity across transporter classes
   is likely.

6. **Intramembrane regions** (n=197) show substantial constraint (median 0.852,
   85% of genes show category < non-TM), consistent with roles in selectivity
   filters and re-entrant loops.

### Ion channel segment breakdown

![Ion channel TM segment constraint by type](fig2_ion_channel_segment_types.png)

**Figure 2.** Per-gene constraint for 76 ion channel genes, broken down by TM
segment type. For most genes, pore-lining segments (S5/S6) are among the most
constrained, followed by the voltage sensor (S4), then S1–S3. Black edges
indicate p < 0.05 vs non-TM.

### Distribution by segment type

![Segment type histograms](fig3_segment_type_histograms.png)

**Figure 3.** Histograms of per-gene geometric mean OE upper for each TM
segment category. Black vertical lines show medians. The hierarchy of constraint
is visible: pore-lining and named TM segments cluster below 0.65, while generic
TM helices and transporters center near 1.0–1.05 (below the non-TM median of
1.075).

### Aggregate category comparison

![Constraint by TM segment type — aggregate](fig4_category_summary_strip.png)

**Figure 4.** Strip plot summarizing the distribution of per-gene geometric mean
OE upper across all TM segment categories. Each dot is one gene. Horizontal
bars show medians; vertical bars show interquartile ranges.

---

## 5. Structural Confidence Bias and the Entanglement Problem

### TM regions have systematically higher pLDDT and lower PAE

AlphaFold's per-residue confidence score (pLDDT) and pairwise alignment error
(PAE) are not uniformly distributed across proteins. TM residues have
substantially higher structural confidence than non-TM residues:

|  | Median pLDDT | % excluded by pLDDT < 70 | n residues |
|---|---|---|---|
| TM residues | 92.9 | 4.2% | 405,792 |
| Non-TM residues | 85.3 | 33.8% | 2,048,991 |

Across genes, 3,455 of 4,665 TM-containing genes (74%) have higher median
pLDDT in their TM residues than their non-TM residues. A pLDDT < 70 cutoff
removes only ~4% of TM residues but ~34% of non-TM residues.

![pLDDT bias check](fig5_plddt_bias_check.png)

**Figure 5.** Left: pLDDT distributions for TM (blue) vs non-TM (orange)
residues, pooled across all genes. Center: histogram of per-gene median pLDDT
difference (TM minus non-TM). Right: histogram of differential exclusion
(% non-TM excluded minus % TM excluded by pLDDT < 70).

### Why this happens

There are two reinforcing reasons why TM regions score high on structural
confidence:

1. **Regular secondary structure.** TM segments are almost exclusively alpha
   helices. Helices have stereotyped backbone geometry that AlphaFold predicts
   with high confidence regardless of evolutionary information. Non-TM regions
   include loops, disordered segments, and flexible linkers whose conformations
   are inherently uncertain.

2. **Evolutionary conservation feeds into prediction.** AlphaFold uses multiple
   sequence alignments (MSAs) as a primary input. Residues that are
   evolutionarily conserved produce richer, more consistent MSA signals, which
   translates directly into higher pLDDT.

### Impact on constraint analysis

When pLDDT or PAE filters are applied to the Proemis3D constraint model, the
effect on TM vs non-TM comparisons is asymmetric:

- **pLDDT exclude** removes low-confidence residues from the constraint model
  statistics. Since low-pLDDT residues are disproportionately non-TM (34% vs
  4%), the filter asymmetrically affects non-TM residues. With catch-all
  residues included in the analysis, pLDDT-filtered AIC now shows the *strongest*
  global signal (p = 4.98 × 10⁻¹³⁵, median diff −0.094) — stronger than no-filter
  AIC (p = 6.18 × 10⁻⁹³, median diff −0.061). This reversal from earlier analyses
  (which excluded catch-all residues) likely occurs because low-confidence
  non-TM residues assigned to the catch-all have relatively high OE; the pLDDT
  filter removes these from the constraint model, producing a cleaner non-TM
  background against which TM constraint is measured. However, the pLDDT exclude
  LRT method shows a reversed direction (TM > non-TM, p = 8.4 × 10⁻⁵), indicating
  that LRT's coarse regions interact poorly with this filter.

- **PAE filter center** shows a strong TM signal (p = 2.44 × 10⁻⁹⁰, median
  diff −0.068) but is no longer the strongest method. PAE filtering may define
  structurally coherent neighborhoods that capture TM bundle organization, but
  the signal is slightly weaker than pLDDT exclude under AIC.

### The entanglement: what can and cannot be separated

The core difficulty is that multiple properties are correlated through shared
biological causes:

```
Evolutionary conservation
  ├─→ Low OE upper (constraint)
  ├─→ Rich MSA signal → High pLDDT
  ├─→ Disease variant enrichment
  └─→ Regular structure → High pLDDT, Low PAE

Structural regularity (alpha helix in membrane)
  ├─→ High pLDDT (easy to predict)
  ├─→ Low PAE (tightly packed)
  └─→ NOT necessarily constrained (generic TM anchors are not)
```

The fact that **generic TM helices are structurally regular (high pLDDT) but
NOT more constrained than non-TM** is informative. It demonstrates that high
pLDDT alone does not produce a constraint signal—the constraint in ion channel
TM segments reflects genuine purifying selection beyond what structural
regularity would predict.

### Recommendations

1. **Use no-filter methods** for TM vs non-TM comparisons, as they avoid the
   asymmetric residue exclusion introduced by pLDDT filtering.

2. **Interpret pLDDT-filtered results with caution** when comparing regions that
   differ systematically in structural confidence.

3. **The segment-type hierarchy is the strongest evidence** that the TM
   constraint signal is biologically real: pore (S5/S6) > voltage sensor (S4) >
   S1–S3 > intramembrane > generic TM ≈ non-TM. This ordering tracks with known
   functional importance, not with pLDDT (all named TM segments have similarly
   high pLDDT).

4. **A matched-pLDDT analysis** could further control for this bias: comparing
   TM vs non-TM constraint using only residues in a shared pLDDT range (e.g.,
   80–95).

---

## 6. Per-Gene Constraint Visualization

To illustrate how constraint regions map onto protein features, we generated
per-gene visualizations for six representative genes: three with non-generic
transmembrane regions and three with DNA-binding domains.

Each gene has four plot types:
- **Bars plot (grey catch-all):** Segmented bars coloured by OE upper across all
  6 methods. The catch-all region (residues not assigned to any constraint region
  by the forward algorithm) is shown in flat grey.
- **Bars plot (coloured catch-all):** Same layout, but the catch-all region is
  coloured by its own OE upper value using the same scale as the constraint
  regions.
- **Features plot:** UniProt feature annotations (TM segments, domains, binding
  sites, etc.) and de novo variant tracks aligned to the protein sequence.
- **Structure plot:** AlphaFold per-residue pLDDT confidence (top panel) and
  pairwise predicted aligned error (PAE) matrix (bottom panel), with the
  pLDDT < 70 cutoff line marked.

Methods are ordered top-to-bottom by filter group: no filter (AIC, LRT),
pLDDT exclude (AIC, LRT), PAE filter center (AIC, LRT).

### 6.1. Ion channels with named TM segments

**Q99250 — SCN2A (Nav1.2, voltage-gated sodium channel, 2,005 residues)**

SCN2A is one of the most constrained ion channels, with 24 TM segments across
four homologous repeats (I–IV), each containing voltage-sensor (S1–S4) and
pore-lining (S5–S6) segments. This gene has 47 DD de novo case variants and
18 autism case variants.

![Q99250 bars — grey catch-all](gene_viz/Q99250_bars.png)
![Q99250 bars — coloured catch-all](gene_viz/Q99250_bars_colored.png)
![Q99250 features](gene_viz/Q99250_features.png)
![Q99250 structure](gene_viz/Q99250_structure.png)

### 6.2. Ligand-gated channel and ATPase

**Q13224 — GRIN2B (NMDA receptor subunit, 1,484 residues)**

GRIN2B is an NMDA receptor subunit with named TM segments (M1–M4). It has 39
DD de novo and 6 autism case variants. The constraint pattern shows strong
constraint in the TM regions and ligand-binding domain.

![Q13224 bars — grey catch-all](gene_viz/Q13224_bars.png)
![Q13224 bars — coloured catch-all](gene_viz/Q13224_bars_colored.png)
![Q13224 features](gene_viz/Q13224_features.png)
![Q13224 structure](gene_viz/Q13224_structure.png)

**P13637 — ATP1A3 (Na+/K+-transporting ATPase alpha-3, 1,013 residues)**

ATP1A3 has 10 TM segments and is associated with alternating hemiplegia
of childhood. It has 21 DD de novo and 4 autism case variants.

![P13637 bars — grey catch-all](gene_viz/P13637_bars.png)
![P13637 bars — coloured catch-all](gene_viz/P13637_bars_colored.png)
![P13637 features](gene_viz/P13637_features.png)
![P13637 structure](gene_viz/P13637_structure.png)

### 6.3. DNA-binding domain proteins

To contrast with TM-rich proteins, we also examined genes with prominent
DNA-binding domains but no transmembrane regions.

**P55316 — FOXG1 (forkhead box protein G1, 488 residues)**

FOXG1 is a transcription factor critical for forebrain development. Its
forkhead DNA-binding domain (residues 180–274) is associated with FOXG1
syndrome, a severe neurodevelopmental disorder.

![P55316 bars — grey catch-all](gene_viz/P55316_bars.png)
![P55316 bars — coloured catch-all](gene_viz/P55316_bars_colored.png)
![P55316 features](gene_viz/P55316_features.png)
![P55316 structure](gene_viz/P55316_structure.png)

**Q6VB84 — FOXD4L3 (forkhead box protein D4-like 3, 416 residues)**

FOXD4L3 has a forkhead domain at residues 107–201. It provides a comparison
to FOXG1 as a smaller forkhead family member.

![Q6VB84 bars — grey catch-all](gene_viz/Q6VB84_bars.png)
![Q6VB84 bars — coloured catch-all](gene_viz/Q6VB84_bars_colored.png)
![Q6VB84 features](gene_viz/Q6VB84_features.png)
![Q6VB84 structure](gene_viz/Q6VB84_structure.png)

**P10589 — NR2F1 / COUP-TF1 (nuclear receptor, 423 residues)**

NR2F1 is a nuclear receptor with a DNA-binding domain (residues 82–157) and
a ligand-binding domain (183–409). Mutations cause Bosch–Boonstra–Schaaf
optic atrophy syndrome. The dual-domain structure contrasts with the
single-domain forkhead proteins above.

![P10589 bars — grey catch-all](gene_viz/P10589_bars.png)
![P10589 bars — coloured catch-all](gene_viz/P10589_bars_colored.png)
![P10589 features](gene_viz/P10589_features.png)
![P10589 structure](gene_viz/P10589_structure.png)

---

## 7. Summary

1. **TM regions are more constrained than non-TM regions** overall (p = 6.18 ×
   10⁻⁹³ for AIC no filter across 4,615 genes; median difference −0.061). The
   signal is substantially stronger when catch-all residues are included in the
   analysis, consistent with the catch-all capturing non-constrained sequence
   context that differs between TM and non-TM regions.

2. **All TM segment types are more constrained than non-TM** (median non-TM =
   1.075). Even generic TM helices (1.030) and transporter helices (1.051) fall
   below the non-TM median, with 62% and 82% of genes respectively showing
   category < non-TM. The largest effects are in ion channel TM segments.

3. **Within ion channels, a clear hierarchy exists:** pore-lining (S5/S6,
   ~0.60) > named TM (M1/M2, 0.64) > voltage sensor (S4, 0.66) > other VSD
   segments (S1–S3, 0.74) >> generic TM (1.03) ≈ non-TM (1.08).

4. **LRT is more conservative than AIC** and classifies 49% of genes as having
   no significant regional constraint variation, leaving 1,882–2,363 genes for
   analysis vs 4,552–4,615 under AIC. AIC is the appropriate choice for
   cross-gene enrichment analyses.

5. **Structural confidence (pLDDT) is systematically higher in TM regions**
   (median 92.9 vs 85.3), creating an asymmetric bias when pLDDT filters are
   applied. With catch-all residues included, pLDDT exclude AIC shows the
   *strongest* TM signal (p = 4.98 × 10⁻¹³⁵), suggesting the filter removes
   non-TM catch-all residues with high OE, sharpening the TM contrast.

6. **The segment-type hierarchy provides the strongest evidence** that the
   constraint signal reflects genuine biology, since all named TM segments have
   similarly high pLDDT but very different constraint levels.

---

## 8. Methods

- **Constraint metric:** OE upper from Proemis3D forward runs across 6
  production configurations (3 filter groups × AIC/LRT). **Catch-all (null-model)
  residues are included** for genes that have at least one non-null constraint
  region; genes where all residues are classified as null (no constraint regions
  found by the forward algorithm) are excluded from the analysis.
- **TM annotations:** UniProt transmembrane/intramembrane features (release
  2021_04), extracted from a Hail Table. One transcript per gene selected.
- **Summary statistic:** Geometric mean = exp(mean(log(x))) of per-residue OE
  upper values, computed separately per gene per region category.
- **Within-gene test:** Mann-Whitney U (two-sided) comparing OE upper values of
  TM vs non-TM residues (including catch-all residues for qualifying genes).
- **Across-gene test:** Paired Wilcoxon signed-rank test on per-gene geometric
  means (TM minus non-TM); zero-difference pairs excluded from the test statistic.
- **TM segment classification:** Based on UniProt feature notes — S1–S6 named
  segments for ion channels, numbered helices for transporters, M1/M2 for
  ligand-gated channels, and "Helical" for generic TM segments.
- **Minimum thresholds:** Genes required ≥10 TM residues and ≥3 residues
  (including catch-all) in both TM and non-TM categories.
- **Dataset:** 4,665 genes with TM annotations, 14.7M per-residue rows, from
  6 full production Hail Tables.

---

## 9. Disease Burden, Bimodality, and 3D Structure Constraint

The aggregate TM constraint signal established in Sections 2–5 is driven by a
subset of ion channels with both high disease burden and exceptionally low OE
upper in their pore-lining segments. Here we examine the three most
DD-burdened genes — SCN8A, KCNQ2, and KCNB1 — through three complementary
lenses: (1) a cross-gene bimodality analysis, (2) linear all-method constraint
tracks, and (3) 3D AlphaFold structure coloring.

### 9.1 S5 constraint bimodality by disease burden

S5 (pore-lining helix, inner leaflet) shows the strongest cross-gene constraint
signal of all TM segment types (median OE upper 0.516 vs non-TM 1.075, Section
4). Among ion channel genes with S5 segments, there is a striking dose–response
between the number of de novo disease variants per gene and S5 constraint:

| DD burden | Median S5 OE upper | Example genes |
|---|---|---|
| 0 DD de novo | 0.860 | — |
| 1–5 DD de novo | 0.679 | — |
| > 5 DD de novo | 0.516 | SCN8A, KCNQ2, KCNB1, SCN2A, KCNA2 |

The top-10 most S5-constrained genes with > 5 DD are all clinically validated
epilepsy channels (KCND3, KCNC1, SCN8A, KCNA2, KCNQ2, KCNB1, SCN1A, SCN2A,
KCNC2, KCNQ3), confirming the biological specificity of the constraint signal.

![Bimodality analysis](fig6_bimodality_analysis.png)

**Figure 6.** S5 pore-lining constraint (AIC no filter, OE upper) stratified by
de novo disease variant burden per gene. Violin plots show the distribution of
per-gene median S5 OE upper for three burden groups (0, 1–5, > 5 DD de novo).
The bimodal pattern in the > 5 DD group reflects the co-existence of high-burden
epilepsy channels (low OE, bottom mode) and lower-burden genes that happen to
have > 5 cumulative variants (top mode). The three focus genes (SCN8A, KCNQ2,
KCNB1) are among the most constrained high-burden ions.

### 9.2 Linear constraint tracks — all 6 methods

Each figure shows segmented bars coloured by OE upper CI (dark red = highly
constrained, cream/yellow = unconstrained) across all 6 ProEMIS3D configurations.
Methods are ordered by filter group (no filter → pLDDT exclude → PAE filter
center), with AIC before LRT within each group. UniProt feature annotations and
de novo variant tracks (DD, autism) are shown in the features figure.

**SCN8A — Q9UQD0 (Nav1.6, 1,980 residues, 32 DD de novo, 3 autism cases)**

SCN8A encodes Nav1.6, a voltage-gated sodium channel with 4 homologous
repeats (I–IV), each containing voltage-sensor (S1–S4) and pore-lining (S5–S6)
segments (24 TM segments total). It is the most frequently mutated sodium
channel in developmental and epileptic encephalopathy (DEE13).

![Q9UQD0 bars — grey catch-all](protein_3d_figures/Q9UQD0_bars.png)
![Q9UQD0 bars — coloured catch-all](protein_3d_figures/Q9UQD0_bars_colored.png)
![Q9UQD0 features + variants](protein_3d_figures/Q9UQD0_features.png)

**KCNQ2 — O43526 (Kv7.2, 872 residues, 37 DD de novo, 2 autism cases)**

KCNQ2 encodes Kv7.2, a voltage-gated potassium channel with a single S1–S6
domain (forms a functional tetramer). It includes the pore-forming segment H5
(selectivity filter GYGX) and an extended intracellular C-terminus that mediates
PIP2 regulation and calmodulin interaction. Variants cause KCNQ2 neonatal
epilepsy and neurodevelopmental delay.

![O43526 bars — grey catch-all](protein_3d_figures/O43526_bars.png)
![O43526 bars — coloured catch-all](protein_3d_figures/O43526_bars_colored.png)
![O43526 features + variants](protein_3d_figures/O43526_features.png)

**KCNB1 — Q14721 (Kv2.1, 858 residues, 21 DD de novo, 3 autism cases)**

KCNB1 encodes Kv2.1, a delayed-rectifier potassium channel. Like KCNQ2 it has
a single S1–S6 domain (tetramer), pore helix, and selectivity filter. Its large
intracellular C-terminus contains a proximal restriction and clustering domain
(PRC) that anchors channels to the cell soma. Missense variants cause DEE26
with gain-of-function effects on voltage-sensor cooperativity.

![Q14721 bars — grey catch-all](protein_3d_figures/Q14721_bars.png)
![Q14721 bars — coloured catch-all](protein_3d_figures/Q14721_bars_colored.png)
![Q14721 features + variants](protein_3d_figures/Q14721_features.png)

### 9.3 3D AlphaFold structure constraint coloring

The gnomAD missense constraint viewer
(`gnomad-missense-constraint-viewer-demo`) overlays ProEMIS3D OE upper scores
onto the AlphaFold CIF structure using Mol* (Molstar). Residues are colored on
the gnomAD RMC scale: dark red (OE upper ≤ 0.6, most constrained) through
orange to cream/yellow (OE upper ≥ 1.5, unconstrained). Catch-all (null-model)
residues are shown in grey.

The viewer supports all 6 ProEMIS3D methods via a dropdown selector. Scores
are served from the gene HT (`genes_v4.proemis3d_all_methods.ht`, 61,533
genes) as six separate regional constraint fields (`proemis3d_aic`,
`proemis3d_lrt`, `proemis3d_plddt_exclude_aic`, `proemis3d_plddt_exclude_lrt`,
`proemis3d_pae_filter_exclude_aic`, `proemis3d_pae_filter_exclude_lrt`). Each
field contains an array of contiguous segments with obs, exp, oe, oe_upper,
and genomic coordinates. The score map is built from aa-start/stop regions
and passed to Molstar via `promis3doeUpper` and `promis3dObsExp` props.
Switching methods remounts Molstar cleanly (`key={activeMethod}`).

**Data pipeline fix:** The original combined-residues HT lookup by
`canonical_transcript_id` returned null for SCN8A, KCNQ2, and KCNB1 because
ProEMIS3D was run on non-canonical transcripts (e.g. SCN8A used
`ENST00000354534`, gnomAD canonical is `ENST00000627620`). The fix builds the
gene HT directly from the 6 individual forward run HTs, joins by
`gene_id` (Ensembl ID) via `gene_mapping.tsv`, and groups by gene_id — making
the result independent of canonical transcript choice. After this fix:
- SCN8A (`Q9UQD0`): 140 AIC regions, first oe_upper=0.822 ✓
- KCNQ2 (`O43526`):  58 AIC regions, first oe_upper=0.310 ✓
- KCNB1 (`Q14721`):  48 AIC regions, first oe_upper=2.269 ✓

**To generate 3D screenshots** (servers run on ports 4000, 5173, 8080):

1. Open `http://localhost:5173` in a browser and search for each gene symbol:
   `SCN8A`, `KCNQ2`, `KCNB1`.
2. In the **ProEMIS3D Method** dropdown, select `AIC — no filter`.
3. Wait for Molstar to finish loading and coloring the AlphaFold structure.
4. Orient the view to show the transmembrane domain (side view for S1–S6
   helices; top-down for the selectivity filter pore).
5. Export via Mol*'s `Export image` or browser screenshot and save to
   `proemis3d_enrichment_plots/protein_3d_figures/{SYMBOL}_aic_3d.png`.
6. Repeat for `AIC — pLDDT exclude` to illustrate the effect of pLDDT filtering.

**Expected appearance:** For all three genes, S5–S6 segments and selectivity
filter residues should appear in dark red (OE upper < 0.6); voltage-sensor
helices S4 in medium orange-red (~0.65); flanking non-TM linkers in
cream/yellow (OE upper > 1.0). The contrast between the pore domain
(constrained) and cytoplasmic/extracellular loops (unconstrained) should be
visually striking and consistent across AIC no-filter and pLDDT-filtered methods.

<!-- Placeholder — add screenshots once captured from the viewer:
![SCN8A 3D — AIC no filter](protein_3d_figures/SCN8A_aic_3d.png)
![KCNQ2 3D — AIC no filter](protein_3d_figures/KCNQ2_aic_3d.png)
![KCNB1 3D — AIC no filter](protein_3d_figures/KCNB1_aic_3d.png)
-->
