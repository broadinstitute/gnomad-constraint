# Matched pLoF o/e: Detailed Numbers

## Per-Percentile Aggregates (AlphaMissense, MANE Select, All Genes)

| Pctl | mis_pos | mis_obs | mis_exp | mis_exp_adjr | mis_oe | plof_obs_xE | plof_exp_xE | plof_oe_E | plof_obs_xP | plof_exp_xP | plof_oe_P | n_w_plof |
|-----:|--------:|--------:|--------:|-------------:|-------:|------------:|------------:|----------:|------------:|------------:|----------:|---------:|
| 1 | 661,701 | 187,172 | 148,750 | 163,620 | 1.1439 | 9,434,500 | 19,528,000 | 0.4831 | 36,986,470 | 75,998,000 | 0.4867 | 650,313 |
| 10 | 649,206 | 155,181 | 147,260 | 160,160 | 0.9689 | 10,158,000 | 19,718,000 | 0.5152 | 41,141,614 | 78,463,000 | 0.5243 | 633,377 |
| 25 | 629,492 | 129,472 | 129,340 | 140,450 | 0.9218 | 9,143,700 | 17,624,000 | 0.5188 | 41,667,508 | 78,297,000 | 0.5322 | 609,209 |
| 50 | 646,331 | 116,575 | 124,170 | 133,950 | 0.8703 | 8,655,800 | 17,482,000 | 0.4951 | 42,812,771 | 83,975,000 | 0.5098 | 626,873 |
| 75 | 693,041 | 108,147 | 131,850 | 141,030 | 0.7668 | 8,366,500 | 17,955,000 | 0.4660 | 42,050,609 | 88,245,000 | 0.4765 | 681,665 |
| 90 | 713,061 | 88,934 | 133,990 | 143,180 | 0.6211 | 7,822,300 | 18,061,000 | 0.4331 | 39,524,003 | 89,738,000 | 0.4404 | 708,148 |
| 99 | 722,381 | 53,004 | 136,620 | 145,770 | 0.3636 | 5,683,300 | 16,981,000 | 0.3347 | 28,583,904 | 84,984,000 | 0.3363 | 721,705 |
| 100 | 723,060 | 41,867 | 140,340 | 151,010 | 0.2772 | 4,555,500 | 15,663,000 | 0.2908 | 21,961,601 | 75,436,000 | 0.2911 | 722,711 |

### Column key

| Column | Description |
|--------|-------------|
| mis_pos | Possible missense variants in this bin |
| mis_obs | Observed missense count in this bin |
| mis_exp | Expected missense (raw, no adj_r) |
| mis_exp_adjr | Expected missense (with adj_r correction) |
| mis_oe | mis_obs / mis_exp_adjr (solid line in plots) |
| plof_obs_xE | sum_genes(gene_plof_obs * gene_exp_adjr_missense_in_bin) |
| plof_exp_xE | sum_genes(gene_plof_exp * gene_exp_adjr_missense_in_bin) |
| plof_oe_E | plof_obs_xE / plof_exp_xE (dashed line, weight = expected) |
| plof_obs_xP | sum_genes(gene_plof_obs * gene_possible_missense_in_bin) |
| plof_exp_xP | sum_genes(gene_plof_exp * gene_possible_missense_in_bin) |
| plof_oe_P | plof_obs_xP / plof_exp_xP (dashed line, weight = possible) |
| n_w_plof | Variants with defined gene-level pLoF data |

## Example Genes (constraint_metrics, no adj_r)

| Gene | pLoF obs | pLoF exp | pLoF o/e | mis obs | mis exp | mis o/e | ~mis pos |
|-----:|---------:|---------:|---------:|--------:|--------:|--------:|---------:|
| TTN | 1,461 | 3,180 | 0.459 | 40,245 | 44,300 | 0.908 | ~88K |
| BRCA2 | 188 | 264 | 0.711 | 3,781 | 3,960 | 0.956 | ~8K |
| SCN1A | 9 | 160 | 0.056 | 1,291 | 2,290 | 0.564 | ~4K |

## How the Matched pLoF Works

For each percentile bin, each gene *g* contributes missense variants. The gene's **weight** in the matched pLoF calculation is its adj_r-corrected expected missense count in that bin (or the number of possible missense variants, depending on the weighting approach).

```
plof_obs_xE(bin) = sum over genes g of: gene_plof_obs_g * weight_g
plof_exp_xE(bin) = sum over genes g of: gene_plof_exp_g * weight_g

matched_plof_oe(bin) = plof_obs_xE / plof_exp_xE
```

### Worked example for P100

If **SCN1A** contributes weight = 20 to P100:
- Adds 9 * 20 = 180 to plof_obs_xE
- Adds 160 * 20 = 3,200 to plof_exp_xE
- SCN1A's contribution ratio: 180 / 3,200 = 0.056 (its gene-level pLoF o/e)

If **TTN** contributes weight = 5 to P100:
- Adds 1,461 * 5 = 7,305 to plof_obs_xE
- Adds 3,180 * 5 = 15,900 to plof_exp_xE
- TTN's contribution ratio: 7,305 / 15,900 = 0.459 (its gene-level pLoF o/e)

The final matched pLoF o/e for the bin is the ratio of the summed numerator to the summed denominator across **all** genes. Highly constrained genes like SCN1A (pLoF o/e = 0.06) pull it down; unconstrained genes like TTN (pLoF o/e = 0.46) pull it up.

## Per-Gene Per-Percentile Breakdown (AlphaMissense, MANE Select)

pLoF obs/exp below are computed from the per-SNV table (LOFTEE HC, `possible_variants == 1`, adj_r-corrected exp), which differs slightly from constraint_metrics values above.

### TTN (ENST00000589042)
- pLoF obs = 1,525, pLoF exp = 3,214.5, pLoF o/e = 0.474

| Pctl | AM_pos | AM_obs | AM_exp_adjr | plof_obs_xE | plof_exp_xE | plof_obs_xP | plof_exp_xP | AM_oe | plof_oe_E | plof_oe_P |
|-----:|-------:|-------:|------------:|------------:|------------:|------------:|------------:|------:|----------:|----------:|
| 1 | 983 | 304 | 272 | 414,205 | 873,090 | 1,499,075 | 3,159,854 | 1.119 | 0.474 | 0.474 |
| 10 | 1,298 | 311 | 314 | 479,064 | 1,009,803 | 1,979,450 | 4,172,421 | 0.990 | 0.474 | 0.474 |
| 25 | 1,879 | 396 | 399 | 607,774 | 1,281,107 | 2,865,475 | 6,040,046 | 0.994 | 0.474 | 0.474 |
| 50 | 3,005 | 506 | 569 | 867,908 | 1,829,436 | 4,582,625 | 9,659,572 | 0.889 | 0.474 | 0.474 |
| 75 | 2,616 | 432 | 469 | 715,301 | 1,507,761 | 3,989,400 | 8,409,132 | 0.921 | 0.474 | 0.474 |
| 90 | 2,378 | 385 | 455 | 693,341 | 1,461,472 | 3,626,450 | 7,644,081 | 0.847 | 0.474 | 0.474 |
| 99 | 1,536 | 217 | 274 | 417,576 | 880,194 | 2,342,400 | 4,937,472 | 0.792 | 0.474 | 0.474 |
| 100 | 647 | 76 | 138 | 210,267 | 443,215 | 986,675 | 2,079,782 | 0.551 | 0.474 | 0.474 |

TTN is a very large gene (~88K possible missense) with moderate pLoF constraint. Missense depletion is mild across most percentile bins (o/e ~0.85-1.0), with only the most pathogenic bin (P100) showing substantial depletion (0.55). TTN dominates the aggregate matched pLoF calculation due to its large weight (e.g., plof_exp_xP ~10M at P50 vs SCN1A's ~19K).

### BRCA2 (ENST00000380152)
- pLoF obs = 203, pLoF exp = 256.1, pLoF o/e = 0.793

| Pctl | AM_pos | AM_obs | AM_exp_adjr | plof_obs_xE | plof_exp_xE | plof_obs_xP | plof_exp_xP | AM_oe | plof_oe_E | plof_oe_P |
|-----:|-------:|-------:|------------:|------------:|------------:|------------:|------------:|------:|----------:|----------:|
| 1 | 220 | 70 | 53.9 | 10,945 | 13,811 | 44,660 | 56,353 | 1.298 | 0.793 | 0.793 |
| 10 | 431 | 95 | 86.5 | 17,567 | 22,166 | 87,493 | 110,401 | 1.098 | 0.793 | 0.793 |
| 25 | 436 | 77 | 70.2 | 14,253 | 17,985 | 88,508 | 111,681 | 1.097 | 0.793 | 0.793 |
| 50 | 172 | 27 | 26.2 | 5,319 | 6,712 | 34,916 | 44,058 | 1.030 | 0.793 | 0.793 |
| 75 | 95 | 21 | 17.9 | 3,627 | 4,576 | 19,285 | 24,334 | 1.175 | 0.793 | 0.793 |
| 90 | 44 | 8 | 6.6 | 1,342 | 1,694 | 8,932 | 11,271 | 1.210 | 0.793 | 0.793 |
| 99 | 0 | 0 | 0.0 | 0 | 0 | 0 | 0 | — | — | — |
| 100 | 0 | 0 | 0.0 | 0 | 0 | 0 | 0 | — | — | — |

BRCA2 shows essentially no missense depletion at any percentile (o/e >= 1.0 everywhere). Its constraint is primarily LoF-specific. P99-P100 have zero AlphaMissense-scored variants, so BRCA2 contributes nothing to the aggregate matched pLoF in those bins. Its moderate pLoF o/e (0.79) pulls the aggregate toward neutrality in bins where it has weight.

### SCN1A (ENST00000674923)
- pLoF obs = 13, pLoF exp = 183.5, pLoF o/e = 0.071

| Pctl | AM_pos | AM_obs | AM_exp_adjr | plof_obs_xE | plof_exp_xE | plof_obs_xP | plof_exp_xP | AM_oe | plof_oe_E | plof_oe_P |
|-----:|-------:|-------:|------------:|------------:|------------:|------------:|------------:|------:|----------:|----------:|
| 1 | 16 | 2 | 3.5 | 45 | 637 | 208 | 2,936 | 0.576 | 0.071 | 0.071 |
| 10 | 36 | 7 | 8.7 | 113 | 1,599 | 468 | 6,606 | 0.803 | 0.071 | 0.071 |
| 25 | 40 | 10 | 9.5 | 123 | 1,742 | 520 | 7,340 | 1.054 | 0.071 | 0.071 |
| 50 | 101 | 13 | 15.9 | 207 | 2,922 | 1,313 | 18,534 | 0.816 | 0.071 | 0.071 |
| 75 | 131 | 13 | 24.8 | 322 | 4,551 | 1,703 | 24,038 | 0.524 | 0.071 | 0.071 |
| 90 | 290 | 23 | 54.7 | 712 | 10,045 | 3,770 | 53,215 | 0.420 | 0.071 | 0.071 |
| 99 | 594 | 24 | 104 | 1,347 | 19,007 | 7,722 | 108,999 | 0.232 | 0.071 | 0.071 |
| 100 | 518 | 5 | 91.0 | 1,183 | 16,704 | 6,734 | 95,053 | 0.055 | 0.071 | 0.071 |

SCN1A is the classic example: extreme pLoF constraint (o/e = 0.07) paired with dramatic missense depletion that tracks AlphaMissense percentile perfectly — from o/e ~0.8-1.0 at benign bins down to 0.055 at P100. Despite its strong signal, SCN1A's weight in the aggregate is tiny compared to TTN (e.g., plof_exp_xP ~19K at P50 vs TTN's ~10M), illustrating how large genes dominate the matched pLoF calculation.

## Fraction-Weighted Matched pLoF (AlphaMissense, MANE Select, All Genes)

Instead of weighting each gene by its raw missense count in the bin, weight by the **fraction** of the gene's total missense that falls in the bin. This removes gene-size bias: a small gene that concentrates variants in P100 gets proportionally more influence.

```
w_g(bin) = AM_exp_adjr_g(bin) / Σ_all_bins(AM_exp_adjr_g)     [frac_E]
w_g(bin) = AM_pos_g(bin) / Σ_all_bins(AM_pos_g)               [frac_P]
```

| Pctl | mis_oe | plof_oe_E (raw) | plof_oe_FE (frac) | plof_oe_P (raw) | plof_oe_FP (frac) |
|-----:|-------:|----------------:|------------------:|----------------:|------------------:|
| 1 | 1.144 | 0.4831 | 0.5360 | 0.4867 | 0.5301 |
| 10 | 0.969 | 0.5152 | 0.5718 | 0.5243 | 0.5705 |
| 25 | 0.922 | 0.5188 | 0.5761 | 0.5322 | 0.5792 |
| 50 | 0.870 | 0.4951 | 0.5519 | 0.5098 | 0.5571 |
| 75 | 0.767 | 0.4660 | 0.5191 | 0.4765 | 0.5221 |
| 90 | 0.621 | 0.4331 | 0.4805 | 0.4404 | 0.4817 |
| 99 | 0.364 | 0.3347 | 0.3619 | 0.3363 | 0.3601 |
| 100 | 0.277 | 0.2908 | 0.3121 | 0.2911 | 0.3112 |

| Weighting | P1 | P100 | Delta |
|-----------|------|------|-------|
| raw_E | 0.483 | 0.291 | -0.192 |
| frac_E | 0.536 | 0.312 | -0.224 |
| raw_P | 0.487 | 0.291 | -0.196 |
| frac_P | 0.530 | 0.311 | -0.219 |

Fraction weighting starts higher (~0.53 vs ~0.48) and drops further (-0.22 vs -0.19). Removing gene-size bias gives constrained small genes (like SCN1A) more relative influence in high-pathogenicity bins, steepening the matched pLoF decline.

## Global Reference Lines

| Metric | Source | Value |
|--------|--------|------:|
| syn o/e | per-SNV table, adj_r-corrected, MANE Select | 0.91 |
| pLoF o/e | per-SNV table, adj_r-corrected, MANE Select | 0.53 |
| syn o/e | constraint_metrics, no adj_r, MANE Select | 1.02 |
| mis o/e | constraint_metrics, no adj_r, MANE Select | 0.90 |
| pLoF o/e | constraint_metrics, no adj_r, MANE Select | 0.58 |
