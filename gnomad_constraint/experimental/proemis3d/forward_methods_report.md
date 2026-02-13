# Proemis3D forward method comparison (158 test genes)

## Summary table (all runs)

| Model comparison method | pLDDT filter | PAE filter | % Assigned missing | Const. regions (mean / median / min / max) | Region length (mean / median / min / max) | Null frac (mean) | Tx all-null % |
|-------------------------|--------------|------------|--------------------|-----------------------------------------|------------------------------------------|------------------|----------------|
| AIC | – | – | 0.00% | 9.2 / 9 / 1 / 29 | 135.7 / 42 / 8 / 2591 | 6.9% | 0.0% |
| AIC | – | truncate_center | 9.33% | 26.4 / 22 / 2 / 80 | 29.5 / 17 / 8 / 930 | 30.9% | 0.0% |
| AIC | – | filter_center | 0.22% | 29.2 / 24 / 2 / 100 | 34.0 / 18 / 8 / 1081 | 24.2% | 0.0% |
| AIC | – | filter_region | 12.56% | 21.8 / 17 / 1 / 82 | 33.1 / 17 / 8 / 655 | 30.9% | 0.0% |
| AIC | – | exclude_center | 0.03% | 6.7 / 5 / 1 / 24 | 192.3 / 53 / 10 / 2115 | 3.0% | 0.0% |
| AIC | – | exclude_region | 3.09% | 7.9 / 6 / 1 / 33 | 150.2 / 37 / 8 / 1814 | 6.3% | 0.0% |
| AIC | exclude | – | 0.00% | 7.2 / 6 / 1 / 28 | 152.0 / 53 / 8 / 2028 | 17.7% | 0.0% |
| AIC | exclude | truncate_center | 35.84% | 16.6 / 14 / 1 / 70 | 37.0 / 24 / 8 / 512 | 17.9% | 0.0% |
| AIC | exclude | filter_center | 32.72% | 12.0 / 10 / 1 / 53 | 57.8 / 30 / 8 / 1007 | 15.0% | 0.0% |
| AIC | exclude | filter_region | 35.80% | 14.0 / 12 / 1 / 67 | 42.8 / 26 / 8 / 508 | 18.6% | 0.0% |
| AIC | exclude | exclude_center | 10.04% | 3.3 / 3 / 1 / 9 | 308.2 / 168 / 10 / 2115 | 11.5% | 0.0% |
| AIC | exclude | exclude_region | 0.11% | 3.0 / 2 / 1 / 18 | 366.1 / 206 / 11 / 2139 | 15.7% | 0.0% |
| AIC | remove | – | 46.66% | 7.2 / 6 / 1 / 28 | 92.0 / 32 / 8 / 1599 | 3.6% | 0.0% |
| AIC | remove | truncate_center | 46.94% | 14.0 / 11 / 1 / 68 | 41.2 / 23 / 8 / 984 | 9.7% | 0.0% |
| AIC | remove | filter_center | 46.82% | 12.0 / 10 / 1 / 53 | 49.8 / 24 / 8 / 971 | 8.1% | 0.0% |
| AIC | remove | filter_region | 46.98% | 14.4 / 12 / 1 / 67 | 39.7 / 23 / 8 / 498 | 10.1% | 0.0% |
| AIC | remove | exclude_center | 46.68% | 3.3 / 3 / 1 / 9 | 207.4 / 96 / 10 / 1693 | 1.1% | 0.0% |
| AIC | remove | exclude_region | 46.67% | 3.0 / 2 / 1 / 18 | 229.4 / 102 / 11 / 1457 | 1.0% | 0.0% |
| AIC weight | – | – | 0.00% | 4.3 / 4 / 1 / 15 | 179.3 / 59 / 8 / 2591 | 40.4% | 0.0% |
| AIC weight | – | truncate_center | 9.33% | 6.1 / 5 / 0 / 21 | 56.5 / 27 / 8 / 930 | 65.3% | 4.4% |
| AIC weight | – | filter_center | 0.22% | 5.7 / 4 / 0 / 25 | 77.7 / 27 / 8 / 1081 | 64.5% | 2.5% |
| AIC weight | – | filter_region | 12.56% | 5.6 / 4 / 0 / 25 | 61.9 / 29 / 8 / 655 | 59.7% | 3.8% |
| AIC weight | – | exclude_center | 0.03% | 6.1 / 4 / 1 / 22 | 210.2 / 66 / 11 / 2115 | 4.1% | 0.0% |
| AIC weight | – | exclude_region | 3.09% | 6.6 / 4 / 1 / 31 | 177.3 / 54 / 12 / 1814 | 7.8% | 0.0% |
| AIC weight | exclude | – | 0.00% | 3.3 / 3 / 0 / 14 | 178.1 / 65 / 8 / 1666 | 56.0% | 4.4% |
| AIC weight | exclude | truncate_center | 35.84% | 4.4 / 4 / 0 / 22 | 57.5 / 34 / 8 / 512 | 46.6% | 7.6% |
| AIC weight | exclude | filter_center | 32.72% | 3.8 / 3 / 0 / 20 | 87.8 / 40 / 8 / 1007 | 43.9% | 6.3% |
| AIC weight | exclude | filter_region | 35.80% | 4.1 / 3 / 0 / 22 | 64.7 / 34 / 8 / 508 | 45.3% | 6.3% |
| AIC weight | exclude | exclude_center | 10.04% | 3.1 / 3 / 1 / 9 | 324.8 / 188 / 14 / 2115 | 12.2% | 0.0% |
| AIC weight | exclude | exclude_region | 0.11% | 2.8 / 2 / 1 / 16 | 385.6 / 237 / 11 / 2139 | 16.6% | 0.0% |
| AIC weight | remove | – | 46.66% | 3.3 / 3 / 0 / 14 | 119.0 / 46 / 8 / 1100 | 25.5% | 4.4% |
| AIC weight | remove | truncate_center | 46.94% | 4.0 / 3 / 0 / 19 | 68.8 / 34 / 8 / 984 | 34.1% | 7.0% |
| AIC weight | remove | filter_center | 46.82% | 3.8 / 3 / 0 / 20 | 77.6 / 35 / 8 / 971 | 32.8% | 6.3% |
| AIC weight | remove | filter_region | 46.98% | 4.2 / 3 / 0 / 22 | 63.4 / 32 / 8 / 498 | 34.2% | 5.1% |
| AIC weight | remove | exclude_center | 46.68% | 3.1 / 3 / 1 / 9 | 219.5 / 115 / 14 / 1693 | 1.5% | 0.0% |
| AIC weight | remove | exclude_region | 46.67% | 2.8 / 2 / 1 / 16 | 242.4 / 120 / 11 / 1457 | 1.4% | 0.0% |
| LRT | – | – | 0.00% | 3.0 / 3 / 0 / 11 | 204.3 / 90 / 11 / 2591 | 52.3% | 8.2% |
| LRT | – | truncate_center | 9.33% | 4.1 / 3 / 0 / 17 | 68.6 / 35 / 10 / 930 | 70.3% | 11.4% |
| LRT | – | filter_center | 0.22% | 3.7 / 3 / 0 / 15 | 101.6 / 38 / 10 / 1081 | 69.9% | 8.9% |
| LRT | – | filter_region | 12.56% | 3.7 / 3 / 0 / 14 | 77.0 / 40 / 10 / 655 | 64.7% | 12.0% |
| LRT | – | exclude_center | 0.03% | 5.8 / 4 / 1 / 21 | 219.0 / 71 / 12 / 2115 | 4.6% | 0.0% |
| LRT | – | exclude_region | 3.09% | 6.2 / 4 / 1 / 29 | 186.9 / 59 / 12 / 1814 | 8.3% | 0.0% |
| LRT | exclude | – | 0.00% | 2.2 / 2 / 0 / 12 | 211.5 / 84 / 12 / 1666 | 65.8% | 15.2% |
| LRT | exclude | truncate_center | 35.84% | 2.4 / 2 / 0 / 15 | 70.0 / 43 / 12 / 512 | 53.1% | 22.2% |
| LRT | exclude | filter_center | 32.72% | 2.3 / 2 / 0 / 10 | 113.2 / 54 / 12 / 1007 | 49.8% | 18.4% |
| LRT | exclude | filter_region | 35.80% | 2.4 / 2 / 0 / 12 | 81.6 / 44 / 12 / 508 | 50.4% | 19.0% |
| LRT | exclude | exclude_center | 10.04% | 3.0 / 3 / 0 / 8 | 333.1 / 193 / 14 / 2115 | 12.7% | 0.6% |
| LRT | exclude | exclude_region | 0.11% | 2.7 / 2 / 0 / 15 | 398.2 / 244 / 16 / 2139 | 17.7% | 1.3% |
| LRT | remove | – | 46.66% | 2.2 / 2 / 0 / 12 | 145.4 / 62 / 12 / 1100 | 31.8% | 15.2% |
| LRT | remove | truncate_center | 46.94% | 2.4 / 2 / 0 / 12 | 85.5 / 42 / 12 / 984 | 39.9% | 19.0% |
| LRT | remove | filter_center | 46.82% | 2.3 / 2 / 0 / 10 | 100.5 / 46 / 12 / 971 | 38.0% | 18.4% |
| LRT | remove | filter_region | 46.98% | 2.4 / 2 / 0 / 11 | 79.9 / 42 / 11 / 498 | 39.9% | 17.7% |
| LRT | remove | exclude_center | 46.68% | 3.0 / 3 / 0 / 8 | 225.5 / 121 / 14 / 1693 | 1.8% | 0.6% |
| LRT | remove | exclude_region | 46.67% | 2.7 / 2 / 0 / 15 | 253.0 / 129 / 16 / 1457 | 1.8% | 1.3% |

## Summary table (one transcript per gene, 144 genes)

| Model comparison method | pLDDT filter | PAE filter | % Assigned missing | Const. regions (mean / median / min / max) | Region length (mean / median / min / max) | Null frac (mean) | Tx all-null % |
|-------------------------|--------------|------------|--------------------|-----------------------------------------|------------------------------------------|------------------|----------------|
| AIC | – | – | 0.00% | 9.4 / 9 / 1 / 29 | 137.6 / 42 / 8 / 2591 | 6.8% | 0.0% |
| AIC | – | truncate_center | 9.73% | 27.3 / 23 / 2 / 80 | 29.4 / 17 / 8 / 930 | 31.1% | 0.0% |
| AIC | – | filter_center | 0.23% | 30.5 / 24 / 2 / 100 | 33.8 / 18 / 8 / 1081 | 24.2% | 0.0% |
| AIC | – | filter_region | 12.86% | 22.8 / 18 / 1 / 82 | 32.7 / 17 / 8 / 655 | 31.1% | 0.0% |
| AIC | – | exclude_center | 0.03% | 6.9 / 5 / 1 / 24 | 192.8 / 52 / 10 / 2115 | 3.0% | 0.0% |
| AIC | – | exclude_region | 3.20% | 8.3 / 6 / 1 / 33 | 148.5 / 37 / 8 / 1814 | 6.5% | 0.0% |
| AIC | exclude | – | 0.00% | 7.3 / 6 / 1 / 28 | 156.0 / 54 / 8 / 2028 | 18.2% | 0.0% |
| AIC | exclude | truncate_center | 36.87% | 17.0 / 15 / 1 / 70 | 36.9 / 24 / 8 / 512 | 17.7% | 0.0% |
| AIC | exclude | filter_center | 33.71% | 12.2 / 10 / 1 / 53 | 58.0 / 30 / 8 / 1007 | 14.9% | 0.0% |
| AIC | exclude | filter_region | 36.76% | 14.3 / 12 / 1 / 67 | 42.8 / 26 / 8 / 508 | 18.4% | 0.0% |
| AIC | exclude | exclude_center | 10.43% | 3.4 / 3 / 1 / 9 | 311.4 / 169 / 10 / 2115 | 11.7% | 0.0% |
| AIC | exclude | exclude_region | 0.12% | 3.0 / 2 / 1 / 18 | 374.9 / 210 / 12 / 2139 | 16.1% | 0.0% |
| AIC | remove | – | 47.59% | 7.3 / 6 / 1 / 28 | 92.9 / 32 / 8 / 1599 | 3.6% | 0.0% |
| AIC | remove | truncate_center | 47.87% | 14.2 / 11 / 1 / 68 | 41.5 / 23 / 8 / 984 | 9.7% | 0.0% |
| AIC | remove | filter_center | 47.74% | 12.2 / 10 / 1 / 53 | 50.0 / 24 / 8 / 971 | 7.9% | 0.0% |
| AIC | remove | filter_region | 47.91% | 14.9 / 12 / 1 / 67 | 39.3 / 22 / 8 / 498 | 9.9% | 0.0% |
| AIC | remove | exclude_center | 47.61% | 3.4 / 3 / 1 / 9 | 207.3 / 94 / 10 / 1693 | 1.0% | 0.0% |
| AIC | remove | exclude_region | 47.60% | 3.0 / 2 / 1 / 18 | 231.7 / 100 / 12 / 1457 | 1.0% | 0.0% |
| AIC weight | – | – | 0.00% | 4.5 / 4 / 1 / 15 | 180.9 / 58 / 8 / 2591 | 40.0% | 0.0% |
| AIC weight | – | truncate_center | 9.73% | 6.3 / 5 / 0 / 21 | 56.5 / 27 / 8 / 930 | 64.8% | 4.2% |
| AIC weight | – | filter_center | 0.23% | 5.8 / 5 / 0 / 25 | 78.2 / 26 / 8 / 1081 | 64.4% | 2.1% |
| AIC weight | – | filter_region | 12.86% | 5.8 / 4 / 0 / 25 | 61.8 / 29 / 8 / 655 | 59.5% | 3.5% |
| AIC weight | – | exclude_center | 0.03% | 6.3 / 5 / 1 / 22 | 210.7 / 66 / 11 / 2115 | 4.2% | 0.0% |
| AIC weight | – | exclude_region | 3.20% | 6.9 / 5 / 1 / 31 | 175.4 / 51 / 12 / 1814 | 8.1% | 0.0% |
| AIC weight | exclude | – | 0.00% | 3.3 / 3 / 0 / 14 | 182.0 / 65 / 8 / 1666 | 55.8% | 4.9% |
| AIC weight | exclude | truncate_center | 36.87% | 4.5 / 4 / 0 / 22 | 56.8 / 34 / 8 / 512 | 45.6% | 8.3% |
| AIC weight | exclude | filter_center | 33.71% | 3.9 / 3 / 0 / 20 | 87.4 / 39 / 8 / 1007 | 42.9% | 6.9% |
| AIC weight | exclude | filter_region | 36.76% | 4.2 / 3 / 0 / 22 | 64.1 / 33 / 8 / 508 | 44.1% | 6.9% |
| AIC weight | exclude | exclude_center | 10.43% | 3.2 / 3 / 1 / 9 | 327.5 / 187 / 14 / 2115 | 12.3% | 0.0% |
| AIC weight | exclude | exclude_region | 0.12% | 2.9 / 2 / 1 / 16 | 395.1 / 238 / 13 / 2139 | 17.1% | 0.0% |
| AIC weight | remove | – | 47.59% | 3.3 / 3 / 0 / 14 | 120.3 / 46 / 8 / 1100 | 24.4% | 4.9% |
| AIC weight | remove | truncate_center | 47.87% | 4.1 / 4 / 0 / 19 | 68.8 / 33 / 8 / 984 | 33.0% | 7.6% |
| AIC weight | remove | filter_center | 47.74% | 3.9 / 3 / 0 / 20 | 77.2 / 34 / 8 / 971 | 31.9% | 6.9% |
| AIC weight | remove | filter_region | 47.91% | 4.4 / 4 / 0 / 22 | 62.6 / 32 / 8 / 498 | 33.5% | 5.6% |
| AIC weight | remove | exclude_center | 47.61% | 3.2 / 3 / 1 / 9 | 219.1 / 115 / 14 / 1693 | 1.4% | 0.0% |
| AIC weight | remove | exclude_region | 47.60% | 2.9 / 2 / 1 / 16 | 245.0 / 116 / 13 / 1457 | 1.4% | 0.0% |
| LRT | – | – | 0.00% | 3.1 / 3 / 0 / 11 | 207.1 / 90 / 11 / 2591 | 51.8% | 9.0% |
| LRT | – | truncate_center | 9.73% | 4.2 / 3 / 0 / 17 | 68.7 / 35 / 10 / 930 | 69.8% | 11.8% |
| LRT | – | filter_center | 0.23% | 3.7 / 3 / 0 / 15 | 103.7 / 38 / 10 / 1081 | 69.8% | 9.0% |
| LRT | – | filter_region | 12.86% | 3.8 / 3 / 0 / 14 | 77.5 / 40 / 10 / 655 | 64.3% | 11.8% |
| LRT | – | exclude_center | 0.03% | 6.0 / 4 / 1 / 21 | 219.6 / 70 / 12 / 2115 | 4.7% | 0.0% |
| LRT | – | exclude_region | 3.20% | 6.5 / 5 / 1 / 29 | 185.2 / 57 / 12 / 1814 | 8.6% | 0.0% |
| LRT | exclude | – | 0.00% | 2.2 / 2 / 0 / 12 | 216.8 / 83 / 12 / 1666 | 65.8% | 16.0% |
| LRT | exclude | truncate_center | 36.87% | 2.5 / 2 / 0 / 15 | 69.9 / 43 / 12 / 512 | 51.6% | 22.2% |
| LRT | exclude | filter_center | 33.71% | 2.3 / 2 / 0 / 10 | 115.1 / 54 / 12 / 1007 | 48.6% | 19.4% |
| LRT | exclude | filter_region | 36.76% | 2.5 / 2 / 0 / 12 | 82.0 / 44 / 12 / 508 | 49.2% | 20.1% |
| LRT | exclude | exclude_center | 10.43% | 3.1 / 3 / 0 / 8 | 336.6 / 193 / 14 / 2115 | 12.9% | 0.7% |
| LRT | exclude | exclude_region | 0.12% | 2.7 / 2 / 0 / 15 | 407.3 / 245 / 16 / 2139 | 18.3% | 1.4% |
| LRT | remove | – | 47.59% | 2.2 / 2 / 0 / 12 | 147.4 / 61 / 12 / 1100 | 30.8% | 16.0% |
| LRT | remove | truncate_center | 47.87% | 2.5 / 2 / 0 / 12 | 86.6 / 42 / 12 / 984 | 38.6% | 20.1% |
| LRT | remove | filter_center | 47.74% | 2.3 / 2 / 0 / 10 | 102.0 / 46 / 12 / 971 | 37.0% | 19.4% |
| LRT | remove | filter_region | 47.91% | 2.5 / 2 / 0 / 11 | 80.0 / 41 / 11 / 498 | 38.8% | 18.8% |
| LRT | remove | exclude_center | 47.61% | 3.1 / 3 / 0 / 8 | 225.5 / 121 / 14 / 1693 | 1.8% | 0.7% |
| LRT | remove | exclude_region | 47.60% | 2.7 / 2 / 0 / 15 | 255.3 / 124 / 16 / 1457 | 1.8% | 1.4% |

## Metrics explained

- **% Assigned missing**: When outputs have one row per residue (with NA for filtered): fraction of residues with NA for obs/exp in this run. When only assigned residues are emitted: fraction of residues absent from this run vs the unfiltered run (same model).
- **Const. regions**: Number of constraint regions (non-null) per transcript: mean, median, min, max.
- **Region length**: Residues per constraint region: mean, median, min, max.
- **Null frac (mean)**: Mean over transcripts of (residues in null region / total residues).
- **Tx all-null %**: Fraction of transcripts with zero constraint regions (all residues in null).

## Expectation for missing values (NA) by filtering method

**% Assigned missing** is the fraction of residues with **NA** for obs/exp in the per-residue output.

The forward algorithm builds a set of **valid_residues**: the union of all residues that appear in at least one candidate region. Only valid_residues get a region assignment (either a constraint region or the null catch-all). The final output left-joins all residues (from the OE array) against the forward results, so residues **not** in valid_residues get **NA** for region_index, obs, exp, oe, etc.

Whether a residue ends up in valid_residues depends on whether it was **hard-filtered** (physically removed from the distance matrix) or **soft-excluded** (kept in the distance matrix but marked `exclude_from_stats`).

### pLDDT filter

- **–** (none): All residues are in valid_residues. **0% NA**.
- **exclude** (`exclude_low_plddt_from_stats`): Low-pLDDT residues stay in the distance matrix and appear in candidate regions → they are **in valid_residues** and get a region assignment (constraint or null). Their obs/exp are excluded from region-level nLL/OE calculations via `excluded_residues`, but the per-residue export row is **not NA**. **0% NA** from pLDDT exclude alone.
- **remove** (`remove_low_plddt_residues`): Low-pLDDT residues are physically removed from the distance matrix → they never appear in any candidate region → **not in valid_residues** → **NA** in the export. Typically **~40–47% NA**.

### PAE filter

- **–** (none): No residues removed by PAE. **0% NA** from PAE.
- **truncate_center**: All residues **after** the first with PAE(center, neighbor) &gt; cutoff are removed from the distance matrix (sequential cutoff). Residues that are truncated from **every** center's candidate region are not in valid_residues → **NA**.
- **filter_center**: Only residues with PAE(center, neighbor) &gt; cutoff are removed from the distance matrix. Residues filtered from **every** center's candidate region are not in valid_residues → **NA**.
- **filter_region**: Residues whose maximum pairwise PAE to any residue in the region &gt; cutoff are removed from the distance matrix. Residues filtered from **every** center are not in valid_residues → **NA**.
- **exclude_center** (`exclude_on_pairwise_pae_with_center`): High-PAE residues stay in the distance matrix and appear in candidate regions → **in valid_residues**. They are marked `exclude_from_stats` (excluded from region nLL/OE), but their export row is **not NA**. Very few residues end up NA (only those that happen to not appear in any candidate region for other reasons).
- **exclude_region** (`exclude_on_pairwise_pae_in_region`): Same as exclude_center but based on pairwise PAE to any residue in the region. Residues stay in valid_residues → export row is **not NA**. Slightly more residues may be excluded from stats than exclude_center, but still very low NA.

### Combined

Runs that use both pLDDT and PAE filters can have more NA residues. For example, **remove** + **truncate_center** removes low-pLDDT residues and truncates high-PAE residues, so both sets are not in valid_residues → NA. Conversely, **exclude** + **exclude_center** keeps all residues in valid_residues (0% NA from either) but excludes them from region statistics.

The **% Assigned missing** column summarizes the fraction of residues with NA in each run.

## What to look for when choosing a method

- **More constraint regions** (higher mean/median) → more granular; **fewer** → more conservative.
- **Larger null fraction** → more residues in catch-all; smaller → more residues assigned to constraint regions.
- **High Tx all-null %** → many genes end up with no constraint regions (may be too strict).
- **% Assigned missing** → residues with NA (or absent when only assigned residues are emitted); baseline = unfiltered run for that model.
- **Region length**: very short regions may be noise; very long may be overmerged.
- **LRT** is usually most conservative (fewest regions); **AIC** most permissive; **AIC weight** tunable.

## Example genes to characterize

Genes below are candidates for full characterization: one stable across methods, and several outliers for different reasons. For each gene we show summary stats, **constraint regions as residue ranges** (for contrasting runs), and **what differs** across methods.

### Stable (low variation in region count across methods)

#### **Q5SQ80** (stable)

Low variation: std(mean regions across runs) = 0.76; rarely all-null.

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **1** to **3** across runs.
- Mean constraint regions per transcript ranges from **1.0** to **3.0**.
- Null fraction (residues in catch-all region) ranges from **13.1%** to **44.1%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, PAE exclude_center** — 1 constraint region(s) total, 32.7% in null.
  - **Transcript ENST00000377601**: Region 1: residues 0–663 (554 residues); Null (catch-all): residues 370–822 (269 residues)

- **Run: AIC (no pLDDT/PAE)** — 3 constraint region(s) total, 13.6% in null.
  - **Transcript ENST00000377601**: Region 1: residues 46–822 (275 residues); Region 2: residues 0–600 (372 residues); Region 3: residues 28–618 (64 residues); Null (catch-all): residues 71–615 (112 residues)

- **Run: AIC weight, PAE filter_region** — 1 constraint region(s) total, 23.5% in null.
  - **Transcript ENST00000377601**: Region 1: residues 38–262 (209 residues); Null (catch-all): residues 3–822 (193 residues)

**What is different across methods:**

- **Region count:** **AIC weight, PAE exclude_center** has the fewest constraint regions (1); **AIC (no pLDDT/PAE)** has the most (3). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC weight, pLDDT remove, PAE filter_region** (13.1%); highest in **AIC weight, PAE filter_center** (44.1%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.


### Outlier: region count varies most across methods

#### **Q13813** (region count variance)

Region count varies a lot: std(mean regions) = 20.94, median = 14.0.

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **6** to **82** across runs.
- Mean constraint regions per transcript ranges from **6.0** to **82.0**.
- Null fraction (residues in catch-all region) ranges from **0.4%** to **87.3%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, pLDDT exclude, PAE exclude_center** — 6 constraint region(s) total, 3.6% in null.
  - **Transcript ENST00000372731**: Region 1: residues 179–2304 (1099 residues); Region 2: residues 475–2471 (761 residues); Region 3: residues 596–2456 (208 residues); Region 4: residues 174–1522 (66 residues); Region 5: residues 48–245 (174 residues); Region 6: residues 837–1186 (74 residues); Null (catch-all): residues 0–1524 (90 residues)

- **Run: AIC, PAE filter_region** — 82 constraint region(s) total, 21.1% in null.
  - **Transcript ENST00000372731**: Region 1: residues 2216–2316 (69 residues); Region 2: residues 1686–1767 (15 residues); Region 3: residues 1455–1540 (33 residues); Region 4: residues 1660–1738 (21 residues); Region 5: residues 175–254 (27 residues); Region 6: residues 1234–1315 (25 residues); Region 7: residues 1578–1657 (21 residues); Region 8: residues 1791–1869 (25 residues); Region 9: residues 1985–2067 (15 residues); Region 10: residues 150–232 (22 residues); Region 11: residues 1246–1323 (15 residues); Region 12: residues 2318–2401 (17 residues); Region 13: residues 1183–1197 (13 residues); Region 14: residues 2108–2192 (22 residues); Region 15: residues 1881–1957 (18 residues); Region 16: residues 71–145 (19 residues); Region 17: residues 599–673 (18 residues); Region 18: residues 281–370 (23 residues); Region 19: residues 261–345 (21 residues); Region 20: residues 920–1161 (28 residues); Region 21: residues 1105–1232 (49 residues); Region 22: residues 386–466 (24 residues); Region 23: residues 1548–1628 (19 residues); Region 24: residues 810–889 (20 residues); Region 25: residues 2205–2282 (14 residues); Region 26: residues 727–738 (12 residues); Region 27: residues 1562–1646 (36 residues); Region 28: residues 899–1076 (25 residues); Region 29: residues 708–871 (42 residues); Region 30: residues 2328–2391 (16 residues); Region 31: residues 1772–1857 (27 residues); Region 32: residues 295–309 (14 residues); Region 33: residues 493–572 (20 residues); Region 34: residues 1365–1444 (18 residues); Region 35: residues 2092–2174 (32 residues); Region 36: residues 1273–1287 (12 residues); Region 37: residues 600–750 (28 residues); Region 38: residues 1888–1976 (34 residues); Region 39: residues 1789–1946 (26 residues); Region 40: residues 2–21 (16 residues); Region 41: residues 178–334 (20 residues); Region 42: residues 926–1083 (18 residues); Region 43: residues 1673–1752 (48 residues); Region 44: residues 809–1064 (20 residues); Region 45: residues 1260–1411 (32 residues); Region 46: residues 570–668 (76 residues); Region 47: residues 971–1015 (16 residues); Region 48: residues 375–455 (17 residues); Region 49: residues 2004–2089 (14 residues); Region 50: residues 212–235 (16 residues); Region 51: residues 1347–1433 (34 residues); Region 52: residues 1346–1442 (38 residues); Region 53: residues 1683–1837 (29 residues); Region 54: residues 1452–1529 (12 residues); Region 55: residues 1467–1547 (19 residues); Region 56: residues 368–447 (27 residues); Region 57: residues 1980–2064 (18 residues); Region 58: residues 478–554 (11 residues); Region 59: residues 1099–1214 (31 residues); Region 60: residues 2330–2383 (37 residues); Region 61: residues 55–149 (43 residues); Region 62: residues 34–119 (29 residues); Region 63: residues 457–547 (28 residues); Region 64: residues 1120–1316 (13 residues); Region 65: residues 2122–2202 (16 residues); Region 66: residues 2214–2285 (13 residues); Region 67: residues 689–768 (13 residues); Region 68: residues 1035–1047 (13 residues); Region 69: residues 701–779 (15 residues); Region 70: residues 966–1026 (36 residues); Region 71: residues 2412–2459 (22 residues); Region 72: residues 2402–2471 (16 residues); Region 73: residues 163–243 (18 residues); Region 74: residues 1998–2079 (25 residues); Region 75: residues 828–875 (16 residues); Region 76: residues 1883–2055 (26 residues); Region 77: residues 276–352 (16 residues); Region 78: residues 2091–2183 (28 residues); Region 79: residues 1770–1848 (17 residues); Region 80: residues 1859–1951 (15 residues); Region 81: residues 1474–1631 (27 residues); Region 82: residues 2445–2460 (13 residues); Null (catch-all): residues 11–2467 (522 residues)

**What is different across methods:**

- **Region count:** **AIC weight, pLDDT exclude, PAE exclude_center** has the fewest constraint regions (6); **AIC, PAE filter_region** has the most (82). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC weight, PAE exclude_center** (0.4%); highest in **LRT, PAE filter_region** (87.3%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.

#### **P02751** (region count variance)

Region count varies a lot: std(mean regions) = 18.79, median = 7.0.

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **1** to **81** across runs.
- Mean constraint regions per transcript ranges from **1.0** to **81.0**.
- Null fraction (residues in catch-all region) ranges from **0.3%** to **95.9%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: LRT, pLDDT exclude** — 1 constraint region(s) total, 80.1% in null.
  - **Transcript ENST00000354785**: Region 1: residues 120–2417 (492 residues); Null (catch-all): residues 0–2476 (1985 residues)

- **Run: AIC, PAE filter_center** — 81 constraint region(s) total, 19.6% in null.
  - **Transcript ENST00000354785**: Region 1: residues 730–888 (36 residues); Region 2: residues 407–461 (26 residues); Region 3: residues 349–398 (25 residues); Region 4: residues 229–267 (15 residues); Region 5: residues 527–556 (17 residues); Region 6: residues 118–174 (16 residues); Region 7: residues 825–879 (27 residues); Region 8: residues 1297–1349 (12 residues); Region 9: residues 1727–1892 (48 residues); Region 10: residues 559–598 (17 residues); Region 11: residues 1903–1986 (32 residues); Region 12: residues 1454–1518 (24 residues); Region 13: residues 1460–1621 (24 residues); Region 14: residues 725–802 (29 residues); Region 15: residues 2291–2335 (27 residues); Region 16: residues 1093–1152 (29 residues); Region 17: residues 2352–2419 (38 residues); Region 18: residues 909–980 (33 residues); Region 19: residues 185–226 (17 residues); Region 20: residues 307–341 (13 residues); Region 21: residues 1999–2079 (16 residues); Region 22: residues 1001–1083 (36 residues); Region 23: residues 1643–1802 (17 residues); Region 24: residues 350–432 (20 residues); Region 25: residues 469–507 (14 residues); Region 26: residues 53–86 (17 residues); Region 27: residues 1099–1257 (42 residues); Region 28: residues 641–681 (29 residues); Region 29: residues 917–994 (18 residues); Region 30: residues 911–1088 (83 residues); Region 31: residues 1914–1992 (16 residues); Region 32: residues 142–219 (39 residues); Region 33: residues 1634–1686 (16 residues); Region 34: residues 96–133 (19 residues); Region 35: residues 1185–1345 (12 residues); Region 36: residues 2121–2134 (14 residues); Region 37: residues 1847–1899 (13 residues); Region 38: residues 2213–2264 (19 residues); Region 39: residues 903–988 (15 residues); Region 40: residues 1366–1416 (13 residues); Region 41: residues 114–218 (19 residues); Region 42: residues 2005–2081 (13 residues); Region 43: residues 2056–2095 (18 residues); Region 44: residues 1815–1893 (14 residues); Region 45: residues 1816–1987 (72 residues); Region 46: residues 574–611 (22 residues); Region 47: residues 728–894 (24 residues); Region 48: residues 610–694 (30 residues); Region 49: residues 2428–2444 (17 residues); Region 50: residues 1725–1794 (15 residues); Region 51: residues 1660–1705 (22 residues); Region 52: residues 1540–1626 (49 residues); Region 53: residues 1421–1523 (13 residues); Region 54: residues 1481–1533 (16 residues); Region 55: residues 1420–1538 (36 residues); Region 56: residues 207–270 (12 residues); Region 57: residues 2165–2182 (18 residues); Region 58: residues 2289–2349 (27 residues); Region 59: residues 1996–2075 (40 residues); Region 60: residues 1965–2076 (14 residues); Region 61: residues 1911–2031 (27 residues); Region 62: residues 622–714 (20 residues); Region 63: residues 13–27 (13 residues); Region 64: residues 2224–2261 (23 residues); Region 65: residues 346–481 (18 residues); Region 66: residues 520–597 (25 residues); Region 67: residues 2190–2202 (13 residues); Region 68: residues 2096–2112 (17 residues); Region 69: residues 2136–2150 (13 residues); Region 70: residues 1188–1356 (82 residues); Region 71: residues 1364–1525 (24 residues); Region 72: residues 1090–1259 (40 residues); Region 73: residues 2343–2426 (32 residues); Region 74: residues 271–283 (12 residues); Region 75: residues 2449–2463 (15 residues); Region 76: residues 315–455 (17 residues); Region 77: residues 33–45 (13 residues); Region 78: residues 51–93 (19 residues); Region 79: residues 722–804 (33 residues); Region 80: residues 484–515 (19 residues); Region 81: residues 1549–1798 (53 residues); Null (catch-all): residues 0–2476 (485 residues)

- **Run: AIC, pLDDT exclude, PAE truncate_center** — 43 constraint region(s) total, 38.3% in null.
  - **Transcript ENST00000354785**: Region 1: residues 407–462 (23 residues); Region 2: residues 730–888 (36 residues); Region 3: residues 349–400 (29 residues); Region 4: residues 200–269 (47 residues); Region 5: residues 842–879 (23 residues); Region 6: residues 1297–1349 (12 residues); Region 7: residues 112–175 (27 residues); Region 8: residues 1727–1892 (48 residues); Region 9: residues 1012–1083 (24 residues); Region 10: residues 1643–1802 (17 residues); Region 11: residues 1093–1152 (29 residues); Region 12: residues 1454–1518 (24 residues); Region 13: residues 1905–1982 (14 residues); Region 14: residues 914–994 (15 residues); Region 15: residues 605–696 (57 residues); Region 16: residues 905–987 (32 residues); Region 17: residues 1999–2081 (24 residues); Region 18: residues 1460–1623 (32 residues); Region 19: residues 1099–1257 (43 residues); Region 20: residues 739–803 (24 residues); Region 21: residues 1914–1992 (16 residues); Region 22: residues 2341–2418 (54 residues); Region 23: residues 307–448 (48 residues); Region 24: residues 1634–1686 (16 residues); Region 25: residues 1366–1416 (13 residues); Region 26: residues 51–93 (37 residues); Region 27: residues 1569–1607 (27 residues); Region 28: residues 1815–1893 (14 residues); Region 29: residues 2224–2278 (21 residues); Region 30: residues 1421–1527 (16 residues); Region 31: residues 1660–1705 (22 residues); Region 32: residues 1279–1438 (15 residues); Region 33: residues 998–1078 (17 residues); Region 34: residues 1846–1900 (20 residues); Region 35: residues 1481–1533 (16 residues); Region 36: residues 1450–1535 (25 residues); Region 37: residues 96–145 (29 residues); Region 38: residues 1996–2077 (27 residues); Region 39: residues 531–598 (48 residues); Region 40: residues 812–894 (13 residues); Region 41: residues 1750–1774 (15 residues); Region 42: residues 1545–1594 (22 residues); Region 43: residues 1543–1714 (26 residues); Null (catch-all): residues 43–2374 (949 residues)

**What is different across methods:**

- **Region count:** **LRT, pLDDT exclude** has the fewest constraint regions (1); **AIC, PAE filter_center** has the most (81). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC, pLDDT remove** (0.3%); highest in **LRT, PAE filter_center** (95.9%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.


### Outlier: often all-null (no constraint regions) across runs

#### **Q5D862** (often all-null)

In 44% of runs this gene has no constraint regions (all residues in null).

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **0** to **100** across runs.
- Mean constraint regions per transcript ranges from **0.0** to **100.0**.
- Null fraction (residues in catch-all region) ranges from **0.1%** to **100.0%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, PAE filter_center** — 0 constraint region(s) total, 100.0% in null.
  - **Transcript ENST00000388718**: Null (catch-all): residues 0–2390 (2391 residues)

- **Run: AIC, PAE filter_center** — 100 constraint region(s) total, 31.4% in null.
  - **Transcript ENST00000388718**: Region 1: residues 881–897 (17 residues); Region 2: residues 1846–1867 (20 residues); Region 3: residues 2024–2042 (18 residues); Region 4: residues 1940–1956 (17 residues); Region 5: residues 2110–2130 (20 residues); Region 6: residues 858–875 (15 residues); Region 7: residues 26–69 (13 residues); Region 8: residues 996–1013 (18 residues); Region 9: residues 1570–1583 (14 residues); Region 10: residues 299–312 (14 residues); Region 11: residues 358–377 (20 residues); Region 12: residues 1468–1481 (14 residues); Region 13: residues 615–633 (19 residues); Region 14: residues 1714–1729 (15 residues); Region 15: residues 1598–1610 (13 residues); Region 16: residues 1828–1842 (13 residues); Region 17: residues 2253–2267 (15 residues); Region 18: residues 1686–1700 (15 residues); Region 19: residues 1302–1315 (14 residues); Region 20: residues 2170–2184 (15 residues); Region 21: residues 777–791 (15 residues); Region 22: residues 1084–1102 (19 residues); Region 23: residues 85–102 (18 residues); Region 24: residues 2000–2017 (14 residues); Region 25: residues 1987–1999 (13 residues); Region 26: residues 727–746 (20 residues); Region 27: residues 1493–1506 (14 residues); Region 28: residues 1366–1383 (17 residues); Region 29: residues 1149–1162 (14 residues); Region 30: residues 2150–2165 (15 residues); Region 31: residues 2336–2356 (21 residues); Region 32: residues 477–490 (14 residues); Region 33: residues 573–586 (14 residues); Region 34: residues 1616–1628 (13 residues); Region 35: residues 700–716 (15 residues); Region 36: residues 1783–1796 (14 residues); Region 37: residues 1910–1926 (17 residues); Region 38: residues 1451–1464 (14 residues); Region 39: residues 1546–1561 (16 residues); Region 40: residues 2095–2109 (15 residues); Region 41: residues 1029–1047 (19 residues); Region 42: residues 913–929 (17 residues); Region 43: residues 751–768 (18 residues); Region 44: residues 226–239 (14 residues); Region 45: residues 389–401 (13 residues); Region 46: residues 1414–1430 (15 residues); Region 47: residues 959–974 (16 residues); Region 48: residues 833–851 (19 residues); Region 49: residues 1266–1279 (14 residues); Region 50: residues 548–562 (15 residues); Region 51: residues 2045–2057 (13 residues); Region 52: residues 1969–1982 (14 residues); Region 53: residues 803–820 (18 residues); Region 54: residues 2318–2332 (15 residues); Region 55: residues 2133–2149 (17 residues); Region 56: residues 1124–1143 (20 residues); Region 57: residues 279–294 (15 residues); Region 58: residues 931–945 (15 residues); Region 59: residues 1767–1780 (14 residues); Region 60: residues 186–203 (18 residues); Region 61: residues 441–454 (14 residues); Region 62: residues 672–686 (15 residues); Region 63: residues 133–148 (16 residues); Region 64: residues 1221–1237 (17 residues); Region 65: residues 461–475 (15 residues); Region 66: residues 2302–2314 (13 residues); Region 67: residues 1166–1181 (16 residues); Region 68: residues 1244–1258 (15 residues); Region 69: residues 2066–2080 (15 residues); Region 70: residues 337–357 (20 residues); Region 71: residues 1651–1668 (18 residues); Region 72: residues 406–431 (26 residues); Region 73: residues 2201–2220 (20 residues); Region 74: residues 0–84 (39 residues); Region 75: residues 1068–1082 (15 residues); Region 76: residues 1202–1218 (17 residues); Region 77: residues 530–544 (15 residues); Region 78: residues 1335–1350 (14 residues); Region 79: residues 246–264 (19 residues); Region 80: residues 1877–1900 (24 residues); Region 81: residues 2228–2240 (13 residues); Region 82: residues 158–180 (21 residues); Region 83: residues 591–607 (17 residues); Region 84: residues 108–123 (16 residues); Region 85: residues 1014–1028 (15 residues); Region 86: residues 649–664 (16 residues); Region 87: residues 498–518 (21 residues); Region 88: residues 1533–1545 (13 residues); Region 89: residues 1106–1119 (14 residues); Region 90: residues 1050–1065 (16 residues); Region 91: residues 1807–1823 (17 residues); Region 92: residues 2372–2387 (16 residues); Region 93: residues 1387–1410 (24 residues); Region 94: residues 1434–1450 (17 residues); Region 95: residues 634–647 (14 residues); Region 96: residues 1509–1523 (15 residues); Region 97: residues 1316–1331 (16 residues); Region 98: residues 1743–1756 (14 residues); Region 99: residues 1632–1648 (17 residues); Region 100: residues 2274–2289 (15 residues); Null (catch-all): residues 19–2390 (751 residues)

- **Run: LRT, pLDDT exclude, PAE filter_center** — 0 constraint region(s) total, 4.9% in null.
  - **Transcript ENST00000388718**: Null (catch-all): residues 0–139 (117 residues)

**What is different across methods:**

- **Region count:** **AIC weight, PAE filter_center** has the fewest constraint regions (0); **AIC, PAE filter_center** has the most (100). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC, pLDDT remove, PAE filter_region** (0.1%); highest in **AIC weight, PAE filter_center** (100.0%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.
- **All-null runs:** In **24** run(s) this gene has *no* constraint regions (all residues in null). Example: AIC weight, PAE filter_center

#### **Q6EMK4** (often all-null)

In 43% of runs this gene has no constraint regions (all residues in null).

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **0** to **18** across runs.
- Mean constraint regions per transcript ranges from **0.0** to **18.0**.
- Null fraction (residues in catch-all region) ranges from **0.0%** to **100.0%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, PAE filter_center** — 0 constraint region(s) total, 100.0% in null.
  - **Transcript ENST00000304735**: Null (catch-all): residues 0–672 (673 residues)

- **Run: AIC, PAE filter_center** — 18 constraint region(s) total, 14.4% in null.
  - **Transcript ENST00000304735**: Region 1: residues 21–174 (119 residues); Region 2: residues 331–519 (14 residues); Region 3: residues 582–603 (22 residues); Region 4: residues 390–403 (14 residues); Region 5: residues 498–555 (19 residues); Region 6: residues 330–548 (37 residues); Region 7: residues 344–360 (14 residues); Region 8: residues 413–440 (19 residues); Region 9: residues 64–561 (167 residues); Region 10: residues 636–650 (15 residues); Region 11: residues 614–627 (14 residues); Region 12: residues 658–670 (13 residues); Region 13: residues 454–550 (14 residues); Region 14: residues 404–441 (17 residues); Region 15: residues 323–567 (14 residues); Region 16: residues 363–386 (24 residues); Region 17: residues 4–18 (15 residues); Region 18: residues 41–523 (25 residues); Null (catch-all): residues 0–672 (97 residues)

- **Run: AIC, pLDDT remove, PAE filter_region** — 9 constraint region(s) total, 25.1% in null.
  - **Transcript ENST00000304735**: Region 1: residues 54–151 (34 residues); Region 2: residues 44–94 (15 residues); Region 3: residues 331–519 (14 residues); Region 4: residues 579–603 (22 residues); Region 5: residues 465–552 (22 residues); Region 6: residues 152–322 (114 residues); Region 7: residues 413–440 (19 residues); Region 8: residues 26–157 (37 residues); Region 9: residues 498–555 (15 residues); Null (catch-all): residues 22–585 (169 residues)

**What is different across methods:**

- **Region count:** **AIC weight, PAE filter_center** has the fewest constraint regions (0); **AIC, PAE filter_center** has the most (18). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC weight, pLDDT remove, PAE exclude_center** (0.0%); highest in **AIC weight, PAE filter_center** (100.0%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.
- **All-null runs:** In **23** run(s) this gene has *no* constraint regions (all residues in null). Example: AIC weight, PAE filter_center


### Outlier: null fraction varies most across runs

#### **Q96RG2** (null fraction range)

Null fraction across runs has range 1.00 (full spread).

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **0** to **40** across runs.
- Mean constraint regions per transcript ranges from **0.0** to **40.0**.
- Null fraction (residues in catch-all region) ranges from **0.0%** to **100.0%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: LRT, PAE filter_center** — 0 constraint region(s) total, 100.0% in null.
  - **Transcript ENST00000234040**: Null (catch-all): residues 0–1322 (1323 residues)

- **Run: AIC, PAE filter_center** — 40 constraint region(s) total, 26.4% in null.
  - **Transcript ENST00000234040**: Region 1: residues 1004–1138 (20 residues); Region 2: residues 279–369 (17 residues); Region 3: residues 1276–1296 (21 residues); Region 4: residues 350–924 (16 residues); Region 5: residues 1007–1244 (109 residues); Region 6: residues 75–363 (64 residues); Region 7: residues 549–568 (20 residues); Region 8: residues 1099–1264 (14 residues); Region 9: residues 788–802 (15 residues); Region 10: residues 904–941 (17 residues); Region 11: residues 125–928 (86 residues); Region 12: residues 237–338 (24 residues); Region 13: residues 725–741 (17 residues); Region 14: residues 1303–1321 (19 residues); Region 15: residues 93–101 (9 residues); Region 16: residues 485–504 (20 residues); Region 17: residues 858–875 (16 residues); Region 18: residues 753–767 (15 residues); Region 19: residues 976–1075 (27 residues); Region 20: residues 52–73 (22 residues); Region 21: residues 650–663 (14 residues); Region 22: residues 956–972 (17 residues); Region 23: residues 694–714 (21 residues); Region 24: residues 154–210 (15 residues); Region 25: residues 604–622 (19 residues); Region 26: residues 1108–1250 (17 residues); Region 27: residues 812–825 (14 residues); Region 28: residues 1093–1274 (17 residues); Region 29: residues 463–483 (21 residues); Region 30: residues 242–335 (33 residues); Region 31: residues 108–123 (16 residues); Region 32: residues 74–225 (34 residues); Region 33: residues 672–682 (11 residues); Region 34: residues 18–40 (23 residues); Region 35: residues 839–852 (14 residues); Region 36: residues 402–413 (12 residues); Region 37: residues 1009–1253 (60 residues); Region 38: residues 435–454 (20 residues); Region 39: residues 524–537 (14 residues); Region 40: residues 505–518 (14 residues); Null (catch-all): residues 0–1322 (349 residues)

- **Run: AIC, pLDDT exclude, PAE truncate_center** — 13 constraint region(s) total, 23.2% in null.
  - **Transcript ENST00000234040**: Region 1: residues 1004–1138 (20 residues); Region 2: residues 345–930 (71 residues); Region 3: residues 199–366 (49 residues); Region 4: residues 252–305 (16 residues); Region 5: residues 1099–1264 (14 residues); Region 6: residues 1115–1184 (20 residues); Region 7: residues 132–230 (12 residues); Region 8: residues 1166–1233 (20 residues); Region 9: residues 154–210 (14 residues); Region 10: residues 1105–1252 (21 residues); Region 11: residues 979–1077 (29 residues); Region 12: residues 1037–1181 (28 residues); Region 13: residues 1089–1266 (13 residues); Null (catch-all): residues 74–1273 (307 residues)

**What is different across methods:**

- **Region count:** **LRT, PAE filter_center** has the fewest constraint regions (0); **AIC, PAE filter_center** has the most (40). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC weight, pLDDT remove, PAE exclude_region** (0.0%); highest in **LRT, PAE filter_center** (100.0%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.
- **All-null runs:** In **9** run(s) this gene has *no* constraint regions (all residues in null). Example: LRT, PAE filter_center

#### **Q8IWN7** (null fraction range)

Null fraction across runs has range 1.00 (full spread).

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **0** to **79** across runs.
- Mean constraint regions per transcript ranges from **0.0** to **79.0**.
- Null fraction (residues in catch-all region) ranges from **0.0%** to **100.0%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, PAE filter_region** — 0 constraint region(s) total, 63.1% in null.
  - **Transcript ENST00000382483**: Null (catch-all): residues 1–2399 (1515 residues)

- **Run: AIC, PAE filter_center** — 79 constraint region(s) total, 45.0% in null.
  - **Transcript ENST00000382483**: Region 1: residues 1324–1345 (22 residues); Region 2: residues 1346–1359 (14 residues); Region 3: residues 1934–1951 (18 residues); Region 4: residues 2052–2072 (21 residues); Region 5: residues 1393–1408 (16 residues); Region 6: residues 1458–1469 (12 residues); Region 7: residues 847–858 (12 residues); Region 8: residues 1522–1567 (32 residues); Region 9: residues 1855–1867 (13 residues); Region 10: residues 685–695 (11 residues); Region 11: residues 2158–2176 (16 residues); Region 12: residues 1477–1490 (14 residues); Region 13: residues 2273–2294 (22 residues); Region 14: residues 1436–1451 (16 residues); Region 15: residues 456–470 (12 residues); Region 16: residues 1305–1323 (19 residues); Region 17: residues 1870–1888 (18 residues); Region 18: residues 1088–1159 (27 residues); Region 19: residues 2255–2272 (18 residues); Region 20: residues 1195–1207 (13 residues); Region 21: residues 2088–2106 (19 residues); Region 22: residues 800–812 (13 residues); Region 23: residues 550–564 (14 residues); Region 24: residues 1599–1609 (11 residues); Region 25: residues 2011–2034 (23 residues); Region 26: residues 389–400 (12 residues); Region 27: residues 1071–1145 (13 residues); Region 28: residues 1276–1289 (14 residues); Region 29: residues 1372–1387 (16 residues); Region 30: residues 1615–1634 (20 residues); Region 31: residues 1961–1975 (15 residues); Region 32: residues 1901–1914 (14 residues); Region 33: residues 1809–1821 (13 residues); Region 34: residues 350–368 (19 residues); Region 35: residues 155–221 (45 residues); Region 36: residues 1219–1233 (15 residues); Region 37: residues 732–748 (17 residues); Region 38: residues 2175–2193 (18 residues); Region 39: residues 1118–1146 (16 residues); Region 40: residues 2382–2393 (12 residues); Region 41: residues 312–322 (10 residues); Region 42: residues 703–716 (14 residues); Region 43: residues 2337–2351 (15 residues); Region 44: residues 30–107 (44 residues); Region 45: residues 567–586 (20 residues); Region 46: residues 330–348 (18 residues); Region 47: residues 1769–1783 (15 residues); Region 48: residues 883–894 (12 residues); Region 49: residues 781–798 (18 residues); Region 50: residues 590–607 (18 residues); Region 51: residues 824–838 (15 residues); Region 52: residues 241–257 (17 residues); Region 53: residues 1979–1993 (15 residues); Region 54: residues 1504–1585 (13 residues); Region 55: residues 428–445 (18 residues); Region 56: residues 646–657 (12 residues); Region 57: residues 1829–1839 (11 residues); Region 58: residues 2109–2125 (17 residues); Region 59: residues 983–995 (13 residues); Region 60: residues 1713–1733 (19 residues); Region 61: residues 510–521 (12 residues); Region 62: residues 260–278 (19 residues); Region 63: residues 2320–2332 (13 residues); Region 64: residues 628–643 (16 residues); Region 65: residues 666–681 (16 residues); Region 66: residues 1581–1595 (13 residues); Region 67: residues 2214–2234 (21 residues); Region 68: residues 2354–2371 (18 residues); Region 69: residues 1698–1714 (16 residues); Region 70: residues 2033–2049 (16 residues); Region 71: residues 610–623 (14 residues); Region 72: residues 1787–1805 (19 residues); Region 73: residues 1739–1750 (12 residues); Region 74: residues 916–927 (12 residues); Region 75: residues 294–306 (13 residues); Region 76: residues 765–780 (14 residues); Region 77: residues 1237–1252 (16 residues); Region 78: residues 144–226 (22 residues); Region 79: residues 932–950 (19 residues); Null (catch-all): residues 0–2399 (1080 residues)

- **Run: LRT, pLDDT exclude, PAE truncate_center** — 0 constraint region(s) total, 14.6% in null.
  - **Transcript ENST00000382483**: Null (catch-all): residues 23–1591 (351 residues)

**What is different across methods:**

- **Region count:** **AIC weight, PAE filter_region** has the fewest constraint regions (0); **AIC, PAE filter_center** has the most (79). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC weight, pLDDT remove, PAE exclude_region** (0.0%); highest in **AIC weight, pLDDT exclude** (100.0%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.
- **All-null runs:** In **21** run(s) this gene has *no* constraint regions (all residues in null). Example: AIC weight, PAE filter_region


### User-specified genes

#### **P13637** (user-specified)

User-requested gene for characterization.

**Summary across 54 runs:**

- Number of constraint regions (total over transcripts) ranges from **2** to **17** across runs.
- Mean constraint regions per transcript ranges from **2.0** to **17.0**.
- Null fraction (residues in catch-all region) ranges from **0.0%** to **34.6%**.

**Constraint regions (residue ranges) — contrasting runs:**

Each *constraint region* is a contiguous block of residues with a shared constraint model; the *null* region is the catch-all for residues not assigned to any constraint region.

- **Run: AIC weight, PAE exclude_region** — 2 constraint region(s) total, 1.7% in null.
  - **Transcript ENST00000648268**: Region 1: residues 0–1012 (974 residues); Region 2: residues 485–559 (22 residues); Null (catch-all): residues 864–896 (17 residues)

- **Run: AIC, pLDDT exclude, PAE truncate_center** — 17 constraint region(s) total, 2.4% in null.
  - **Transcript ENST00000648268**: Region 1: residues 21–1012 (367 residues); Region 2: residues 95–880 (71 residues); Region 3: residues 210–615 (105 residues); Region 4: residues 383–574 (112 residues); Region 5: residues 202–681 (14 residues); Region 6: residues 46–247 (45 residues); Region 7: residues 878–900 (15 residues); Region 8: residues 774–912 (12 residues); Region 9: residues 201–655 (17 residues); Region 10: residues 14–1003 (33 residues); Region 11: residues 10–257 (30 residues); Region 12: residues 109–881 (11 residues); Region 13: residues 798–981 (23 residues); Region 14: residues 222–1005 (17 residues); Region 15: residues 113–979 (57 residues); Region 16: residues 37–242 (41 residues); Region 17: residues 846–1006 (11 residues); Null (catch-all): residues 8–993 (24 residues)

- **Run: AIC, PAE truncate_center** — 10 constraint region(s) total, 14.3% in null.
  - **Transcript ENST00000648268**: Region 1: residues 21–1012 (367 residues); Region 2: residues 95–880 (71 residues); Region 3: residues 210–615 (105 residues); Region 4: residues 24–258 (97 residues); Region 5: residues 109–980 (109 residues); Region 6: residues 526–681 (19 residues); Region 7: residues 46–167 (19 residues); Region 8: residues 396–455 (14 residues); Region 9: residues 10–666 (24 residues); Region 10: residues 221–673 (43 residues); Null (catch-all): residues 0–1006 (145 residues)

**What is different across methods:**

- **Region count:** **AIC weight, PAE exclude_region** has the fewest constraint regions (2); **AIC, pLDDT exclude, PAE truncate_center** has the most (17). Different filtering or model comparison (AIC vs AIC weight vs LRT) can split or merge regions, or move residues into the null.
- **Null fraction:** Lowest in **AIC, PAE exclude_center** (0.0%); highest in **LRT, pLDDT exclude, PAE truncate_center** (34.6%). Runs that exclude more residues from stats (e.g. PAE exclude) often put more residues into the null catch-all.


## Plots

Violin plots show the full distribution of each metric across transcripts or regions (one violin per run).

### Distribution: constraint regions per transcript

![Distribution: constraint regions per transcript](report_plots/dist_regions_per_tx.png)

### Distribution: region length (residues per region)

![Distribution: region length (residues per region)](report_plots/dist_region_length.png)

### Distribution: null fraction per transcript

![Distribution: null fraction per transcript](report_plots/dist_null_frac_per_tx.png)

Per-gene heatmaps show each metric (mean constraint regions or null fraction) for every gene (rows) and run (columns).

### Per-gene heatmap: mean constraint regions per transcript (genes × runs, hierarchically clustered)

![Per-gene heatmap: mean constraint regions per transcript (genes × runs, hierarchically clustered)](report_plots/heatmap_regions_per_gene.png)

### Per-gene heatmap: mean null fraction % (genes × runs, hierarchically clustered)

![Per-gene heatmap: mean null fraction % (genes × runs, hierarchically clustered)](report_plots/heatmap_null_frac_per_gene.png)
