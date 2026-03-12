## Proemis3D Forward Method Recommendation (based on comparison report)

Based on the observed patterns in the generated report (assuming default metrics and a goal of balanced sensitivity and specificity), here's a recommendation for the Proemis3D forward method combination:

---

### Primary Choice: AIC, pLDDT exclude, PAE filter_region

This combination offers a good balance of granularity and data retention based on the report's statistics:

*   **Constraint Regions**: This combination yields a mean of **12.5 regions** per transcript with 0% transcripts showing no constraint regions. This indicates good granularity without excessive fragmentation.
*   **Null Fraction**: The mean null fraction is **13.8%**, which is relatively low compared to other filtering options, suggesting that a good portion of residues are being assigned to constraint regions rather than being put into the catch-all.
*   **Filtered Out %**: Approximately **35.8%** of residues are filtered out compared to the unfiltered baseline. This indicates a significant removal of low-confidence or unreliably modeled residues, improving the quality of the remaining data.
*   **Method Rationale**:
    *   **AIC**: Tends to be more permissive, resulting in more granular regions, which aligns with aiming for a more detailed view of constraint.
    *   **pLDDT exclude**: This method keeps residues in the analysis but excludes those with low pLDDT scores from the statistical calculations. This is crucial for avoiding data loss (which `remove` or `truncate` would cause) while mitigating the influence of low-confidence structural predictions on region definitions.
    *   **PAE filter_region**: This method filters out residues based on their pairwise predicted aligned error (PAE) within a region. It's more stringent than `filter_center` in this report, resulting in more refined regions by removing unreliable structural contacts within potential constraint areas.

### Alternatives

*   **AIC weight, pLDDT exclude, PAE filter_region** (10.5 regions, 19.9% null): If a slightly more conservative (fewer regions) yet still granular model is desired, using AIC weight can provide a more stable model selection than raw AIC.
*   **AIC, no pLDDT/PAE** (9.2 regions, 6.9% null): If you prefer to analyze all residues without any confidence-based filtering, this provides a baseline with high coverage, although potentially including more noise from low-confidence regions.

### Combinations to Generally Avoid

*   **PAE exclude_center / exclude_region**: These methods tend to over-downweight, resulting in approximately **1 constraint region per transcript** and very high null fractions (~80–96%). This essentially collapses most constraint into a single catch-all region, losing granularity.
*   **LRT + pLDDT truncate**: This combination consistently showed a high **Tx all-null % (~29%)**, meaning a large proportion of transcripts ended up with no constraint regions at all, suggesting it is overly conservative for this dataset.

---

This recommendation prioritizes a balance between generating a sufficiently granular set of constraint regions and ensuring the quality of the underlying structural data, while avoiding methods that lead to excessive data loss or overly conservative results.