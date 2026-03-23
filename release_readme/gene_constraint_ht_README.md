# Constraint Field Descriptions

Descriptions of columns in the gnomAD v4.1.1 gene constraint metrics tsv. Descriptions also apply to rows in Hail Tables, where a "." in the field name indicates a struct. All constraint metrics were calculated using the gnomAD v4.1.1 exomes.

## Globals

Pipeline parameters and metadata stored as globals on the Hail Table.

- `version`: gnomAD constraint release version

### `calculate_mu_params`

Parameters used for the genome-based mutation rate calculation.

- `calculate_mu_params.ac_cutoff`: Maximum allele count cutoff for variants included in the mutation rate calculation
- `calculate_mu_params.min_cov`: Minimum mean genome coverage for a site to be included
- `calculate_mu_params.max_cov`: Maximum mean genome coverage for a site to be included
- `calculate_mu_params.gerp_lower_cutoff`: Minimum GERP score for a site to be included (default: -3.9885, the 5th percentile of the genome-wide distribution)
- `calculate_mu_params.gerp_upper_cutoff`: Maximum GERP score for a site to be included (default: 2.6607, the 95th percentile of the genome-wide distribution)
- `calculate_mu_params.downsampling_level`: Downsampling level used for the mutation rate calculation
- `calculate_mu_params.most_severe_consequence`: List of most severe transcript consequences used as neutral sites for the mutation rate calculation

### `build_models_params`

Parameters used when building the plateau and coverage correction models.

- `build_models_params.low_cov_cutoff`: Lower AN% cutoff; sites at or below this value are excluded from model training
- `build_models_params.high_cov_cutoff`: AN% cutoff separating high-coverage (plateau model) from low-coverage (coverage correction model) training sites
- `build_models_params.upper_cov_cutoff`: Upper AN% cutoff; sites above this value are excluded from the high-coverage (plateau) model training set but may still be included in the low-coverage model (null if no upper bound)

### `apply_models_params`

Parameters used when applying the models to compute expected variant counts.

- `apply_models_params.low_cov_cutoff`: Lower AN% cutoff; sites at or below this value are excluded from expected variant calculations
- `apply_models_params.high_cov_cutoff`: AN% cutoff separating high-coverage from low-coverage sites when applying the model
- `apply_models_params.plateau_models`: Dictionary mapping (CpG context, genomic region) pairs to linear regression coefficients. Plateau models relate per-context mutation rates to the proportion of possible synonymous sites at which a variant was observed at high-coverage sites (AN% ≥ 90%). Separate models were fit for CpG transitions vs. all other sites (transversions and non-CpG transitions) and for each genomic region (autosomes/PAR, chrX non-PAR, and chrY non-PAR)
- `apply_models_params.coverage_model`: Array of coverage correction model coefficients used for expected variant count calculations at low-coverage sites (AN% between 20% and 90%)
- `apply_models_params.log10_coverage`: Boolean indicating whether log10 transformation was applied to coverage values in the coverage correction model

### `downsamplings`

Struct mapping each genetic ancestry group to an array of downsampling levels (as integers) corresponding to the indices in `gen_anc_obs` and `gen_anc_exp` arrays in the row fields.

- `downsamplings.global`: Array of downsampling levels for the full gnomAD v4 exomes dataset
- `downsamplings.afr`: Array of downsampling levels for the African/African American genetic ancestry group
- `downsamplings.amr`: Array of downsampling levels for the Admixed American genetic ancestry group
- `downsamplings.eas`: Array of downsampling levels for the East Asian genetic ancestry group
- `downsamplings.nfe`: Array of downsampling levels for the non-Finnish European genetic ancestry group
- `downsamplings.sas`: Array of downsampling levels for the South Asian genetic ancestry group

### Other globals

- `max_af`: Maximum alternate allele frequency cutoff used to define observed variants
- `sd_raw_z`: Struct of standard deviations of raw Z-scores used to normalize `z_raw` → `z_score` for each constraint group
  - `sd_raw_z.syn`: Standard deviation of raw Z-scores for synonymous variants
  - `sd_raw_z.mis`: Standard deviation of raw Z-scores for missense variants
  - `sd_raw_z.lof_hc_lc`: Standard deviation of raw Z-scores for high and low confidence pLoF variants
  - `sd_raw_z.lof`: Standard deviation of raw Z-scores for high confidence pLoF variants
## Key Fields

- `gene`: Gene name
- `gene_id`: Ensembl gene ID
- `transcript`: Ensembl or RefSeq transcript ID (GENCODE v39)
- `canonical`: Boolean indicator as to whether the transcript is the canonical transcript for the gene
- `mane_select`: Boolean indicator as to whether the transcript is the MANE Select transcript for the gene

## General Fields
- `transcript_version`: Ensembl or RefSeq transcript version
- `transcript_type`: Transcript biotype from [Gencode](https://www.gencodegenes.org/pages/biotypes.html)
- `transcript_level`: Transcript level from [Gencode](https://www.gencodegenes.org/pages/data_format.html)
- `chromosome`: Chromosome where gene is located
- `start_position`: Start position of the transcript
- `end_position`: End position of the transcript
- `cds_length`: Length of the coding sequences (CDS) in the transcript
- `num_coding_exons`: Number of coding exons in the transcript
- `gene_quality_metrics.exome_prop_bp_AN90`:  Proportion of coding bases in gene with allele number percent (AN%) (percent of total possible AN observed at site) greater than 90%
- `gene_quality_metrics.exome_mean_AS_MQ`: Mean AS_MQ (allele-specific root mean square of the mapping quality of reads) across SNV sites in the coding sequence of the gene
- `gene_quality_metrics.exome_prop_segdup`: Proportion of coding bases in gene that overlap a segmental duplication
- `gene_quality_metrics.exome_prop_LCR`: Proportion of coding bases in gene that overlap a low-complexity region

## Flags
- `gene_flags`: Quality flags for gene based on gnomAD v4 exome data. One of:
  - `low_exome_coverage`: Less than 10% of coding base pairs in gene have allele number percent (AN%) greater than or equal to 90%
  - `low_exome_mapping_quality`: Mean value of AS_MQ across gene is less than 50.
- `constraint_flags`: Reason transcript is considered an outlier for constraint metrics. One of:
  - `no_variants`: Zero observed synonymous, missense, pLoF variants
  - `no_exp_lof`: Zero expected pLoF variants
  - `outlier_lof`: Number of pLoF variants is significantly different than expectation
  - `no_exp_mis`: Zero expected missense variants
  - `outlier_mis`: Number of missense variants is significantly different than expectation
  - `no_exp_syn`: Zero expected synonymous variants
  - `outlier_syn`: Number of synonymous variants is significantly different than expectation

## Synonymous

- `syn.mu`: Mutation rate summed across all synonymous variants in transcript
- `syn.possible`: Number of possible synonymous variants in transcript
- `syn.obs`: Number of observed synonymous variants in transcript
- `syn.exp`: Number of expected synonymous variants in transcript
- `syn.oe`: Observed to expected ratio for synonymous variants in transcript (`syn.obs` divided by `syn.exp`)
- `syn.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for synonymous variants
- `syn.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for synonymous variants
- `syn.z_raw`: Raw (unnormalized) Z-score for synonymous variants in transcript. Computed as the signed square root of the chi-squared deviation of observed from expected counts. Extreme values indicate likely data quality issues.
- `syn.z_score`: Normalized Z-score for synonymous variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values indicate likely data quality issues.
- `syn.gen_anc_obs.global`: Array of observed synonymous variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `syn.gen_anc_obs.afr`: Array of observed synonymous variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `syn.gen_anc_obs.amr`: Array of observed synonymous variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `syn.gen_anc_obs.eas`: Array of observed synonymous variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `syn.gen_anc_obs.nfe`: Array of observed synonymous variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `syn.gen_anc_obs.sas`: Array of observed synonymous variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`
- `syn.gen_anc_exp.global`: Array of expected synonymous variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `syn.gen_anc_exp.afr`: Array of expected synonymous variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `syn.gen_anc_exp.amr`: Array of expected synonymous variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `syn.gen_anc_exp.eas`: Array of expected synonymous variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `syn.gen_anc_exp.nfe`: Array of expected synonymous variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `syn.gen_anc_exp.sas`: Array of expected synonymous variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`

## Missense

- `mis.mu`: Mutation rate summed across all missense variants in transcript
- `mis.possible`: Number of possible missense variants in transcript
- `mis.obs`: Number of observed missense variants in transcript
- `mis.exp`: Number of expected missense variants in transcript
- `mis.oe`: Observed to expected ratio for missense variants in transcript (`mis.obs` divided by `mis.exp`)
- `mis.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.z_raw`: Raw (unnormalized) Z-score for missense variants in transcript. Computed as the signed square root of the chi-squared deviation of observed from expected counts. Extreme values indicate likely data quality issues.
- `mis.z_score`: Normalized Z-score for missense variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values indicate likely data quality issues.
- `mis.gen_anc_obs.global`: Array of observed missense variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `mis.gen_anc_obs.afr`: Array of observed missense variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `mis.gen_anc_obs.amr`: Array of observed missense variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `mis.gen_anc_obs.eas`: Array of observed missense variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `mis.gen_anc_obs.nfe`: Array of observed missense variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `mis.gen_anc_obs.sas`: Array of observed missense variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`
- `mis.gen_anc_exp.global`: Array of expected missense variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `mis.gen_anc_exp.afr`: Array of expected missense variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `mis.gen_anc_exp.amr`: Array of expected missense variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `mis.gen_anc_exp.eas`: Array of expected missense variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `mis.gen_anc_exp.nfe`: Array of expected missense variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `mis.gen_anc_exp.sas`: Array of expected missense variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`

## Loss-of-Function (High and Low Confidence)

- `lof_hc_lc.mu`: Mutation rate summed across all possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.possible`: Number of possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.obs`: Number of observed high and low confidence predicted loss-of-function (pLoF) variants in transcript
- `lof_hc_lc.exp`: Number of expected high and low confidence pLoF variants in transcript
- `lof_hc_lc.oe`: Observed to expected (`oe`) ratio for high and low confidence pLoF variants (`lof_hc_lc.obs` divided by `lof_hc_lc.exp`)
- `lof_hc_lc.oe_ci.lower`: Lower bound of 90% confidence interval (CI) for observed to expected (`oe`) ratio for high and low confidence pLoF variants
- `lof_hc_lc.oe_ci.upper`: LOEUF: Upper bound of 90% confidence interval for `oe` ratio for high and low confidence pLoF variants (lower values indicate more constrained)
- `lof_hc_lc.oe_ci.upper_rank`: Transcript's rank of upper bound `oe` CI value for pLoF variants compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.oe_ci.upper_bin_percentile`: Transcript percentile associated with LOEUF score for transcript. For example, if the gene percentile is 85, then the transcript is more highly constrained against predicted loss-of-function variation than 85% of transcripts.
- `lof_hc_lc.oe_ci.upper_bin_decile`: Decile bin of upper bound of 90% CI `oe` for pLoF variants for given transcript (lower values indicate more constrained).  This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.oe_ci.upper_bin_sextile`: Sextile bin of upper bound of 90% CI `oe` for pLoF variants for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.z_raw`: Raw (unnormalized) Z-score for high and low confidence pLoF variants in transcript. Computed as the signed square root of the chi-squared deviation of observed from expected counts.
- `lof_hc_lc.z_score`: Normalized Z-score for high and low confidence pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `lof_hc_lc.pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% oe pLoF ratio; computed from high and low confidence pLoF gnomAD data)
- `lof_hc_lc.pRec`: Probability that transcript falls into distribution of recessive genes (~71% oe pLoF ratio; computed from high and low confidence pLoF gnomAD data)
- `lof_hc_lc.pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% oe pLoF ratio; computed from high and low confidence pLoF gnomAD data)
- `lof_hc_lc.gen_anc_obs.global`: Array of observed HC+LC pLoF variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `lof_hc_lc.gen_anc_obs.afr`: Array of observed HC+LC pLoF variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `lof_hc_lc.gen_anc_obs.amr`: Array of observed HC+LC pLoF variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `lof_hc_lc.gen_anc_obs.eas`: Array of observed HC+LC pLoF variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `lof_hc_lc.gen_anc_obs.nfe`: Array of observed HC+LC pLoF variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `lof_hc_lc.gen_anc_obs.sas`: Array of observed HC+LC pLoF variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`
- `lof_hc_lc.gen_anc_exp.global`: Array of expected HC+LC pLoF variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `lof_hc_lc.gen_anc_exp.afr`: Array of expected HC+LC pLoF variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `lof_hc_lc.gen_anc_exp.amr`: Array of expected HC+LC pLoF variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `lof_hc_lc.gen_anc_exp.eas`: Array of expected HC+LC pLoF variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `lof_hc_lc.gen_anc_exp.nfe`: Array of expected HC+LC pLoF variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `lof_hc_lc.gen_anc_exp.sas`: Array of expected HC+LC pLoF variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`

## Loss-of-Function (High Confidence Only)

- `lof.mu`: Mutation rate summed across all possible high confidence pLoF variants in transcript
- `lof.possible`: Number of possible high confidence pLoF variants in transcript
- `lof.obs`: Number of observed high confidence predicted loss-of-function (pLoF) variants in transcript
- `lof.exp`: Number of expected high confidence pLoF variants in transcript
- `lof.oe`: Observed to expected ratio for high confidence pLoF variants in transcript (`lof.obs` divided by `lof.exp`)
- `lof.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants
- `lof.oe_ci.upper`: LOEUF: upper bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants (lower values indicate more constrained)
- `lof.oe_ci.upper_rank`: Transcript's rank of LOEUF value compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.oe_ci.upper_bin_percentile`: Transcript percentile associated with LOEUF score for transcript. For example, if the gene percentile is 85, then the transcript is more highly constrained against predicted loss-of-function variation than 85% of transcripts.
- `lof.oe_ci.upper_bin_decile`: Decile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.oe_ci.upper_bin_sextile`: Sextile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.z_raw`: Raw (unnormalized) Z-score for high confidence pLoF variants in transcript. Computed as the signed square root of the chi-squared deviation of observed from expected counts.
- `lof.z_score`: Normalized Z-score for high confidence pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `lof.pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.pRec`: Probability that transcript falls into distribution of recessive genes (~71% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.gen_anc_obs.global`: Array of observed high confidence pLoF variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `lof.gen_anc_obs.afr`: Array of observed high confidence pLoF variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `lof.gen_anc_obs.amr`: Array of observed high confidence pLoF variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `lof.gen_anc_obs.eas`: Array of observed high confidence pLoF variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `lof.gen_anc_obs.nfe`: Array of observed high confidence pLoF variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `lof.gen_anc_obs.sas`: Array of observed high confidence pLoF variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`
- `lof.gen_anc_exp.global`: Array of expected high confidence pLoF variant counts at each downsampling level for the full gnomAD v4 exomes dataset; index i corresponds to `downsamplings.global[i]`
- `lof.gen_anc_exp.afr`: Array of expected high confidence pLoF variant counts at each downsampling level for the African/African American genetic ancestry group; index i corresponds to `downsamplings.afr[i]`
- `lof.gen_anc_exp.amr`: Array of expected high confidence pLoF variant counts at each downsampling level for the Admixed American genetic ancestry group; index i corresponds to `downsamplings.amr[i]`
- `lof.gen_anc_exp.eas`: Array of expected high confidence pLoF variant counts at each downsampling level for the East Asian genetic ancestry group; index i corresponds to `downsamplings.eas[i]`
- `lof.gen_anc_exp.nfe`: Array of expected high confidence pLoF variant counts at each downsampling level for the non-Finnish European genetic ancestry group; index i corresponds to `downsamplings.nfe[i]`
- `lof.gen_anc_exp.sas`: Array of expected high confidence pLoF variant counts at each downsampling level for the South Asian genetic ancestry group; index i corresponds to `downsamplings.sas[i]`
