# Constraint Field Descriptions

Descriptions of columns in the gnomAD v4.1.1 gene constraint metrics tsv. Descriptions also apply to rows in Hail Tables, where a "." in the field name indicates a struct. All constraint metrics were calculated using the gnomAD v4.1.1 exomes.

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
- `cds_length`: Length of the coding sequences (CDS) in the transcript
transcript
- `num_coding_exons`: Number of coding exons in the transcript
- `gene_quality_metrics.exome_prop_bp_AN90`:  Proportion of coding bases in gene with median allele number percent (AN%) (percent of total possible AN observed at site) greater than 90%
- `gene_quality_metrics.exome_mean_AS_MQ`: Mean value of AS_MQ (allele-specific root mean square of the mapping quality of reads across all samples) across gene
- `gene_quality_metrics.exome_prop_segdup`: Proportion of coding bases in gene that overlap a segmental duplication
- `gene_quality_metrics.exome_prop_LCR`: Proportion of coding bases in gene that overlap a low-complexity region

## Flags
- `gene_flags`: Quality flags for gene based on gnomAD v4 exome data. One of:
  - `low_exome_coverage`: Proportion of coding base pairs in gene with median AN% is less than 10%
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
- `syn.z_score`: Raw Z-score for synonymous variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values of syn.z_score indicate likely data quality issues.
- `syn.z_score`: Z-score for synonymous variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values of syn.z_score indicate likely data quality issues.
- `gen_anc_obs.global`: Array of observed values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_obs.afr`: Array of observed values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.amr`: Array of observed values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.eas`: Array of observed values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.nfe`: Array of observed values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.sas`: Array of observed values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.global`: Array of expected values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_exp.afr`: Array of expected values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.amr`: Array of expected values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.eas`: Array of expected values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.nfe`: Array of expected values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.sas`: Array of expected values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here

## Missense

- `mis.mu`: Mutation rate summed across all missense variants in transcript
- `mis.possible`: Number of possible missense variants in transcript
- `mis.obs`: Number of observed missense variants in transcript
- `mis.exp`: Number of expected missense variants in transcript
- `mis.oe`: Observed to expected ratio for missense variants in transcript (`mis.obs` divided by `mis.exp`)
- `mis.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.z_score`: Raw Z-score for missense variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values of mis.z_score indicate likely data quality issues.
- `mis.z_score`: Z-score for missense variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values of mis.z_score indicate likely data quality issues.
- `gen_anc_obs.global`: Array of observed values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_obs.afr`: Array of observed values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.amr`: Array of observed values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.eas`: Array of observed values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.nfe`: Array of observed values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.sas`: Array of observed values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.global`: Array of expected values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_exp.afr`: Array of expected values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.amr`: Array of expected values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.eas`: Array of expected values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.nfe`: Array of expected values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.sas`: Array of expected values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here

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
- `lof_hc_lc.z_score`: Raw Z-score for high and low confidence pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `lof_hc_lc.z_score`: Z-score for high and low confidence pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `gen_anc_obs.global`: Array of observed values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_obs.afr`: Array of observed values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.amr`: Array of observed values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.eas`: Array of observed values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.nfe`: Array of observed values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.sas`: Array of observed values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.global`: Array of expected values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_exp.afr`: Array of expected values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.amr`: Array of expected values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.eas`: Array of expected values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.nfe`: Array of expected values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.sas`: Array of expected values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here

## Loss-of-Function (High Confidence Only)

- `lof.mu`: Mutation rate summed across all possible high confidence pLoF variants in transcript
- `lof.possible`: Number of possible high confidence pLoF variants in transcript
- `lof.obs`: Number of observed high confidence predicted loss-of-function (pLoF) variants in transcript
- `lof.exp`: Number of expected high confidence pLoF variants in transcript
- `lof.oe`: Observed to expected ratio for high confidence pLoF variants in transcript (`lof.obs` divided by `lof.exp`)
- `lof.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants
- `lof.oe_ci.upper`: LOEUF: upper bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants (lower values indicate more constrained)
- `lof.oe_ci.upper_rank`: Transcript's rank of LOEUF value compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.oe_ci.upper_bin_percentile`: Transcript percentile associated with LOEUF score for transcript. For example, if the gene percentile is 85, then the transcript is more highly constrained against predicted loss-of-function variation than 85% of transcripts.
- `lof.oe_ci.upper_bin_decile`: Decile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.oe_ci.upper_bin_sextile`: Sextile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.z_score`: Raw Z-score for pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `lof.z_score`: Z-score for pLoF variants in transcript. Higher (more positive) Z-scores indicate that the transcript is more intolerant of variation (more constrained).
- `pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `pRec`: Probability that transcript falls into distribution of recessive genes (~71% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `gen_anc_obs.global`: Array of observed values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_obs.afr`: Array of observed values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.amr`: Array of observed values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.eas`: Array of observed values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.nfe`: Array of observed values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_obs.sas`: Array of observed values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.global`: Array of expected values across different strata using the full gnomAD v4 exomes dataset (global). Indices in array are described in TODO: add global field here
- `gen_anc_exp.afr`: Array of expected values across different strata using the African/African American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.amr`: Array of expected values across different strata using the Admixed American genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.eas`: Array of expected values across different strata using the East Asian genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.nfe`: Array of expected values across different strata using the non-Finnish European genetic ancestry group. Indices in array are described in TODO: add global field here
- `gen_anc_exp.sas`: Array of expected values across different strata using the South Asian genetic ancestry group. Indices in array are described in TODO: add global field here
