# Description of fields in release mutation rate table

Table containing per-base mutation rates estimated from intronic and intergenic variants in the gnomAD v4 genomes, grouped by trinucleotide context, reference/alternate allele, and methylation level. For more information on how these mutation rates are calculated, see [Karczewski et al. *Nature* 2020](https://www.nature.com/articles/s41586-020-2308-7).

## Globals

| Global | Description |
| :--- | :--- |
| **version** | gnomAD constraint release version. |

### `calculate_mu_params`

Parameters used when computing the mutation rate. Only autosomal sites with Ensembl VEP most severe consequence in `most_severe_consequence` and GERP scores between `gerp_lower_cutoff` and `gerp_upper_cutoff` are included. Sites must also have genome mean coverage between `min_cov` and `max_cov`, and observed variants must have AC <= `ac_cutoff` at the specified `downsampling_level`.

| Global | Description |
| :--- | :--- |
| **calculate_mu_params.ac_cutoff** | Maximum allele count for a variant to be counted as observed in the mutation rate calculation (default: 5). |
| **calculate_mu_params.min_cov** | Minimum genome mean coverage required for a site to be included (default: 15). |
| **calculate_mu_params.max_cov** | Maximum genome mean coverage allowed for a site to be included (default: 60). |
| **calculate_mu_params.gerp_lower_cutoff** | Minimum GERP score for a site to be included. Default of -3.9885 is the 5th percentile of the genome-wide GERP score distribution, precalculated on the GRCh37 context Table. |
| **calculate_mu_params.gerp_upper_cutoff** | Maximum GERP score for a site to be included. Default of 2.6607 is the 95th percentile of the genome-wide GERP score distribution, precalculated on the GRCh37 context Table. |
| **calculate_mu_params.downsampling_level** | Number of genome samples in the downsampling used to compute `mu` (default: 1000). |
| **calculate_mu_params.most_severe_consequence** | Ensembl VEP most severe consequences used to identify neutral sites for the mutation rate calculation (intronic and intergenic variants). |

## Key fields

| Field Name | Description |
| :--- | :--- |
| **context** | Trinucleotide sequence context centered on the mutated base (e.g., "ACG"). |
| **ref** | Reference allele at the central position. |
| **alt** | Alternate allele at the central position. |
| **methylation_level** | Discrete methylation category for the site (integer). Used to stratify mutation rates because CpG mutation rates vary with methylation status. |

## Row fields

| Field Name | Description |
| :--- | :--- |
| **mu** | Per-base SNP mutation rate for this trinucleotide context, allele, and methylation class. Computed at the `calculate_mu_params.downsampling_level`-sample genome downsampling level by scaling the fraction of possible variants observed by a correction factor derived from the total per-generation mutation rate (1.2 x 10^-8). |
| **cpg** | Whether the mutation occurs at a CpG dinucleotide. |
| **transition** | Whether the mutation is a transition (A<->G or C<->T) as opposed to a transversion. |
| **mutation_type** | Categorical mutation class: one of "CpG", "non-CpG transition", or "transversion". |
