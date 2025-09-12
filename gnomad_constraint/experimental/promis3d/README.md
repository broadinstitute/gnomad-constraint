# Proemis3D Pipeline

The Proemis3D (Protein Missense Constraint in 3D) pipeline is an experimental module for analyzing protein missense constraint using 3D structural information from AlphaFold2. This pipeline integrates genomic variant data with protein structural data to identify regions of high missense constraint in 3D space.

## Overview

The Proemis3D pipeline combines:
- **Genomic variant data** from gnomAD and other sources
- **Protein structural data** from AlphaFold2
- **Functional annotations** from various databases
- **Statistical methods** for constraint analysis

To identify regions of proteins that are highly intolerant to missense variation in 3D space.

## Key Features

- **3D Constraint Analysis**: Identifies regions of high missense constraint based on 3D protein structure
- **AlphaFold2 Integration**: Uses AlphaFold2 structural predictions and confidence scores
- **Multiple Annotation Sources**: Integrates data from COSMIS, InterPro, ClinVar, and other databases
- **Statistical Rigor**: Implements both greedy and forward algorithms for constraint region identification
- **Scalable Processing**: Built on Hail for large-scale genomic data processing

## Directory Structure

```
promis3d/
├── README.md                 # This file
├── __init__.py              # Package initialization
├── constants.py             # Pipeline constants and gene lists
├── data_import.py           # Data import and processing functions
├── proemis_3d.py           # Main pipeline script
├── resources.py             # Resource definitions and paths
└── utils.py                 # Utility functions and algorithms
```

## Core Components

### 1. Constants (`constants.py`)
- **MIN_EXP_MIS**: Minimum expected missense variants for constraint calculation (16)
- **Gene Lists**: HI genes, severe HI genes, and gene categories
- **Configuration**: Pipeline parameters and thresholds

### 2. Data Import (`data_import.py`)
- **COSMIS Scores**: Protein constraint scores from multiple structure sources
- **Varity Data**: Variant effect predictions
- **MTR3D Data**: Missense tolerance ratio in 3D
- **InterPro Annotations**: Protein domain and functional annotations
- **Kaplanis Variants**: Developmental delay de novo missense variants
- **Fu Variants**: Autism spectrum disorder de novo variants
- **ClinVar Data**: Clinical variant interpretations and significance
- **Constraint Metrics**: Gene-level constraint statistics (synonymous, missense, LoF)
- **MTR Data**: Missense tolerance ratio annotations
- **RMC Data**: Regional missense constraint metrics
- **Context Data**: Variant context and coverage information
- **Genetics Gym Scores**: AI-based missense prediction scores
- **REVEL Scores**: Rare Exome Variant Ensemble Learner predictions

#### Key Functions:
- `import_cosmis_score_data()`: Import COSMIS constraint scores from AlphaFold2, PDB, and Swiss Model structures
- `import_varity_data()`: Import VariTY variant effect prediction scores
- `import_mtr3d_data()`: Import MTR3D (Missense Tolerance Ratio in 3D) scores
- `import_mtr_data()`: Import MTR (Missense Tolerance Ratio) annotations
- `import_kaplanis_variants()`: Process developmental delay de novo variants with GRCh37/GRCh38 liftover
- `get_kaplanis_sig_gene_annotations()`: Get significant gene set annotations from Kaplanis study
- `import_fu_variants()`: Import autism spectrum disorder de novo variants
- `import_interpro_annotations()`: Import InterPro protein domain and functional annotations
- `process_clinvar_ht()`: Process ClinVar data and filter to missense variants
- `process_gnomad_site_ht()`: Process gnomAD site-level variant data
- `process_pext_base_ht()`: Process base-level expression data
- `process_pext_annotation_ht()`: Process annotation-level expression data
- `process_gnomad_de_novo_ht()`: Process gnomAD de novo variant data
- `process_rmc_ht()`: Process Regional Missense Constraint data with p-values and confidence intervals
- `process_constraint_metrics_ht()`: Process gene-level constraint metrics (syn, mis, lof)
- `process_context_ht()`: Process variant context data with coverage and frequency information
- `process_genetics_gym_missense_scores_ht()`: Process AI-based missense prediction scores with percentiles
- `import_revel_ht()`: Import REVEL (Rare Exome Variant Ensemble Learner) scores

### 4. Utils (`utils.py`)
- **FASTA Processing**: Convert GENCODE FASTA files to Hail tables
- **AlphaFold2 Processing**: Extract sequences, distance matrices, and confidence scores
- **Constraint Algorithms**: Greedy and forward algorithms for region identification
- **Annotation Functions**: Combine variant and residue-level annotations

### 5. Main Pipeline (`proemis_3d.py`)
- **Pipeline Orchestration**: Coordinates all pipeline steps
- **Resource Management**: Handles data dependencies and checkpointing
- **Command Line Interface**: Provides CLI for pipeline execution

## Key Algorithms

### Greedy Algorithm
Identifies the most intolerant region by iteratively selecting residues with the lowest upper bound of the observed/expected (OE) confidence interval.

### Forward Algorithm
Uses a forward selection approach with Akaike Information Criterion (AIC) to identify optimal constraint regions.

### 3D Constraint Calculation
- Calculates distance matrices from AlphaFold2 structures
- Identifies spatially proximal residues
- Computes constraint metrics for 3D regions

## Usage

### Basic Pipeline Execution

```python
import hail as hl
from gnomad_constraint.experimental.proemis3d import get_proemis3d_resources

# Initialize Hail
hl.init()

# Get pipeline resources
resources = get_proemis3d_resources(
    version="4.1",
    overwrite=False,
    test=False
)

# Run the pipeline
resources.run()
```

### Command Line Interface

```bash
python proemis_3d.py --version 4.1 --test
```

### Key Parameters

- `version`: gnomAD version (2.1.1 or 4.1)
- `test`: Run in test mode with smaller datasets
- `overwrite`: Overwrite existing outputs

## Command Line Parameters

The Proemis3D pipeline supports extensive command-line parameters for fine-grained control over execution. Here's a complete reference:

### Basic Parameters

#### `--version`
- **Type**: String
- **Default**: `4.1` (current version)
- **Description**: Which version of the resource Tables will be used
- **Options**: `2.1.1`, `4.1`

#### `--test`
- **Type**: Boolean flag
- **Description**: Whether to run a test instead of the full pipeline
- **Effect**: Uses smaller test datasets and filters to specific test transcript (ENST00000372435)

#### `--overwrite`
- **Type**: Boolean flag
- **Description**: Whether to overwrite existing output files
- **Effect**: Forces regeneration of existing intermediate and output files

### Data Processing Steps

#### `--convert-gencode-fastn-to-ht`
- **Type**: Boolean flag
- **Description**: Import and pre-process GENCODE transcripts FASTA file as a Hail Table
- **Output**: GENCODE transcripts Hail Table with sequence data

#### `--convert-gencode-fasta-to-ht`
- **Type**: Boolean flag
- **Description**: Import and pre-process GENCODE translations FASTA file as a Hail Table
- **Output**: GENCODE translations Hail Table with protein sequence data

#### `--read-af2-sequences`
- **Type**: Boolean flag
- **Description**: Process AlphaFold2 structures from GCS bucket into a Hail Table
- **Output**: AlphaFold2 sequences Hail Table
- **Mode**: `sequence`

#### `--compute-af2-distance-matrices`
- **Type**: Boolean flag
- **Description**: Compute distance matrices for AlphaFold2 structures
- **Output**: AlphaFold2 distance matrices Hail Table
- **Mode**: `distance_matrix`

#### `--extract-af2-plddt`
- **Type**: Boolean flag
- **Description**: Extract pLDDT (per-residue confidence) scores from AlphaFold2 structures
- **Output**: AlphaFold2 pLDDT scores Hail Table
- **Mode**: `plddt`

#### `--extract-af2-pae`
- **Type**: Boolean flag
- **Description**: Extract pAE (predicted aligned error) scores from AlphaFold2 structures
- **Output**: AlphaFold2 PAE matrices Hail Table
- **Mode**: `pae`

#### `--gencode-alignment`
- **Type**: Boolean flag
- **Description**: Join GENCODE translations and AlphaFold2 structures based on sequence
- **Output**: Matched GENCODE-AlphaFold2 Hail Table

#### `--get-gencode-positions`
- **Type**: Boolean flag
- **Description**: Create GENCODE positions Hail Table with genomic coordinates
- **Output**: GENCODE positions Hail Table

### Constraint Analysis

#### `--run-greedy`
- **Type**: Boolean flag
- **Description**: Execute the greedy algorithm for constraint region identification
- **Output**: Greedy algorithm results Hail Table

#### `--run-forward`
- **Type**: Boolean flag
- **Description**: Execute the forward algorithm for constraint region identification
- **Output**: Forward algorithm results Hail Table

#### `--min-exp-mis`
- **Type**: Integer
- **Default**: `16`
- **Description**: Minimum expected number of missense variants to consider for constraint algorithms
- **Effect**: Filters regions with insufficient expected missense variants

### Output Generation

#### `--write-per-variant`
- **Type**: Boolean flag
- **Description**: Generate per-variant annotated Hail Table with comprehensive annotations
- **Output**: Fully annotated per-variant Hail Table
- **Dependencies**: Requires forward algorithm results

#### `--write-per-missense-variant`
- **Type**: Boolean flag
- **Description**: Generate per-variant annotated Hail Table filtered to missense variants only
- **Output**: Missense-only per-variant Hail Table
- **Dependencies**: Requires per-variant Hail Table

#### `--write-per-residue`
- **Type**: Boolean flag
- **Description**: Generate per-residue Hail Table from per-variant data
- **Output**: Per-residue annotated Hail Table
- **Dependencies**: Requires per-variant Hail Table

#### `--write-per-region`
- **Type**: Boolean flag
- **Description**: Generate per-region Hail Table from per-residue data
- **Output**: Per-region annotated Hail Table
- **Dependencies**: Requires per-residue Hail Table

#### `--create-missense-viewer-input-ht`
- **Type**: Boolean flag
- **Description**: Create missense viewer input Hail Table for visualization
- **Output**: Formatted data for web-based visualization
- **Dependencies**: Requires forward algorithm results and GENCODE positions

### Performance Tuning

#### `--all-snv-n-partitions`
- **Type**: Integer
- **Default**: `5000`
- **Description**: Number of partitions to use for the all possible SNVs Hail Table
- **Effect**: Controls memory usage and processing speed for large datasets

### Usage Examples

#### Complete Pipeline
```bash
python proemis_3d.py --version 4.1 --overwrite
```

#### Test Run
```bash
python proemis_3d.py --version 4.1 --test --overwrite
```

#### Individual Steps
```bash
# Process GENCODE data
python proemis_3d.py --convert-gencode-fastn-to-ht --convert-gencode-fasta-to-ht

# Process AlphaFold2 data
python proemis_3d.py --read-af2-sequences --compute-af2-distance-matrices --extract-af2-plddt --extract-af2-pae

# Run constraint analysis
python proemis_3d.py --run-greedy --run-forward --min-exp-mis 20

# Generate outputs
python proemis_3d.py --write-per-variant --write-per-residue --write-per-region
```

#### Custom Configuration
```bash
python proemis_3d.py \
    --version 4.1 \
    --test \
    --overwrite \
    --min-exp-mis 25 \
    --all-snv-n-partitions 10000 \
    --run-forward \
    --write-per-variant \
    --write-per-residue
```

### Pipeline Dependencies

The pipeline steps have specific dependencies that must be satisfied:

1. **GENCODE Processing**: `--convert-gencode-fastn-to-ht` and `--convert-gencode-fasta-to-ht` can run independently
2. **AlphaFold2 Processing**: All AF2 steps can run independently
3. **Alignment**: `--gencode-alignment` requires both GENCODE and AF2 sequence data
4. **Positions**: `--get-gencode-positions` requires GENCODE data and alignment
5. **Constraint Analysis**: `--run-greedy` and `--run-forward` require positions and RMC data
6. **Output Generation**: Each output step depends on its prerequisite data

### Resource Management

- **Checkpointing**: Intermediate results are automatically checkpointed
- **Temporary Files**: Uses `gs://gnomad-tmp-4day` for temporary storage
- **Logging**: Pipeline logs are written to `/proemis_3d.log`
- **Memory Management**: Large operations use repartitioning for memory efficiency

## Data Requirements

### Input Data
- **GENCODE**: Transcript and translation FASTA files
- **AlphaFold2**: Protein structures and confidence scores
- **gnomAD**: Variant data and constraint metrics
- **Annotations**: COSMIS, InterPro, ClinVar, and other functional annotations

### Output Data
- **Per-SNV Tables**: Annotated variant-level data
- **Per-Residue Tables**: Residue-level constraint metrics
- **Per-Region Tables**: 3D constraint regions
- **Viewer Input**: Formatted data for visualization

## Key Functions

### Data Processing
- `convert_fasta_to_table()`: Convert FASTA files to Hail tables
- `process_af2_structures()`: Process AlphaFold2 structural data
- `get_gencode_positions()`: Extract genomic positions for transcripts

### Constraint Analysis
- `run_greedy()`: Execute greedy constraint algorithm
- `run_forward()`: Execute forward constraint algorithm
- `determine_regions_with_min_oe_upper()`: Identify constraint regions

### Annotation
- `annotate_snvs_with_variant_level_data()`: Add variant-level annotations
- `annotate_proemis3d_with_af2_metrics()`: Add AlphaFold2 metrics
- `create_per_snv_combined_ht()`: Create comprehensive annotation tables

## Dependencies

- **Hail**: Genomic data processing
- **Biopython**: Protein structure analysis
- **NumPy/Pandas**: Numerical computing
- **PySpark**: Distributed computing

## Output Structure

The pipeline generates several key outputs:

1. **Variant-Level Data**: All possible SNVs with comprehensive annotations
2. **Residue-Level Data**: Per-residue constraint metrics and AlphaFold2 scores
3. **Region-Level Data**: 3D constraint regions with statistical significance
4. **Viewer Data**: Formatted data for web-based visualization

## Statistical Methods

- **Observed/Expected Ratios**: Compare observed to expected missense variants
- **Confidence Intervals**: Chi-squared based confidence intervals
- **AIC Selection**: Model selection using Akaike Information Criterion
- **3D Distance Metrics**: Spatial proximity calculations

## Performance Considerations

- **Checkpointing**: Intermediate results are checkpointed for fault tolerance
- **Partitioning**: Data is partitioned for efficient processing
- **Caching**: Frequently accessed data is cached in memory
- **Parallel Processing**: Uses Spark for distributed computation

## Testing

The pipeline includes test mode with smaller datasets:
- **Test Transcript**: ENST00000372435
- **Test UniProt**: P60891
- **Reduced Data**: Smaller subsets for development and testing

## Future Development

- **Additional Algorithms**: New constraint identification methods
- **Enhanced Annotations**: Integration with additional data sources
- **Visualization Tools**: Interactive 3D constraint visualization
- **Performance Optimization**: Improved scalability and efficiency

## Citation

If you use this pipeline in your research, please cite the relevant gnomAD and AlphaFold2 papers, as well as any specific Proemis3D methodology papers.

## Contact

For questions or issues with the Proemis3D pipeline, please contact the gnomAD team or create an issue in the repository.

## Detailed Output Data Structure

### Fully Annotated Hail Table Schema

The main output of the Proemis3D pipeline is a comprehensive Hail Table with the following structure:

#### Global Fields
- **None** (no global fields)

#### Row Fields

##### Basic Identifiers
- `locus`: Genomic position (locus<GRCh38>)
- `alleles`: Variant alleles (array<str>)
- `transcript_id`: GENCODE transcript ID (str)
- `uniprot_id`: UniProt protein ID (str)
- `gene_id`: GENCODE gene ID (str)
- `gene_symbol`: Gene symbol (str)

##### Transcript and Gene Metadata
- `canonical`: Whether transcript is canonical (bool)
- `mane_select`: Whether transcript is MANE select (bool)
- `transcript_biotype`: Transcript biotype (str)
- `most_severe_consequence`: Most severe consequence (str)
- `cds_len_mismatch`: CDS length mismatch flag (bool)
- `cds_len_not_div_by_3`: CDS length not divisible by 3 flag (bool)

##### Gene Classification Flags
- `is_phaplo_gene`: Haploinsufficiency gene flag (bool)
- `is_ptriplo_gene`: Triplosensitivity gene flag (bool)
- `is_hi_gene`: Haploinsufficiency gene flag (bool)
- `hi_gene_category`: HI gene category (str)
- `one_uniprot_per_transcript`: One UniProt per transcript flag (bool)
- `one_transcript_per_gene`: One transcript per gene flag (bool)

##### Variant Level Annotations (`variant_level_annotations`)

**Context and Basic Info:**
- `context`: Sequence context (str)
- `ref`: Reference allele (str)
- `alt`: Alternative allele (str)
- `was_flipped`: Whether variant was flipped (bool)
- `transition`: Whether variant is a transition (bool)
- `cpg`: Whether variant is in CpG context (bool)
- `mutation_type`: Mutation type (str)
- `methylation_level`: Methylation level (int32)

**gnomAD Exomes Data:**
- `gnomad_exomes_filters`: gnomAD exomes filters (set<str>)
- `gnomad_exomes_coverage`: Coverage statistics (struct with mean, median_approx, AN, percent_AN)
- `gnomad_exomes_freq`: Frequency data by population (struct with total, afr, amr, eas, nfe, sas)
- `gnomad_exomes_flags`: gnomAD exomes flags (set<str>)

**Functional Predictions:**
- `sift_score`: SIFT score (float64)
- `polyphen_score`: PolyPhen score (float64)
- `vep_domains`: VEP domains (array<struct>)
- `revel`: REVEL score (float64)
- `cadd`: CADD scores (struct with phred, raw_score)
- `phylop`: PhyloP score (float64)

**Genetics Gym Missense Scores:**
- `genetics_gym_missense_scores`: Comprehensive missense prediction scores (struct with esm_score, proteinmpnn_llr, am_pathogenicity, rasp_score, MisFit_S, MisFit_D, popeve, eve, esm1_v, mpc, and negative controls)

**Disease Associations:**
- `autism`: Autism-related annotations (struct with role)
- `dd_denovo`: Developmental delay de novo variants (struct with grch37_locus, grch37_alleles, case_control, gene flags)
- `dd_denovo_no_transcript_match`: DD de novo without transcript match
- `gnomad_de_novo`: gnomAD de novo variants (struct with de_novo_AC, p_de_novo_stats)
- `clinvar`: ClinVar annotations (comprehensive struct with rsid, frequencies, clinical significance, etc.)

**Expression Data:**
- `base_level_pext`: Base-level expression data across 49 tissues (struct with exp_prop_mean and tissue-specific data)
- `annotation_level_pext`: Annotation-level expression data (similar structure to base_level_pext)

**Constraint Metrics:**
- `mtr`: Missense tolerance ratio (struct with mtr, synExp, misExp, expMTR, synObs, misObs, obsMTR, adj_rate, pvalue, qvalue, proteinLength, percentiles)
- `rmc`: Regional missense constraint (struct with section_obs, section_exp, section_oe, section_chisq, interval, section_oe_ci, section_p_value)

##### Residue Level Annotations (`residue_level_annotations`)

**Basic Residue Info:**
- `residue_index`: Residue position in protein (int32)
- `residue_ref`: Reference amino acid (str)
- `residue_alt`: Alternative amino acid (str)

**Functional Annotations:**
- `interpro`: InterPro domain annotations (struct with interpro_id, interpro_short_description, interpro_description)
- `varity`: VariTY scores (struct with varity_r, varity_er, varity_r_loo, varity_er_loo)
- `mtr3d`: MTR3D scores (struct with mean_pLDDT, mtr3daf2_5a, mtr3daf2_8a, mtr3daf2_11a, mtr3daf2_14a)

**COSMIS Scores (Multiple Sources):**
- `cosmis_alphafold`: COSMIS scores from AlphaFold2 structures
- `cosmis_pdb`: COSMIS scores from PDB structures
- `cosmis_swiss_model`: COSMIS scores from Swiss Model structures
- `cosmis`: Combined COSMIS scores (struct with alphafold, pdb, swiss_model)

**Proemis3D Annotations:**
- `promis3d`: Proemis3D constraint annotations (struct with residue_level_annotations and region_level_annotations)

**Proemis3D Residue Level:**
- `residue_to_region_aa_dist_stats`: Amino acid distance statistics within region
- `alphafold2_info`: AlphaFold2 metrics (residue_plddt, residue_to_region_pae_stats, residue_to_region_dist_stats)

**Proemis3D Region Level:**
- `region_index`: Region index (int32)
- `region_residues`: Array of residue indices in region (array<int32>)
- `region_length`: Length of region (int32)
- `obs`: Observed missense variants (int64)
- `exp`: Expected missense variants (float64)
- `oe`: Observed/expected ratio (float64)
- `oe_upper`: Upper bound of OE confidence interval (float64)
- `oe_ci`: OE confidence interval (struct with lower, upper)
- `chisq`: Chi-squared statistic (float64)
- `p_value`: P-value (float64)
- `is_null`: Whether region is null (bool)
- `region_aa_dist_stats`: Region amino acid distance statistics
- `alphafold2_info`: Region-level AlphaFold2 metrics (region_plddt, region_pae_stats, region_dist_stats)

##### Gene Level Annotations (`gene_level_annotations`)

**Basic Gene Info:**
- `strand`: Gene strand (str)
- `cds_length`: CDS length (int32)
- `cds_len_mismatch`: CDS length mismatch flag (bool)
- `cds_len_not_div_by_3`: CDS length not divisible by 3 flag (bool)
- `aminoacid_length`: Amino acid length (int32)
- `gene`: Gene symbol (str)
- `canonical`: Canonical transcript flag (bool)
- `mane_select`: MANE select flag (bool)

**Constraint Metrics (syn, mis, lof):**
Each variant type (synonymous, missense, loss-of-function) includes:
- `mu_snp`: Mutation rate (float64)
- `mu`: Mutation rate (float64)
- `possible_variants`: Number of possible variants (int64)
- `coverage_correction`: Coverage correction factor (float64)
- `flags`: Quality flags (set<str>)
- `z_score`: Z-score (float64)
- `upper_rank`: Upper rank (int64)
- `upper_bin_sextile`: Upper bin sextile (int32)
- `upper_bin_decile`: Upper bin decile (int32)
- `observed_variants`: Number of observed variants (int64)
- `predicted_proportion_observed`: Predicted proportion observed (float64)
- `expected_variants`: Number of expected variants (float64)
- `oe`: Observed/expected ratio (float64)
- `oe_ci`: OE confidence interval (struct with lower, upper)
- `z_raw`: Raw Z-score (float64)

#### Key Structure
The table is keyed by: `['locus', 'alleles', 'transcript_id', 'uniprot_id', 'gene_id']`

This comprehensive schema provides detailed annotations at multiple levels (variant, residue, gene) with extensive functional predictions, constraint metrics, and structural information from AlphaFold2.
