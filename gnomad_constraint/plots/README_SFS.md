# Site Frequency Spectrum (SFS) Computation

This directory contains a Jupyter notebook for computing Site Frequency Spectrum tables from constraint pipeline data.

## Overview

The `compute_sfs_tables.ipynb` notebook generates SFS tables for each combination of:
- **Functional categories**: synonymous (syn), missense (mis), loss-of-function (lof)
- **Mutation type categories**: transversion, transition, CpG

## Output Format

For each combination, the notebook generates a compressed tab-separated file named:
```
SFS_{functional_category}_{mutation_type_category}.txt.gz
```

### Example Output Files
- `SFS_syn_transversion.txt.gz`
- `SFS_syn_transition.txt.gz`
- `SFS_syn_CpG.txt.gz`
- `SFS_mis_transversion.txt.gz`
- `SFS_mis_transition.txt.gz`
- `SFS_mis_CpG.txt.gz`
- `SFS_lof_transversion.txt.gz`
- `SFS_lof_transition.txt.gz`
- `SFS_lof_CpG.txt.gz`

### Table Structure

Each SFS table contains:
- **Rows**: Different downsampling levels (e.g., AC_ds10, AC_ds100, AC_ds500, etc.)
- **Columns**: Allele count bins (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10-100, 100-1000, 1000-10000, >10000)
- **Values**: Number of variants with that allele count in that downsampling

### Example Table Structure
```
Sample_size    0    1    2    3    4    5    6    7    8    9    10   10-100  100-1000  1000-10000  >10000
AC_ds1461892   357  455  234  123  89   67   45   34   23   18   735  17199    3783       887         397
AC_ds10        28627 148  89   45   23   12   8    5    3    2    5    4        0          0           0
AC_ds100       15234 234  156  89   56   34   23   15   9    6    12   23       2          0           0
...
```

## Usage

### Prerequisites
1. **Hail**: Make sure Hail is installed and configured
2. **Data Access**: Ensure you have access to the constraint pipeline data tables
3. **Python Dependencies**: pandas, numpy, gzip

### Running the Notebook

1. **Update Data Paths**: Modify the data paths in the notebook to point to your actual constraint data:
   ```python
   constraint_ht_path = "gs://your-bucket/constraint/constraint_metrics.ht"
   downsampling_ht_path = "gs://your-bucket/constraint/downsampling_constraint_metrics.ht"
   ```

2. **Adapt Data Extraction Logic**: The notebook provides a framework, but you'll need to adapt the `extract_sfs_from_constraint_data()` function based on your actual data structure.

3. **Run the Notebook**: Execute the cells in order to:
   - Load and investigate the data structure
   - Extract SFS data for each combination
   - Save the results as compressed files

### Key Functions

#### `get_ac_bins()`
Defines the allele count bins used for SFS computation:
- Individual bins: 0, 1, 2, ..., 10
- Range bins: 10-100, 100-1000, 1000-10000, >10000

#### `extract_sfs_from_constraint_data()`
Main function that extracts SFS data from constraint pipeline tables. This function:
1. Filters data by functional category and mutation type
2. Extracts frequency data for each downsampling level
3. Counts variants in each AC bin
4. Returns a pandas DataFrame with the SFS data

#### `create_ac_bin_expression()`
Creates a Hail expression to assign variants to AC bins based on their allele count.

## Data Structure Requirements

The notebook expects your constraint data to have the following structure:

### Constraint Table
- `constraint_group_meta`: Metadata about constraint groups
- `constraint_groups`: Array of constraint group data
- `mutation_type`: Mutation type annotation (transversion, non-CpG transition, CpG)

### Downsampling Table
- `annotation`: Functional annotation (synonymous_variant, missense_variant, etc.)
- `mutation_type`: Mutation type annotation
- `most_severe_consequence`: Most severe consequence annotation
- Frequency data for different downsampling levels

## Customization

### Adding New Functional Categories
To add new functional categories, modify the `get_functional_categories()` function and update the filtering logic in `extract_sfs_from_constraint_data()`.

### Adding New Mutation Types
To add new mutation types, modify the `get_mutation_type_categories()` function and update the filtering logic.

### Modifying AC Bins
To change the allele count bins, modify the `get_ac_bins()` function and update the `create_ac_bin_expression()` function accordingly.

## Troubleshooting

### Common Issues

1. **No Data Found**: If no variants are found for a combination, check:
   - Data filtering logic
   - Data availability for that combination
   - Constraint group metadata structure

2. **Memory Issues**: For large datasets, consider:
   - Processing subsets of data
   - Using checkpointing
   - Increasing cluster resources

3. **Data Structure Mismatch**: If the notebook doesn't work with your data:
   - Investigate the data structure using the provided functions
   - Adapt the filtering and extraction logic
   - Test with a small subset first

### Data Investigation

Use the data structure investigation section in the notebook to understand your data:
```python
# Show table schemas
constraint_ht.describe()
downsampling_ht.describe()

# Show sample data
constraint_ht.show(3)
downsampling_ht.show(3)
```

## Output

The notebook generates:
1. **SFS Tables**: Compressed tab-separated files for each combination
2. **Console Output**: Progress information and data statistics
3. **Error Messages**: Detailed error information for troubleshooting

## Example Usage

```python
# Load data
constraint_ht = hl.read_table("path/to/constraint_metrics.ht")
downsampling_ht = hl.read_table("path/to/downsampling_constraint_metrics.ht")

# Extract SFS for synonymous transversions
sfs_df = extract_sfs_from_constraint_data(
    constraint_ht, 
    downsampling_ht, 
    "syn", 
    "transversion"
)

# Save result
with gzip.open("SFS_syn_transversion.txt.gz", 'wt') as f:
    sfs_df.to_csv(f, sep='\t')
```

## Support

For questions or issues with the SFS computation:
1. Check the data structure investigation section
2. Verify your data paths and access permissions
3. Test with a small subset of data first
4. Review the error messages for specific issues

