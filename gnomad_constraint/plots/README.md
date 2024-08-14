## Install R
Before making plots, make sure you have R installed on your machine. You can download
and install R from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

## Install `renv`
This project uses `renv` for dependency management. `renv` helps manage
project-specific R dependencies, ensuring that the R packages are consistent across all
project collaborators. To install `renv`, run the following command in your R console:
```R
if (!require(renv)) install.packages("renv")
```

## Restore Dependencies
Once `renv` is installed, restore the project dependencies to ensure you're using the
same package versions as the rest of the team:
```R
renv::restore()
```
This command should be run with the working directory in R set to the directory that contains the 'renv.lock' file. For more detailed information on `renv`, please refer to the official `renv`
documentation: [Using `renv` with R](https://rstudio.github.io/renv/articles/renv.html).

## Generating Constraint Plots
### Set-up Google Cloud Storage (GCS) authentication
To set up authentication for Google Cloud Storage (GCS) access, go to
https://console.cloud.google.com/apis/credentials and download the JSON file for
the `tidyverse-gcs-access` under `OAuth 2.0 Client IDs`. Then pass that file to the
`interactive_authenticate_gcs` function in an interactive R session. The saved token
can be passed `generate_constraint_plots.R`.

If you prefer not to set up GCS authentication, just download the required input files
and store in the location you plan to use as output for `generate_constraint_plots.R`
under a `data` directory.

### Running R script
From the command line:
```commandline
Rscript gnomad_constraint/plots/generate_constraint_plots.R \
  -o "/Users/my_name/Documents/Constraint/plots/" \
  -w "/Users/my_name/git_projects_location/gnomad-constraint/gnomad_constraint/plots/"
  -t "/path/to/token.rds"
```
Running this script will create subdirectories in the output directory if they don't already exist. For example, the above command would create the following directories:
- /Users/my_name/Documents/Constraint/plots/plots
- /Users/my_name/Documents/Constraint/plots/gene_lists/lists
- /Users/my_name/Documents/Constraint/plots/data


If not using a token to download the files from the Google Cloud bucket, the necessary files should be placed in the "data" and "gene_lists/lists" subdirectories. For example, "/Users/my_name/Documents/Constraint/plots/data" would need to contain the following:

- gnomad.v4.1.constraint_metrics.tsv
- gnomad.v4.1.downsampling_constraint_metrics.tsv.bgz
- gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
- gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz

## Contributing
If you're adding new dependencies to the project, update the `renv` lockfile by running:
```R
renv::snapshot()
```

Before pushing your changes, make sure to run the pre-commit hooks to ensure your code
adheres to the project's style guidelines and passes all checks:
```commandline
pre-commit run --all-files
```

The above pre-commit hooks run [`styler`](https://styler.r-lib.org/) and
[`lintr`](https://lintr.r-lib.org/). They use the
[tidyverse style guide](https://style.tidyverse.org/) as the code format. A few things
to note about their use:
 - Unfortunately, at the time this was last updated, `styler` does not handle reformatting long lines. If you add the newlines yourself, you don't need to worry about correcting the tabs because `styler` can handle that for you.
 - `dplyr` uses non-standard evaluation to interpret the names of columns inside its verbs like `mutate` and `select`, and this causes `lintr` to issue warnings of "no visible binding for global variable" because the linter does not recognize column names in your dataframe when used within `dplyr` functions.
   - With most `dplyr` statements, this can be fixed with the use of '.data' before the column name.
   - The `select` statement requires the use of quotes around the column names where possible.
   - The `all_of()` function in the `select` statement is typically the best approach for cases where column names are dynamic or external but still need to be recognized as originating from within the data frame.
