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
For more detailed information on `renv`, please refer to the official `renv`
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


## Contributing
If you're adding new dependencies to the project, update the `renv` lockfile by running:
```R
renv::snapshot()
```

Before pushing your changes, make sure to run the pre-commit hooks to ensure your code adheres to the project's style guidelines and passes all checks:
```commandline
pre-commit run --all-files
```
