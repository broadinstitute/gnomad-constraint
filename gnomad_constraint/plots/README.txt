## Install R
Before making plots, make sure you have R installed on your machine. You can download and install R from The Comprehensive R Archive Network (CRAN).

## Install renv
This project uses renv for dependency management. renv helps manage project-specific R dependencies, ensuring that the R packages are consistent across all project collaborators. To install renv, run the following command in your R console:
```
if (!require(renv)) install.packages("renv")
```
## Restore Dependencies
Once renv is installed, restore the project dependencies to ensure you're using the same package versions as the rest of the team:
```
renv::restore()
```
