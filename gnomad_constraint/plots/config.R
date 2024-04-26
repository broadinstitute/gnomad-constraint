library(dplyr)
library(ggplot2)
library(grid)
library(tidyr)
library(purrr, include.only = c("map_df"))
library(readr, include.only = c("read_tsv", "cols", "col_double"))
library(forcats, include.only = c("fct_recode", "fct_relevel"))
library(stringr, include.only = c("str_sub"))
library(pROC, include.only = c("roc")) # nolint
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("%>%", "dplyr")

data_versions <- c("v2.1.1", "v4.0")

####################################################################
# Define colors for plots
####################################################################
color_lof <- "#9D1309"
color_dominant <- "#F0810F"
color_recessive <- "#FFBB00"
color_benign <- "#87A4DC"
gene_list_colors <- c(
  "Haploinsufficient" = color_lof,
  "Essential Genes" = "#CB0000",
  "Autosomal Dominant" = color_dominant,
  "Autosomal Recessive" = color_recessive,
  "Olfactory Genes" = color_benign,
  "Background" = "lightgray",
  "Universe" = "lightgray"
)

####################################################################
# Constants and functions related to data loading
####################################################################
default_output_path <- "."
default_data_dir <- "data"
default_plot_dir <- "plots"
default_gene_lists <- "data/gene_lists"
bucket <- "gs://gnomad"

setup_directories <- function(
    wd_path,
    output_path = default_output_path,
    data_dir = default_data_dir,
    plot_dir = default_plot_dir,
    gene_lists = default_gene_lists) {
  setwd(wd_path)
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, data_dir)))
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, plot_dir)))
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, gene_lists)))
}

get_data_url <- function(version = "v4.0") {
  # TODO: Replace with sprintf("gs://gnomad/%s/constraint/", version)
  if (version == "v2.1.1") {
    data_path <-
      sprintf("gs://gcp-public-data--gnomad/release/2.1.1/constraint/")
  } else if (version == "v4.0") {
    data_path <- "gs://gnomad-kristen/constraint/"
  }

  return(data_path)
}

get_or_download_file <- function(
    base_fname,
    version = "v4.0",
    output_path = default_output_path,
    data_dir = default_data_dir,
    subfolder = "",
    local_name = "") {
  local_path <- paste0(
    output_path,
    "/",
    data_dir,
    "/",
    ifelse(local_name != "", local_name, base_fname)
  )
  print(local_path)
  if (version == "v2.1.1") {
    data_bucket <- "gs://gcp-public-data--gnomad"
  } else {
    data_bucket <- "gs://gnomad-kristen"
  }
  if (!file.exists(local_path)) {
    remote_path <- paste0(get_data_url(version), subfolder, base_fname)
    print(paste0("Downloading ", remote_path, " to ", local_path))
    print(strsplit(remote_path, paste0(data_bucket, "/"), fixed = TRUE)[[1]][2])
    gcs_get_object(
      strsplit(
        remote_path,
        paste0(data_bucket, "/"),
        fixed = TRUE
      )[[1]][2],
      data_bucket,
      saveToDisk = local_path
    )
  }

  return(local_path)
}

get_plot_path <- function(
    base_fname,
    version = "",
    output_path = default_output_path,
    plot_dir = default_plot_dir,
    extension = ".png") {
  if (version == "") {
    format_path <- "%s/%s/%s%s%s"
  } else {
    format_path <- "%s/%s/%s.%s%s"
  }
  return(sprintf(format_path, output_path, plot_dir, base_fname, version, extension))
}

load_constraint_metrics <- function(
    version = "v4.0",
    output_path = default_output_path,
    data_dir = default_data_dir) {
  # TODO: Change to sprintf('gnomad.%s.lof_metrics.by_gene.txt.bgz', version) when ready
  if (version == "v2.1.1") {
    data_path <- "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
  } else if (version == "v4.0") {
    data_path <- "extra_annotations.tsv"
  }
  fname <- get_or_download_file(data_path, version, output_path, data_dir)

  return(read.delim(fname))
}

load_all_gene_list_data <- function(
    output_path = default_output_path,
    gene_lists = default_gene_lists) {
  list_dir <- sprintf("%s/%s/lists", output_path, gene_lists)
  all_files <- list.files(list_dir, ".+tsv")
  if (length(all_files) == 0) {
    if (length(all_files) == 0) {
      system(
        paste0("git clone https://github.com/macarthur-lab/gene_lists.git ", list_dir)
      )
    }
    all_files <- list.files(list_dir, ".+tsv")
  }
  gene_lists <- map_df(
    all_files,
    function(x) {
      read_tsv(
        paste(list_dir, x, sep = "/"),
        col_names = F,
        col_types = cols()
      ) %>% transmute(gene = toupper(X1), gene_list = str_sub(x, end = -5))
    }
  )

  # Define olfactory genes based on gene names
  or_genes <- constraint_data %>%
    filter(grepl("^OR", gene)) %>%
    transmute(gene = gene, gene_list = "Olfactory Genes")

  # Add olfactory genes to gene list
  gene_lists <- gene_lists %>% bind_rows(or_genes)

  # Rename gene lists
  gene_lists <- gene_lists %>%
    mutate(
      gene_list = if_else(
        grepl("haploinsufficiency", gene_list),
        "Haploinsufficient",
        gene_list
      )
    ) %>%
    mutate(
      gene_list = fct_recode(
        gene_list,
        "Autosomal Recessive" = "all_ar",
        "Autosomal Dominant" = "all_ad"
      )
    )

  return(gene_lists)
}
