library(dplyr)
library(ggplot2)
library(grid)
library(tidyr)
library(purrr, include.only = c("map_df"))
library(readr, include.only = c("read_tsv", "cols", "col_double"))
library(forcats, include.only = c("fct_recode", "fct_relevel"))
library(stringr, include.only = c("str_sub"))
library(pROC, include.only = c("roc"))
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

setup_directories <-
  function(wd_path,
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

get_or_download_file <-
  function(base_fname,
           version = "v4.0",
           output_path = default_output_path,
           data_dir = default_data_dir,
           subfolder = "",
           local_name = "") {
    local_path <-
      paste0(output_path,
             "/",
             data_dir,
             "/",
             ifelse(local_name != "", local_name, base_fname))
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
      gcs_get_object(strsplit(remote_path, paste0(data_bucket, "/"), fixed = TRUE)[[1]][2],
                     data_bucket,
                     saveToDisk = local_path)
    }

    return(local_path)
  }

get_plot_path <-
  function(base_fname,
           version = "",
           output_path = default_output_path,
           plot_dir = default_plot_dir,
           extension = ".png") {
    if (version == "") {
      format_path <- "%s/%s/%s%s%s"
    } else {
      format_path <- "%s/%s/%s.%s%s"
    }
    return(sprintf(
      format_path,
      output_path,
      plot_dir,
      base_fname,
      version,
      extension
    ))
  }

load_constraint_metrics <-
  function(version = "v4.0",
           output_path = default_output_path,
           data_dir = default_data_dir) {
    # TODO: Change to sprintf('gnomad.%s.lof_metrics.by_gene.txt.bgz', version) when ready
    if (version == "v2.1.1") {
      data_path <- "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    } else if (version == "v4.0") {
      data_path <- "extra_annotations.tsv"
    }
    fname <-
      get_or_download_file(data_path, version, output_path, data_dir)
    return(read.delim(fname))
  }

load_all_gene_list_data <-
  function(output_path = default_output_path,
           gene_lists = default_gene_lists) {
    list_dir <- sprintf("%s/%s/lists", output_path, gene_lists)
    all_files <- list.files(list_dir, ".+tsv")
    if (length(all_files) == 0) {
      if (length(all_files) == 0) {
        system(
          paste0(
            "git clone https://github.com/macarthur-lab/gene_lists.git ",
            list_dir
          )
        )
      }
      all_files <- list.files(list_dir, ".+tsv")
    }
    gene_lists <- map_df(all_files, function(x) {
      read_tsv(paste(list_dir, x, sep = "/"),
               col_names = F,
               col_types = cols()) %>% transmute(gene = toupper(X1), gene_list = str_sub(x, end = -5))
    })

    # Define olfactory genes based on gene names
    or_genes <- constraint_data %>%
      filter(grepl("^OR", gene)) %>%
      transmute(gene = gene, gene_list = "Olfactory Genes")

    # Add olfactory genes to gene list
    gene_lists <- gene_lists %>% bind_rows(or_genes)

    # Rename gene lists
    gene_lists <- gene_lists %>%
      mutate(gene_list = if_else(
        grepl("haploinsufficiency", gene_list),
        "Haploinsufficient",
        gene_list
      )) %>%
      mutate(
        gene_list = fct_recode(
          gene_list,
          "Autosomal Recessive" = "all_ar",
          "Autosomal Dominant" = "all_ad"
        )
      )

    return(gene_lists)
  }

plot_projected_sample_size <- function(df, version) {
  # Get plot of the percent of genes that have a variable expected number of variants
  # across sample sizes
  # df: dataframe consisting of downsampling data per specific genetic ancestry groups
  # (has to include 'global')
  # version: version of gnomAD to use for the plot
  # Returns: plot of the percent of genes that would be expected to have a certain
  # number or variants based on sample size

  # Filter to canoncial/MANE Select transcripts and fit linear models
  if (version == "v2") {
    df <-
      filter(df,
             (canonical == "true") & (pop == "global") & (downsampling >= 100))
  } else {
    df <-
      filter(
        df,
        (mane_select == "true") &
          (grepl("^ENST", transcript)) &
          (gen_anc == "global") & (downsampling > 100) & (!is.na(gene))
      )
  } # Remove 8 genes with missing gene names

  # Filter to rows where a respective gene has at least 1 lof, mis, and syn variant within all its respective rows
  # Convert expected counts and n downsamplings to log scale
  # Fit linear models where expected counts are a function of downsampling size for each gene
  df <- df %>%
    group_by(gene) %>%
    filter(min(exp_lof) > 0 &
             min(exp_mis) > 0 & min(exp_syn) > 0) %>%
    mutate(
      log_exp_lof = log10(exp_lof),
      log_exp_mis = log10(exp_mis),
      log_exp_syn = log10(exp_syn),
      log_n = log10(downsampling)
    ) %>%
    summarize(
      lof_fit = list(lm(log_exp_lof ~ log_n)),
      mis_fit = list(lm(log_exp_mis ~ log_n)),
      syn_fit = list(lm(log_exp_syn ~ log_n))
    )

  # Extract slope and intercept for lof variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
  # why need this step to sum?
  gene_lof_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(lof_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(lof_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(gene) %>%
    summarize(slope = sum(slope), intercept = sum(intercept))

  # Extract slope and intercept for missense variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
  gene_mis_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(mis_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(mis_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(gene) %>%
    summarize(slope = sum(slope), intercept = sum(intercept))

  ####################################################################
  # Process predictions
  ####################################################################
  post_process_predictions <- function(data) {
    # Transform predictions at specific points ('5', '10', '20', '50', '100') using fitted model coefficients
    # 'data' is the dataframe to use (should contain columns 'gene', 'slope' and 'intercept')
    # Calculates the percentile rank of each value in the n_required column within the respective n_variants group
    data %>%
      mutate(
        `5` = 10 ^ ((log10(5) - intercept) / slope),
        `10` = 10 ^ ((log10(10) - intercept) / slope),
        `20` = 10 ^ ((log10(20) - intercept) / slope),
        `50` = 10 ^ ((log10(50) - intercept) / slope),
        `100` = 10 ^ ((log10(100) - intercept) / slope)
      ) %>%
      select(-slope,-intercept) %>%
      gather("n_variants", "n_required",-gene) %>%
      mutate(n_variants = as.numeric(n_variants)) %>%
      group_by(n_variants) %>%
      mutate(rank = percent_rank(n_required)) %>%
      ungroup()
  }

  samples_required_lof <- post_process_predictions(gene_lof_fit)
  samples_required_mis <- post_process_predictions(gene_mis_fit)

  ####################################################################
  # Plot projections
  ####################################################################
  # Define the limits of the x-axis
  xlimits <- c(100, 1e8)

  # Create dataframe of ExAC, gnomAD v2, and gnomAD v4 sample sizes
  lines <-
    data.frame(intercepts = c(60706, 141456, 807162),
               names = c("1", "2", "3"))

  expected_projections <- function(projection_df, label = "pLoF") {
    # Generate plot displaying the percent of genes that would be expected to have a certain number or variants based on sample size
    # 'projection_df is the input dataframe with columns 'n_variants', 'n_required', and 'rank'
    # 'label' is the text that will be included at the top of the plot
    p <- projection_df %>%
      mutate(n_variants = forcats::fct_reorder(as.factor(n_variants), n_variants)) %>%
      ggplot() +
      aes(y = rank, x = n_required, color = n_variants) +
      geom_line(size = 2) +
      theme_classic() +
      scale_x_log10(labels = scales::comma, limits = xlimits) +
      xlab("Sample size required") +
      ylab("Percent of human genes") +
      geom_vline(xintercept = 60706, linetype = "dotted") + # Exac size
      geom_vline(xintercept = 141456, linetype = "dashed") + # gnomAD v2 size
      geom_vline(xintercept = 807162, linetype = "twodash") + # gnomAD v4 size

      # Add manual linetype scale for legend
      p <-
        p + scale_linetype_manual(
          name = "Database size",
          values = c("dotted", "dashed", "twodash"),
          labels = c("ExAC", "gnomADv2", "gnomADv4")
        ) + geom_vline(
          data = lines,
          aes(xintercept = intercepts, linetype = names),
          key_glyph = "path"
        ) + scale_color_discrete(name = ">= N variants\nexpected", h = c(40, 120)) + annotate(
          "text",
          x = xlimits[1],
          y = 1,
          hjust = 0,
          vjust = 1,
          label = label
        )

      return(p)
  }

  lof_projections <-
    expected_projections(samples_required_lof, "pLoF")
  missense_projections <-
    expected_projections(samples_required_mis, "Missense")

  ggsave(
    lof_projections,
    filename = paste("lof_ds_projections_", version, ".png", sep = ""),
    dpi = 300,
    width = 8,
    height = 6,
    units = "in"
  )
  ggsave(
    missense_projections,
    filename = paste("mis_ds_projections_", version, ".png", sep = ""),
    dpi = 300,
    width = 8,
    height = 6,
    units = "in"
  )
}
