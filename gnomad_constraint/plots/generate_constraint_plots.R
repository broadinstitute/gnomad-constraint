# Example use:
# Rscript gnomad_constraint/plots/generate_constraint_plots.R \
#  -o "/Users/my_name/Documents/Constraint/plots/" \
#  -w "/Users/my_name/git_projects_location/gnomad-constraint/gnomad_constraint/plots/"
#  -t "/path/to/token.rds"
#
# To set up authentication for Google Cloud Storage (GCS) access, go to
# https://console.cloud.google.com/apis/credentials and download the JSON file for
# the tidyverse-gcs-access under OAuth 2.0 Client IDs. Then pass that file to the
# interactive_authenticate_gcs function in an interactive session. The saved token
# can be passed to this script.

library(optparse)
library(glue)

# Define input arguments
option_list <- list(
  make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NA,
    help = "Output directory to save plots to."
  ),
  make_option(
    c("-w", "--working_directory"),
    type = "character",
    default = NA,
    help = "Working directory where scripts are located locally."
  ),
  make_option(
    c("-t", "--gcs_auth_token"),
    type = "character",
    default = NA,
    help = "Path to the RDS file containing the GCS authentication token."
  )
)

# Parse the options
opt <- parse_args(OptionParser(option_list = option_list))
output_basedir <- opt$output_directory

# Set the working directory if supplied by user
if (!is.na(opt$working_directory)) {
  setwd(opt$working_directory)
}
print(glue("Working directory is {getwd()}"))

# Source the loading and plotting functions R files
source("loading_functions.R")
source("plotting_functions.R")

# Setup directory structure for output
setup_directories(output_basedir)

# Authenticate google storage if client secrets JSON is supplied
if (!is.na(opt$gcs_auth_token)) {
  authenticate_gcs(opt$gcs_auth_token)
}

####################################################################
####################################################################
# Load and join constraint data from v2 and v4
####################################################################
####################################################################
v2_constraint_data <- load_constraint_metrics(
  version = "v2",
  output_basedir = output_basedir
)

v4_constraint_data <- load_constraint_metrics(
  version = "v4",
  output_basedir = output_basedir,
)

v2_constraint_data <- mutate(v2_constraint_data, v2 = TRUE)
v4_constraint_data <- mutate(v4_constraint_data, v4 = TRUE)

constraint_data <- full_join(
  v2_constraint_data,
  v4_constraint_data,
  by = c("gene", "transcript")
)

v4 <- filter(constraint_data, .data$v4)
v2 <- filter(constraint_data, .data$v2)

####################################################################
####################################################################
# Load in gene lists
####################################################################
####################################################################
gene_lists <- load_all_gene_list_data(output_basedir = output_basedir)

# Define olfactory genes based on gene names
or_genes <- constraint_data %>%
  filter(grepl("^OR", gene)) %>%
  transmute(gene = gene, gene_list = as.factor("Olfactory Genes"))

# Add olfactory genes to gene list
gene_lists <- gene_lists %>% bind_rows(or_genes)

# Obtain just haploinsufficient genes
hi_genes <- gene_lists %>% filter(gene_list == "Haploinsufficient")

# Join gene lists to constraint data
gene_data <- left_join(constraint_data, gene_lists, by = "gene")

####################################################################
####################################################################
# Plot gene list distribution
####################################################################
####################################################################
versions_to_plot <- c("v2", "v4")
lof_upper_bin <- list(
  v2 = "oe_lof_upper_bin",
  v4 = "lof.oe_ci.upper_bin_decile"
)

for (version in versions_to_plot) {
  gene_list_sums <- summarize_gene_lists(gene_data, lof_upper_bin[[version]], version)
  summary_gene_list_per_sums <- gene_list_sums %>% spread(.data$gene_list, .data$count)

  # Write out table of gene list membership
  txt_path <- get_plot_path(
    "gene_list_counts",
    version = version,
    output_basedir = output_basedir,
    extension = ".txt"
  )
  write.table(summary_gene_list_per_sums, file = txt_path, quote = FALSE)

  # Plot gene list distribution
  p <- plot_gene_lists(gene_list_sums, lof_upper_bin[version])
  plot_path <- get_plot_path(
    "gene_list_barplot",
    version = version,
    output_basedir = output_basedir
  )
  ggsave(p, filename = plot_path, dpi = 300, width = 11, height = 6, units = "in")
}

####################################################################
####################################################################
# Plot ROC Curves
####################################################################
####################################################################
metric_by_version <- list(
  loeuf = list(
    v2 = "oe_lof_upper",
    v4 = "lof.oe_ci.upper",
    title_label = "LOEUF"
  ),
  pli = list(
    v2 = "pLI",
    v4 = "lof.pLI",
    title_label = "pLI"
  )
)
for (metric in names(metric_by_version)) {
  v2_metric <- metric_by_version[[metric]]$v2
  v4_metric <- metric_by_version[[metric]]$v4
  metric_title <- metric_by_version[[metric]]$title_label

  # Filter to where the metric is defined in both v2 and v4
  roc_df <- constraint_data %>%
    filter(!is.na(!!sym(v2_metric)) & !is.na(!!sym(v4_metric)))
  v2_roc <- plot_roc(roc_df, hi_genes, v2_metric)
  v4_roc <- plot_roc(roc_df, hi_genes, v4_metric)

  # Get combine ROC curve plot
  roc_plot <- combine_roc_plots(v2_roc, v4_roc, "v2", "v4", metric_title)
  plot_path <- get_plot_path(glue("roc_plot_{metric}"), output_basedir = output_basedir)
  ggsave(roc_plot, filename = plot_path, dpi = 300, width = 6, height = 6, units = "in")
}

####################################################################
####################################################################
# Plot downsampling projections
####################################################################
####################################################################
options(scipen = 50)

# Load in downsampling data for v2 and v4
v2_ds <- load_constraint_metrics(
  version = "v2",
  output_basedir = output_basedir,
  downsamplings = TRUE
)
v4_ds <- load_constraint_metrics(
  version = "v4",
  output_basedir = output_basedir,
  downsamplings = TRUE,
  release = FALSE,
  public = FALSE,
) %>%
  rename(
    exp_syn = syn.exp,
    exp_mis = mis.exp,
    exp_lof = lof.exp,
    obs_syn = syn.obs,
    obs_mis = mis.obs,
    obs_lof = lof.obs
  )

# Filter to canoncial/MANE Select transcripts and fit linear models
v2_ds <- filter(
  v2_ds,
  (.data$canonical == "true") &
    (.data$pop == "global") &
    (.data$downsampling >= 100)
)
v4_ds <- filter(
  v4_ds,
  (.data$mane_select == "true") &
    grepl("^ENST", .data$transcript) &
    (.data$gen_anc == "global") &
    # Note ">" rather than ">=" difference from v2 to avoid dropping too many rows
    (.data$downsampling > 100) &
    # Remove 8 genes with missing gene names
    (!is.na(.data$gene))
)


# Fit linear models for lof, mis, and syn
# Define datasets and versions
datasets <- list(v2 = v2_ds, v4 = v4_ds)

# Iterate over datasets
for (version in names(datasets)) {
  df <- datasets[[version]]

  # Generate plots
  plots <- plot_projected_sample_size(df)

  for (var_type in c("lof", "mis")) {
    p <- plots[[var_type]]

    # Construct output path and save plot
    plot_path <- get_plot_path(
      glue("{var_type}_ds_projections_{version}"),
      output_basedir = output_basedir
    )
    ggsave(p, filename = plot_path, dpi = 300, width = 8, height = 6, units = "in")
  }
}
