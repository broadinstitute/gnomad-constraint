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
  ),
  make_option(
    c("--use_presentation_sizes"),
    type = "logical",
    default = FALSE,
    help = "Whether to use presentation sizes (larger text sizes) when generating plots."
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

setwd(working_directory)

# Source the loading, plotting, and utils functions R files
source("loading_functions.R")
source("plotting_functions.R")
source("plotting_utils.R")

# Setup directory structure for output
setup_directories(output_basedir)

# Authenticate google storage if client secrets JSON is supplied
if (!is.na(opt$gcs_auth_token)) {
  authenticate_gcs(opt$gcs_auth_token)
}

use_presentation_sizes = opt$use_presentation_sizes

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

constraint_data <- full_join(
  v2_constraint_data,
  v4_constraint_data,
  by = c("gene", "transcript")
)

v4 <- filter(constraint_data, .data$v4)
v2 <- filter(constraint_data, .data$v2)

use_presentation_sizes = FALSE
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
gene_lists <- gene_lists %>%
  bind_rows(or_genes) %>%
  distinct()

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

preferred_transcripts <- list(
  v2 = "canonical",
  v4 = "mane_select"
)

for (version in versions_to_plot) {
  # Filter to preferred  Ensembl transcripts (mane select or canonical) without any constraint flags
  filtered_data <- filter_transcripts(gene_data, version)

  # Define metric (lof.oe_ci.upper_bin_decile + version)
  metric <- paste0("lof.oe_ci.upper_bin_decile", ".", version)
  
  gene_list_sums <- summarize_gene_lists(
    filtered_data,
    metric, 
    version
  )
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
  p <- plot_gene_lists(gene_list_sums, metric, use_presentation_sizes=use_presentation_sizes)
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
metric_titles <- list(
  "lof.oe_ci.upper" = "LOEUF",
  "lof.pLI" = "pLI")

for (metric in names(metric_titles)) {
  v2_metric <- glue("{metric}.v2")
  v4_metric <- glue("{metric}.v4")
  metric_title <- metric_titles[metric]

  # Filter to where the metric is defined in both v2 and v4
  roc_df <- constraint_data %>%
    filter(!is.na(!!sym(v2_metric)) & !is.na(!!sym(v4_metric)))
  v2_roc <- plot_roc(roc_df, hi_genes, v2_metric)
  v4_roc <- plot_roc(roc_df, hi_genes, v4_metric)

  # Get combine ROC curve plot
  roc_plot <- combine_roc_plots(v2_roc, v4_roc, "v2", "v4", metric_title, use_presentation_sizes=use_presentation_sizes)
  plot_path <- get_plot_path(glue("roc_plot_{metric_title}"), output_basedir = output_basedir)
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
)


# Filter to canoncial/MANE Select transcripts and fit linear models
v2_ds <- filter(
  v2_ds,
  (.data$canonical) &
    (.data$pop == "global") &
    (.data$downsampling >= 100)
)
v4_ds <- filter(
  v4_ds,
  (.data$mane_select) &
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
  plots <- plot_projected_sample_size(df, use_presentation_sizes=use_presentation_sizes)

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


####################################################################
####################################################################
# Plot LOEUF decile change from v2 to v4
####################################################################
####################################################################
# Filter to MANE Select Ensembl transcripts that are present in both v2 and v4
# and that have no constraint flags
decile_data <- filter_transcripts(constraint_data, "v4") %>% filter(.data$v2 == TRUE)
decile_plot <- plot_decile_change(decile_data, use_presentation_sizes=use_presentation_sizes)

plot_path <- get_plot_path(
  "decile_change",
  output_basedir = output_basedir
)
ggsave(
  decile_plot,
  filename = plot_path,
  dpi = 300,
  width = 8,
  height = 4,
  units = "in"
)


####################################################################
####################################################################
# Plot LOEUF and Z-Score comparison between v2 and v4
####################################################################
####################################################################
# Filter to MANE Select Ensembl transcripts that are present in both v2 and v4
# and that have no constraint flags
comparison_df <- filter_transcripts(constraint_data, "v4") %>% filter(.data$v2 == TRUE)

comparison_df <- comparison_df %>%
  select(
    all_of(
      c(
        "transcript",
        "v2",
        "v4",
        "lof.oe_ci.upper.v2",
        "lof.oe_ci.upper.v4",
        "lof.z_score.v2",
        "lof.z_score.v4"
      )
    )
  )

# Reformat the data for plotting
comparison_df <- reformat_for_metric_comparison(comparison_df)

comparison_plot <- plot_metric_comparison(comparison_df, use_presentation_sizes=use_presentation_sizes)
plot_path <- get_plot_path(
  "v2_v4_compare_values",
  output_basedir = output_basedir
)

ggsave(
  comparsion_plot,
  filename = plot_path,
  dpi = 300,
  width = 8,
  height = 4,
  units = "in"
)


####################################################################
####################################################################
# Plot observed vs expected values
####################################################################
####################################################################
for (version in versions_to_plot) {
  # Filter to preferred Ensembl transcripts (mane select or canonical) without any constraint flags
  filtered_data <- filter_transcripts(
    constraint_data,
    version,
    preferred_transcripts_only = TRUE,
    enst_only = TRUE,
    include_flagged_transcripts = FALSE
  )
  
  # Reformat the data for plotting
  oe_data <- reformat_for_observed_vs_expected(filtered_data, version)
  # Plot observed vs expected values
  oe_plot <- plot_observed_vs_expected(oe_data, version)

  plot_path <- get_plot_path(
    "obs_vs_exp",
    version = version,
    output_basedir = output_basedir
  )
  ggsave(oe_plot, filename = plot_path, dpi = 300, width = 7, height = 6, units = "in")
}


####################################################################
####################################################################
# Plot proportion observed vs mu
####################################################################
####################################################################
setwd("/Users/kristen/code/gnomad-constraint/gnomad_constraint/plots")
source("loading_functions.R")
source("plotting_functions.R")
source("plotting_utils.R")

setwd("/Users/kristen/Desktop/an_cov")
training_data = read.delim("an_po.txt")
name=""
high_coverage_cutoff=90

training_data = read.delim("po_an_raw.txt")
name="_raw"
high_coverage_cutoff=90

training_data = read.delim("po_an_10.txt")
name="_10"
high_coverage_cutoff=9

unique(training_data$exomes_AN_percent)


version="v4"
output_basedir="/Users/kristen/Desktop/an_cov"

data_with_predictions <- get_predicted_proportion_observed(df = training_data,
                                                           coverage_metric="exomes_AN_percent",
                                                           high_coverage_cutoff=high_coverage_cutoff)

po_v_mu <- plot_proportion_observed_vs_mu(df=data_with_predictions,
                                          coverage_metric="exomes_AN_percent",
                                          high_coverage_cutoff=high_coverage_cutoff)

plot_path <- get_plot_path(
  glue("ov_v_mu{name}"),
  version = version,
  output_basedir = output_basedir
)
ggsave(po_v_mu, filename = plot_path, dpi = 300, width = 7, height = 6, units = "in")

####################################################################
####################################################################
# Plot observed to expected ratio vs coverage metric
####################################################################
####################################################################
data_with_predictions <- get_predicted_proportion_observed(df = training_data,
                                                           coverage_metric="exomes_AN_percent",
                                                           high_coverage_cutoff=high_coverage_cutoff)

oe_v_cov <- plot_oe_vs_cov_metric(df=data_with_predictions,
                                  coverage_metric="exomes_AN_percent",
                                  high_coverage_cutoff=high_coverage_cutoff,
                                  add_best_fit=TRUE)

plot_path <- get_plot_path(
  glue("oe_v_cov{name})"),
  version = version,
  output_basedir = output_basedir
)
ggsave(oe_v_cov, filename = plot_path, dpi = 300, width = 7, height = 6, units = "in")

