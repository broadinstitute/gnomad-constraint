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


output_basedir = "/Users/kristen/Desktop/constraint/plots"
working_directory = "/Users/kristen/Desktop/repos/gnomad-constraint/gnomad_constraint/plots"

setwd(working_directory)

# Source the loading and plotting functions R files
source("loading_functions.R")
source("plotting_functions.R")

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

v2_constraint_data <- mutate(v2_constraint_data, v2)
v4_constraint_data <- mutate(v4_constraint_data, v4)

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

test = filter(gene_data, gene== "CHD2")
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

preferred_transcripts <- list(
  v2 = "canonical",
  v4 = "mane_select"
)


filter_transcripts <- function(
    df,
    version,
    preferred_transcripts_only = TRUE,
    enst_only = TRUE,
    include_flagged_transcripts = FALSE) {
  # Filter transcripts to specified criteria
  # df: Dataframe on which filters should be applied
  # version: Version of the dataset, either 'v2' or 'v4. Dataframe will be filtered to
  # variants in this version.
  # preferred_transcripts_only: Whether to filter to only the preferred transcript.
  # Default is TRUE. For 'v4', 'mane_select' are the preferred transcripts. For 'v2',
  # 'canonical' transcripts are the preferred transcripts.
  # enst_only: Whether to filter to only Ensembl transcripts. Default is TRUE.
  # include_flagged_transcripts: Whether to keep transcripts with a constraint flag.
  # Default is FALSE.
  # Returns: Filtered dataframe

  # Filter to transcripts present in specified version
  df <- filter(df, !!sym(version))

  # Filter to only Ensembl transcripts if specified
  if (enst_only) {
    df <- filter(df, grepl("^ENST", .data$transcript))
  }

  # Filter to exclude flagged transcripts if specified
  if (!include_flagged_transcripts) {
    df <- filter(df, .data[[paste0("constraint_flags.", version)]] == "[]")
  }

  # Filter to only preferred transcripts if specified
  if (preferred_transcripts_only) {
    df <- filter(df, !!sym(preferred_transcripts[[version]]) == TRUE)
  }

  return(df)
}


source("plotting_functions.R")

for (version in versions_to_plot) {
  # Filter to preferred  Ensembl transcripts (mane select or canonical) without any constraint flags
  version = "v4"
  filtered_data <- filter_transcripts(gene_data, version)

  #H1
  filtered_data <- rename(filtered_data, lof.oe_ci.upper_bin_decile.v4 = lof.oe_ci.upper_bin_decile)
  #H1
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


constraint_data <- rename(constraint_data, lof.oe_ci.upper.v4 = lof.oe_ci.upper,
                          lof.pLI.v4 = lof.pLI,
                          lof.oe_ci.upper.v2 = oe_lof_upper,
                          lof.pLI.v2 = pLI)


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


v4_ds <- load_constraint_metrics(
  version = "v4",
  output_basedir = output_basedir,
  downsamplings = TRUE,
  release = FALSE,
  public = FALSE,
) 

#H1
v2_ds <- v2_ds %>%
  rename(
    syn.exp = exp_syn,
    mis.exp = exp_mis,
    lof.exp = exp_lof,
    syn.obs = obs_syn,
    mis.obs = obs_mis,
    lof.obs =obs_lof
  )
#H1

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


source("plotting_functions.R")

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

source("plotting_functions.R")

decile_data <- rename(decile_data, lof.oe_ci.upper_bin_decile.v2 = oe_lof_upper_bin)

colnames(constraint_data)
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
comparision_df <- filter_transcripts(constraint_data, "v4") %>% filter(.data$v2 == TRUE)


comparision_df <- rename(comparision_df, lof.oe_ci.upper.v2 = oe_lof_upper,
                         lof.oe_ci.upper.v4 = lof.oe_ci.upper,
                         lof.z_score.v4 = lof.z_score,
                         lof.z_score.v2 = lof_z)



comparision_df <- comparision_df %>%
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

#comparision_df <- comparision_df %>%
#  select(
#    all_of(
#      c(
#        "transcript",
#        "v2",
#        "v4",
#        "oe_lof_upper",
#        "lof.oe_ci.upper",
#        "lof_z",
#        "lof.z_score"
#      )
    )
  )

# Pivot to longer df with metric name and value
comparision_df <- comparision_df %>%
  pivot_longer(
    cols = c(lof.oe_ci.upper.v2, lof.oe_ci.upper.v4, lof.z_score.v2, lof.z_score.v4),
    names_to = "metric_name",
    values_to = "value"
  )

# Pull out version
comparision_df <- comparision_df %>% mutate(
  version = if_else(grepl("v4", metric_name), "v4", "v2")
)

# Rename metrics in df to be the same for v2 and v4
comparision_df$metric_name <- gsub("\\.v(2|4)$", "", comparision_df$metric_name)

comparision_df <- comparision_df %>%
  mutate(metric_name = case_when(
    .data$metric_name == "lof.oe_ci.upper" ~ "LOEUF",
    .data$metric_name == "lof.z_score" ~ "pLoF Z-score",
    TRUE ~ .data$metric_name # Keeps all other values as they are
  ))

comparision_df <- comparision_df %>% select(-v2, -v4)

# Pivot to wider format
comparision_df <- comparision_df %>% pivot_wider(
  names_from = version, values_from = value
)

# Calculate correlations for each metric between v2 and v4
comparision_df %>%
  group_by(metric_name) %>%
  summarize(
    correlation = cor(.data$v2, .data$v4, method = "pearson", use = "complete.obs")
  )

comparsion_plot <- plot_metric_comparison(comparision_df, use_presentation_sizes=use_presentation_sizes)
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

constraint_data <- rename(constraint_data, 
                         syn.exp.v2 = exp_syn,
                         mis.exp.v2 = exp_mis,
                         lof.exp.v2 = exp_lof,
                         syn.obs.v2 = obs_syn,
                         mis.obs.v2 = obs_mis,
                         lof.obs.v2 = obs_lof,
                         syn.exp.v4 = syn.exp,
                         syn.obs.v4 = syn.obs,
                         mis.exp.v4 = mis.exp,
                         mis.obs.v4 = mis.obs,
                         lof.exp.v4 = lof.exp,
                         lof.obs.v4 = lof.obs)


for (version in versions_to_plot) {
  # Filter to preferred Ensembl transcripts (mane select or canonical) without any constraint flags
  filtered_data <- filter_transcripts(
    constraint_data,
    version,
    preferred_transcripts_only = TRUE,
    enst_only = TRUE,
    include_flagged_transcripts = FALSE
  )
  # Reformat column names 
    oe_data <- filtered_data %>% select(
        c(
          "transcript",
          !!sym(version),
          glue("syn.exp.{version}"),
          glue("syn.obs.{version}"),
          glue("mis.exp.{version}"),
          glue("mis.obs.{version}"),
          glue("lof.exp.{version}"),
          glue("lof.obs.{version}")
        )
    )  %>% rename_with(~ gsub("\\.v(4|2)$", "", .))

  # Reformat the data
  oe_data <- oe_data %>%
    pivot_longer(
      # Keep 'transcript' and version, pivot all other columns
      cols = -c(transcript, version),
      names_to = "metric_combo",
      values_to = "counts"
    ) %>%
    # Split metric_combo ("syn.exp") into metric_name ("syn") and count_type ("exp")
    separate(metric_combo, into = c("metric_name", "count_type"), sep = "\\.") %>%
    pivot_wider(
      names_from = count_type,
      values_from = counts,
      names_prefix = ""
    ) %>%
    select(-any_of(c("v2", "v4")))

  # Rname lof to pLoF and reorder the values in the df
  oe_data <- oe_data %>%
    mutate(metric_name = recode(.data$metric_name, "lof" = "pLoF")) %>%
    mutate(metric_name = factor(.data$metric_name, levels = c("syn", "mis", "pLoF")))

  # Find the max value of the observed and expected values for each metric
  max_values <- oe_data %>%
    group_by(metric_name) %>%
    summarise(
      max_exp = max(.data$exp, na.rm = TRUE),
      max_obs = max(.data$obs, na.rm = TRUE)
    ) %>%
    mutate(max_limit = pmax(max_exp, max_obs))

  # Join these max values back to the original data to use in plotting
  oe_data <- oe_data %>%
    left_join(max_values, by = "metric_name")

  oe_plot <- plot_observed_vs_expected(oe_data, version)

  plot_path <- get_plot_path(
    "obs_vs_exp",
    version = version,
    output_basedir = output_basedir
  )
  ggsave(oe_plot, filename = plot_path, dpi = 300, width = 7, height = 6, units = "in")
}
