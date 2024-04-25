library("optparse")
source("config.R")
source("utils.R")

# Define input arguments
option_list <- list(
  make_option(c("-o", "--output-directory"), type = "character", default = "plots", help = "Output directory to save plots to."),
  make_option(c("-w", "--working-directory"), type = "character", help = "Working directory where scripts are located locally."),
)

# Parse the options
opt <- parse_args(OptionParser(option_list = option_list))

# Overwrite default output path
default_output_path <- opt$output_directory
setup_directories(opt$working_directory, default_output_path)

####################################################################
####################################################################
# Load and join constraint data from v2 and v4
####################################################################
####################################################################
v2_constraint_data <- load_constraint_metrics(version = "v2.1.1", output_path = default_output_path)
# TODO: This is the gnomad v4.0 file, will replace when v4.1 is ready
v4_constraint_data <- load_constraint_metrics(version = "v4.0", output_path = default_output_path)

v2_constraint_data <- mutate(v2_constraint_data, v2 = TRUE)
v4_constraint_data <- mutate(v4_constraint_data, v4 = TRUE)

constraint_data <- full_join(v2_constraint_data, v4_constraint_data, by = c("gene", "transcript"))

v4 <- filter(constraint_data, v4)
v2 <- filter(constraint_data, v2)

####################################################################
####################################################################
# Load in gene lists
####################################################################
####################################################################
gene_lists <- load_all_gene_list_data(output_path = default_output_path)

# Obtain just haploinsufficient genes
hi_genes <- gene_lists %>% filter(gene_list == "Haploinsufficient")

# Join gene lists to constraint data
gene_data <- left_join(constraint_data, gene_lists, by = "gene")

####################################################################
####################################################################
# Plot gene list distribution
####################################################################
####################################################################
gene_lists_to_plot <- c("Haploinsufficient", "Autosomal Recessive", "Olfactory Genes")
versions_to_plot <- c("v2", "v4")

for (version in versions_to_plot) {
  plot_gene_lists(
    gene_data,
    version,
    get_plot_path(
      "gene_list_barplot",
      version = version,
      output_path = default_output_path
    )
  )
}

####################################################################
####################################################################
# Plot ROC Curves
####################################################################
####################################################################
plot_roc(constraint_data,"loeuf", get_plot_path(paste("roc_plot_", "loeuf"), output_path = default_output_path))
plot_roc(constraint_data,
         "pli",
         get_plot_path(paste("roc_plot_", "pli"), output_path = default_output_path))

####################################################################
####################################################################
# Plot downsampling projections
####################################################################
####################################################################
# Load in downsampling data for v2 and v4
options(scipen = 50)
v2_ds <- read.delim("gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz")
v4_ds <- read.delim("gnomad.v4.1.downsampling_constraint_metrics.txt.bgz")

# Rename v4 columns
v4_ds <-
  v4_ds %>% rename(
    exp_syn = syn.exp,
    exp_mis = mis.exp,
    exp_lof = lof.exp,
    obs_syn = syn.obs,
    obs_mis = mis.obs,
    obs_lof = lof.obs
  )

####################################################################
####################################################################
# Fit linear models for lof, mis, and syn
####################################################################
####################################################################
plot_projected_sample_size(v2_ds, "v2")
plot_projected_sample_size(v4_ds, "v4")
