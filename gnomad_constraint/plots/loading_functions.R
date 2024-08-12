library(conflicted)
library(dplyr)
library(gargle)
library(glue)
library(googleAuthR)
library(googleCloudStorageR)
library(forcats, include.only = c("fct_recode"))
library(purrr, include.only = c("map_df"))
library(readr, include.only = c("read_tsv", "cols"))
library(stringr, include.only = c("str_sub"))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("%>%", "dplyr")

####################################################################
# Constants and functions related to data loading
####################################################################
data_versions <- c("v2", "v4")
map_data_version <- list(
  v2 = "2.1.1",
  v4 = "4.1"
)

default_output_basedir <- "."
default_data_subdir <- "data"
default_plot_subdir <- "plots"
default_gene_list_subdir <- "data/gene_lists"
bucket <- "gs://gnomad"

setup_directories <- function(
    output_basedir = default_output_basedir,
    data_subdir = default_data_subdir,
    plot_subdir = default_plot_subdir,
    gene_list_subdir = default_gene_list_subdir) {
  # Create the output directory and subdirectories if they don't exist
  # output_basedir: The base directory for the output
  # data_subdir: The subdirectory for data files
  # plot_subdir: The subdirectory for plot files
  # gene_list_subdir: The subdirectory for gene list files
  # Returns: None
  suppressWarnings(dir.create(output_basedir))
  suppressWarnings(dir.create(glue("{output_basedir}/{data_subdir}/")))
  suppressWarnings(dir.create(glue("{output_basedir}/{plot_subdir}/")))
  suppressWarnings(dir.create(glue("{output_basedir}/{gene_list_subdir}/")))
}

get_data_url <- function(version = "v4", release = TRUE, public = TRUE, pipeline_step="metrics") {
  # Get the URL for the data files
  # version: The version of the data
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # pipeline_step: What step of the pipeline to obtain data for. Either "metrics" or "training". Default is "metrics".
  # Returns: The URL for the data files

  full_version <- map_data_version[[version]]
  if (release && public) {
    data_path <- glue("gs://gcp-public-data--gnomad/release/{full_version}/constraint/")
  } else if (release && !public) {
    data_path <- glue("gs://gnomad/release/{full_version}/constraint/")
  } else if (!release && !public) {
    if (pipeline_step == "metrics") {
    data_path <- glue("gs://gnomad/v{full_version}/constraint/metrics/tsv/")}
    if (pipeline_step == "training") {
      data_path <- glue("gs://gnomad/v{full_version}/constraint/training_data/")}
  } else {
    stop("Invalid combination of release and public")
  }

  return(data_path)
}

authenticate_gcs <- function(token_path) {
  # Authenticate with Google Cloud Storage (GCS) using a saved token
  # token_path: The path to the saved token
  # Returns: None
  token <- readRDS(token_path)
  gar_auth(token = token)
}

interactive_authenticate_gcs <- function(
    client_secret_json_path,
    token_output_path = "token.rds") {
  # Authenticate with Google Cloud Storage (GCS) interactively and save the
  # token
  # client_secret_json_path: The path to the client secret JSON file
  # token_output_path: The path to save the token
  # Returns: None
  client <- gargle_oauth_client_from_json(
    path = client_secret_json_path,
    name = "constraint-oauth-client"
  )
  token <- credentials_user_oauth2(
    "https://www.googleapis.com/auth/devstorage.read_only",
    client
  )

  # Save the token to a file
  saveRDS(token, file = token_output_path)
}

get_or_download_file <- function(
    base_fname,
    version = "v4",
    output_basedir = default_output_basedir,
    data_subdir = default_data_subdir,
    subfolder = "",
    local_name = "",
    release = TRUE,
    public = TRUE,
    pipeline_step="metrics",
    gcs_authentication_token = NULL) {
  # Get a local file, and if it doesn't exist, try downloading it from Google
  # Cloud Storage (GCS)
  # base_fname: The base filename of the file
  # version: The version of the data
  # output_basedir: The base directory for the output
  # data_subdir: The subdirectory for data files
  # subfolder: The subfolder in the data directory
  # local_name: The local name of the file
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # pipeline_step: What step of the pipeline to obtain data for. Either "metrics" or "training". Default is "metrics".
  # gcs_authentication_token: The path to the saved token for GCS authentication
  # Returns: The path to the local file for requested data
  local_name <- ifelse(local_name != "", local_name, base_fname)
  local_path <- glue("{output_basedir}/{data_subdir}/{local_name}")

  # To set up authentication for Google Cloud Storage (GCS) access, go to
  # https://console.cloud.google.com/apis/credentials and download the JSON file for
  # the tidyverse-gcs-access under OAuth 2.0 Client IDs. Then pass that file to the
  # interactive_authenticate_gcs function in an interactive session. The saved token
  # can be passed to this function.

  if (public == TRUE) {
    data_bucket <- "gs://gcp-public-data--gnomad"
  } else {
    data_bucket <- "gs://gnomad"
  }
  if (!file.exists(local_path)) {
    if (!missing(gcs_authentication_token)) {
      authenticate_gcs(gcs_authentication_token)
    }
    print(paste0("version", version))
    print(paste0("release", release))
    print(paste0("public", public))
    print(paste0("subfolder", subfolder))
    print(paste0("base_fname", base_fname))
    print(get_data_url(version, release = release, public = public, pipeline_step=pipeline_step))
    remote_path <- paste0(
      get_data_url(version, release = release, public = public),
      subfolder,
      base_fname
    )
    print(glue("Downloading {remote_path} to {local_path}"))
    print(strsplit(
      remote_path,
      paste0(data_bucket, "/"),
      fixed = TRUE
    )[[1]][2])
    print(data_bucket)
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
    output_basedir = default_output_basedir,
    plot_subdir = default_plot_subdir,
    extension = ".png") {
  # Get the path to save a plot
  # base_fname: The base filename of the plot
  # version: The version of the data
  # output_basedir: The base directory for the output
  # plot_subdir: The subdirectory for plots
  # extension: The file extension for the plot
  # Returns: The path to save the plot
  if (version != "") {
    version <- glue(".{version}")
  }
  return(glue("{output_basedir}/{plot_subdir}/{base_fname}{version}{extension}"))
}

load_constraint_metrics <- function(
    version = "v4",
    output_basedir = default_output_basedir,
    data_subdir = default_data_subdir,
    downsamplings = FALSE,
    training_set = FALSE,
    release = TRUE,
    public = TRUE) {
  # Load the constraint metrics data
  # version: The version of the data
  # output_basedir: The base directory for the output
  # data_subdir: The subdirectory for data files
  # downsamplings: Whether to load downsampling data instead of the main data
  # training_set: Whether to load training data instead of the main data
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # Returns: A data frame with constraint metrics data
  if (downsamplings & training_set){
    stop("Only one of downsamplings or training_set can be set to TRUE")
  }
    if (training_set & (release|public)){
    stop("Training set is only available when release and public are both set to FALSE")}

  full_version <- map_data_version[[version]]
  if (!downsamplings) {
    if (version == "v2") {
      base_fname <- glue("gnomad.v{full_version}.lof_metrics.by_gene.txt.bgz")
    } else if (version == "v4") {
      base_fname <- glue("gnomad.v{full_version}.constraint_metrics.tsv")
    }
  } else if (version == "v2") {
    base_fname <- "gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz"
  } else if (version == "v4") {
    base_fname <- glue("gnomad.v{full_version}.downsampling_constraint_metrics.tsv.bgz")
  }
  # Note: v2 training step will be pulled from the reprocessed pipeline
  if (training_set) {
    base_fname <- glue("gnomad.v{full_version}.constraint_training.autosome_par.tsv.bgz")
    pipeline_step = "training_data"
  } else {pipeline_step = "metrics"}

  df <- read.delim(get_or_download_file(
    base_fname,
    version,
    output_basedir,
    data_subdir,
    release = release,
    public = public,
    pipeline_step = pipeline_step
  )
  )
  

  # Convert "true/false" values in "canonical" and/or "mane_select" columns to bools
  for(col_name in c("canonical", "mane_select")) {
    if (col_name %in% colnames(df)) {
      df <- mutate(df, !!col_name := .data[[col_name]] == "true")
    }
  }
  
  # Set absesence of constraint flags to "[]" and append version to constraint flag column name
  if (!downsamplings) {
    if (version == "v2") {
    df <- rename(df, constraint_flags.v2 = constraint_flag) %>%
      mutate(
        constraint_flags.v2 = if_else(
          .data$constraint_flags.v2 == "",
          "[]",
          .data$constraint_flags.v2
        )
      )
    
    # Rename columns in v2 and add version to the column names
    df <- df %>% rename(
           lof.oe_ci.upper.v2 = oe_lof_upper,
           lof.pLI.v2 = pLI,
           lof.z_score.v2 = lof_z,
           lof.oe_ci.upper_bin_decile.v2 = oe_lof_upper_bin,
           syn.exp.v2 = exp_syn,
           mis.exp.v2 = exp_mis,
           lof.exp.v2 = exp_lof,
           syn.obs.v2 = obs_syn,
           mis.obs.v2 = obs_mis,
           lof.obs.v2 = obs_lof,
           lof.oe.v2 = oe_lof,
           mis.oe.v2 = oe_mis,
           syn.oe.v2 = oe_syn)
    }
    
    # Add "v4" version to column names
    if (version == "v4") {
      base_annotations = c("gene", "gene_id", "transcript", "canonical", "mane_select", "level", "transcript_type", "chromosome", "cds_length", "num_coding_exons")
      df <- df %>% rename_with(~ ifelse(. %in% base_annotations, ., paste0(., ".v4")))
    }
  }
    
    # Rename columns in v2 downsamplings
    if (downsamplings && version == "v2") {
      df <- df %>%
            rename(
              syn.exp = exp_syn,
              mis.exp = exp_mis,
              lof.exp = exp_lof,
              syn.obs = obs_syn,
              mis.obs = obs_mis,
              lof.obs = obs_lof
            )
  
  
  
  }
  # Create column for version and set to TRUE
  df <- mutate(df, !!version := TRUE)
  
  return(df)
}

load_all_gene_list_data <- function(
    output_basedir = default_output_basedir,
    gene_list_subdir = default_gene_list_subdir) {
  # Load all gene list data, and clone from GitHub if the data is not available
  # locally
  # output_basedir: The base directory for the output
  # gene_list_subdir: The subdirectory for gene lists
  # Returns: A data frame with gene list data
  gene_list_subdir <- glue("{output_basedir}/{gene_list_subdir}")
  list_dir <- glue("{gene_list_subdir}/lists")
  all_files <- list.files(list_dir, ".+tsv")
  if (length(all_files) == 0) {
    system(
      glue("git clone https://github.com/macarthur-lab/gene_lists.git {gene_list_subdir}")
    )
    all_files <- list.files(list_dir, ".+tsv")
  }
  gene_lists <- map_df(
    all_files,
    function(x) {
      read_tsv(
        paste(list_dir, x, sep = "/"),
        col_names = FALSE,
        col_types = cols()
      ) %>% transmute(gene = toupper(.data$X1), gene_list = str_sub(x, end = -5))
    }
  )

  # Rename gene lists
  gene_lists <- gene_lists %>%
    mutate(
      gene_list = if_else(
        grepl("haploinsufficiency", .data$gene_list),
        "Haploinsufficient",
        .data$gene_list
      )
    ) %>%
    mutate(
      gene_list = fct_recode(
        .data$gene_list,
        "Autosomal Recessive" = "all_ar",
        "Autosomal Dominant" = "all_ad"
      )
    )

  return(gene_lists)
}
