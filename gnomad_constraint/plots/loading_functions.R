library(conflicted)
library(dplyr)
library(gargle)
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
  v4 = "4.1",
  v4_tmp = "4.1"
)

default_output_path <- "."
default_data_dir <- "data"
default_plot_dir <- "plots"
default_gene_lists <- "data/gene_lists"
bucket <- "gs://gnomad"

setup_directories <- function(
    output_path = default_output_path,
    data_dir = default_data_dir,
    plot_dir = default_plot_dir,
    gene_lists = default_gene_lists) {
  # Create the output directory and subdirectories if they don't exist
  # output_path: The base directory for the output
  # data_dir: The subdirectory for data files
  # plot_dir: The subdirectory for plot files
  # gene_lists: The subdirectory for gene list files
  # Returns: None
  suppressWarnings(dir.create(output_path))
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, data_dir)))
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, plot_dir)))
  suppressWarnings(dir.create(sprintf("%s/%s/", output_path, gene_lists)))
}

get_data_url <- function(version = "v4", release = TRUE, public = TRUE) {
  # Get the URL for the data files
  # version: The version of the data
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # Returns: The URL for the data files

  # TODO: Modify if there are changes for v4.1
  full_version <- map_data_version[[version]]
  if (release && public) {
    data_path <- sprintf(
      "gs://gcp-public-data--gnomad/release/%s/constraint/",
      full_version
    )
  } else if (version == "v4_tmp") {
    data_path <- "gs://gnomad-kristen/constraint/"
  } else if (release && !public) {
    data_path <- sprintf("gs://gnomad/release/%s/constraint/", full_version)
  } else if (!release && !public) {
    data_path <- sprintf("gs://gnomad/v%s/constraint/metrics/tsv/", full_version)
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
    output_path = default_output_path,
    data_dir = default_data_dir,
    subfolder = "",
    local_name = "",
    release = TRUE,
    public = TRUE,
    gcs_authentication_token = NULL) {
  # Get a local file, and if it doesn't exist, try downloading it from Google
  # Cloud Storage (GCS)
  # base_fname: The base filename of the file
  # version: The version of the data
  # output_path: The base directory for the output
  # data_dir: The subdirectory for data files
  # subfolder: The subfolder in the data directory
  # local_name: The local name of the file
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # gcs_authentication_token: The path to the saved token for GCS authentication
  # Returns: The path to the local file for requested data
  local_path <- paste0(
    output_path,
    "/",
    data_dir,
    "/",
    ifelse(local_name != "", local_name, base_fname)
  )
  # To set up authentication for Google Cloud Storage (GCS) access, go to
  # https://console.cloud.google.com/apis/credentials and download the JSON file for
  # the tidyverse-gcs-access under OAuth 2.0 Client IDs. Then pass that file to the
  # interactive_authenticate_gcs function in an interactive session. The saved token
  # can be passed to this function.

  # TODO: Modify when v4.1 files are ready
  if (version == "v2") {
    data_bucket <- "gs://gcp-public-data--gnomad"
  } else if (version == "v4") {
    data_bucket <- "gs://gnomad"
  } else if (version == "v4_tmp") {
    data_bucket <- "gs://gnomad-kristen"
  }
  if (!file.exists(local_path)) {
    if (!missing(gcs_authentication_token)) {
      authenticate_gcs(gcs_authentication_token)
    }
    remote_path <- paste0(
      get_data_url(version, release = release, public = public),
      subfolder,
      base_fname
    )
    print(paste0("Downloading ", remote_path, " to ", local_path))
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
  # Get the path to save a plot
  # base_fname: The base filename of the plot
  # version: The version of the data
  # output_path: The base directory for the output
  # plot_dir: The subdirectory for plots
  # extension: The file extension for the plot
  # Returns: The path to save the plot
  if (version == "") {
    format_path <- "%s/%s/%s%s%s"
  } else {
    format_path <- "%s/%s/%s.%s%s"
  }
  return(sprintf(format_path, output_path, plot_dir, base_fname, version, extension))
}

load_constraint_metrics <- function(
    version = "v4",
    output_path = default_output_path,
    data_dir = default_data_dir,
    downsamplings = FALSE,
    release = TRUE,
    public = TRUE) {
  # Load the constraint metrics data
  # version: The version of the data
  # output_path: The base directory for the output
  # data_dir: The subdirectory for data files
  # downsamplings: Whether to load downsampling data instead of the main data
  # release: Whether to use the release or development version
  # public: Whether to use the public or non-public version
  # Returns: A data frame with constraint metrics data

  # TODO: Change if files change or are added for v4.1
  full_version <- map_data_version[[version]]
  if (downsamplings && version == "v2") {
    base_fname <- "gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz"
  } else if (downsamplings) {
    base_fname <- sprintf(
      "gnomad.v%s.downsampling_constraint_metrics.txt.bgz",
      full_version
    )
  } else if (version == "v4_tmp") {
    base_fname <- "extra_annotations.tsv"
  } else {
    base_fname <- sprintf("gnomad.v%s.lof_metrics.by_gene.txt.bgz", full_version)
  }

  fname <- get_or_download_file(
    base_fname,
    version,
    output_path,
    data_dir,
    release = release,
    public = public
  )

  return(read.delim(fname))
}

load_all_gene_list_data <- function(
    output_path = default_output_path,
    gene_lists = default_gene_lists) {
  # Load all gene list data, and clone from GitHub if the data is not available
  # locally
  # output_path: The base directory for the output
  # gene_lists: The subdirectory for gene lists
  # Returns: A data frame with gene list data
  gene_list_dir <- sprintf("%s/%s", output_path, gene_lists)
  list_dir <- paste0(gene_list_dir, "/lists")
  all_files <- list.files(list_dir, ".+tsv")
  if (length(all_files) == 0) {
    system(
      paste0(
        "git clone https://github.com/macarthur-lab/gene_lists.git ",
        gene_list_dir
      )
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
