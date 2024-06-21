library(dplyr)


####################################################################
# Define function for applying transcript filters
####################################################################
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


####################################################################
# Reformat data to use for plotting comparison of metrics between v2 and v4
####################################################################
reformat_for_metric_comparison <- function(comparison_df) {
  # Reformat data to use for plotting comparison of metrics between v2 and v4
  # comparison_df: Dataframe of LOEUF and Z-score values for v2 and v4
  # Returns: Reformatted data for plotting comparisons of metrics between v2 and v4
  
  # Pivot to longer df with metric name and value
  comparison_df <- comparison_df %>%
    pivot_longer(
      cols = c(lof.oe_ci.upper.v2, lof.oe_ci.upper.v4, lof.z_score.v2, lof.z_score.v4),
      names_to = "metric_name",
      values_to = "value"
    )
  
  # Pull out version
  comparison_df <- comparison_df %>% mutate(
    version = if_else(grepl("v4", metric_name), "v4", "v2")
  )
  
  # Rename metrics in df to be the same for v2 and v4
  comparison_df$metric_name <- gsub("\\.v(2|4)$", "", comparison_df$metric_name)
  
  comparison_df <- comparison_df %>%
    mutate(metric_name = case_when(
      .data$metric_name == "lof.oe_ci.upper" ~ "LOEUF",
      .data$metric_name == "lof.z_score" ~ "pLoF Z-score",
      TRUE ~ .data$metric_name # Keeps all other values as they are
    ))
  
  comparison_df <- comparison_df %>% select(-v2, -v4)
  
  # Pivot to wider format
  comparison_df <- comparison_df %>% pivot_wider(
    names_from = version, values_from = value
  )
  
  # Calculate correlations for each metric between v2 and v4
  comparison_df %>%
    group_by(metric_name) %>%
    summarize(
      correlation = cor(.data$v2, .data$v4, method = "pearson", use = "complete.obs")
    )
  return(df)
}
####################################################################
# Reformat data to use for plotting observed vs expected values
####################################################################
reformat_for_observed_vs_expected <- function(df, version) {
  # Reformat data to use for plotting observed vs expected values
  # df: Dataframe of observed and expected values for syn, mis, and lof variants
  # version: Version of gnomAD to which the data has been filtered to (either 'v2' or 'v4')
  # Returns: Reformatted data for plotting observed vs expected values
  
  # Reformat column names 
  df <- df %>% select(
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
  df <- df %>%
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
  df <- df %>%
    mutate(metric_name = recode(.data$metric_name, "lof" = "pLoF")) %>%
    mutate(metric_name = factor(.data$metric_name, levels = c("syn", "mis", "pLoF")))

  # Find the max value of the observed and expected values for each metric
  max_values <- df %>%
    group_by(metric_name) %>%
    summarise(
      max_exp = max(.data$exp, na.rm = TRUE),
      max_obs = max(.data$obs, na.rm = TRUE)
    ) %>%
    mutate(max_limit = pmax(max_exp, max_obs))

  # Join these max values back to the original data to use in plotting
  df <- df %>%
    left_join(max_values, by = "metric_name")
  
  return(df)
}