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
      cols = c(lof.oe_ci.upper.v2, lof.oe_ci.upper.v4, lof.z_score.v2, lof.z_score.v4,lof.oe.v2, lof.oe.v4),
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
      .data$metric_name == "lof.oe" ~ "pLoF o/e",
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
  return(comparison_df)
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


####################################################################
# Get predicted proportion observed
####################################################################
get_predicted_proportion_observed <- function (
    df,
    coverage_metric ="exome_coverage",
    high_coverage_cutoff = 30) {
  # Calculate the predicted proportion of observed variants using well-covered sites
  # df: Dataframe consisting of observed and possible variants per each context (the output of the create_training_set step in the constriant pipeline)
  # coverage_metric: Metric to use to determine well-covered sites (should correspond to column name in the dataframe). Examples: "exome_coverage" or "exomes_AN_percent". Default is 'exome_coverage'. 
  # high_coverage_cutoff: Cutoff for determining well-covered sites in the specified coverage_metric. Default is 30.
  # Returns: Dataframe with the predicted proportion of observed variants in the 'pred_prop_observed' column.
  
  # Calculate proportion observed variants (observed/possible) at well covered sites
  po_data <- df %>%
    filter(!!sym(coverage_metric) >= high_coverage_cutoff) %>%
    group_by(context, ref, alt, methylation_level, mu_snp, mutation_type, cpg) %>%
    summarize(
      obs = sum(observed_variants, na.rm = TRUE),
      poss = sum(possible_variants, na.rm = TRUE),
      prop_observed = obs / poss,
      .groups = 'drop'
    )
  
  # Group by cpg and perform linear regression on the well-covered sites (lm(prop_observed ~ mu_snp)
  high_coverage_models <- po_data %>%
    group_by(cpg) %>%
    do({
      model <- lm(prop_observed ~ mu_snp, data = .)
      data.frame(
        cpg = unique(.$cpg),
        term = c("intercept", "mu_snp"),
        estimate = coef(model)
      )
    }) %>%
    ungroup()
  
  # Join models from well-covered sites to the original dataframe and calculate the predicted proportion observed
  data_with_predictions <- df %>%
    left_join(high_coverage_models, by = "cpg") %>%
    group_by_at(vars(-term, -estimate)) %>%
    summarize(
      pred_prop_observed = sum((term == 'mu_snp') * mu_snp * estimate + (term == 'intercept') * estimate),
      .groups = 'drop'
    ) %>%
    ungroup()
  
  return(data_with_predictions)
  
  
}