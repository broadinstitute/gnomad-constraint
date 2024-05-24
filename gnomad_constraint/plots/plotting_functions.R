library(conflicted)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr, include.only = c("map_dbl"))
library(pROC, include.only = c("roc")) # nolint
library(scales, include.only = c("comma"))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("%>%", "dplyr")

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

default_gene_lists <- c(
  "Haploinsufficient",
  "Autosomal Dominant",
  "Autosomal Recessive",
  "Olfactory Genes"
)
default_gene_lists_to_plot <- c(
  "Haploinsufficient",
  "Autosomal Recessive",
  "Olfactory Genes"
)

# Create dataframe of ExAC, gnomAD v2, and gnomAD v4 sample sizes
dataset_sample_sizes <- data.frame(
  intercepts = c(60706, 141456, 807162),
  labels = c("ExAC", "gnomADv2", "gnomADv4"),
  linetype = c("dotted", "dashed", "twodash")
)

label_function <- function(x) {
  # Function to generate labels for LOEUF or pLI deciles
  # x: Decile to label
  # Returns: label for decile
  paste0(x * 10, "-", x * 10 + 10, "%")
}

summarize_gene_lists <- function(
    df,
    metric,
    version,
    gene_lists_to_summarize = default_gene_lists) {
  # Get table of gene list membership counts
  # df: Dataframe containing gene list membership data
  # metric: Metric to use for the plot (either 'loeuf' or 'pli')
  # version: Version of gnomAD to use for the plot (either 'v2' or 'v4')
  # gene_lists_to_summarize: List of gene lists to summarize
  # Returns: Dataframe with counts of gene list membership

  # Filter gene data to specified version
  df <- filter(df, !!sym(version))

  # Remove rows where gene_list or metric is not defined and filter to gene lists of
  # interest
  df <- df %>%
    filter(!is.na(.data$gene_list) & !is.na(!!sym(metric))) %>%
    mutate(metric = label_function(!!sym(metric))) %>%
    filter(.data$gene_list %in% gene_lists_to_summarize)

  ####################################################################
  # Summarize counts of gene lists by oe_lof_upper_bin
  ####################################################################
  # Generate counts of gene list membership
  gene_list_sums <- df %>%
    group_by(.data$metric, .data$gene_list, .drop = FALSE) %>%
    summarise(count = n()) %>%
    filter(.data$gene_list %in% gene_lists_to_summarize)

  return(gene_list_sums)
}

plot_gene_lists <- function(
    df,
    metric,
    gene_lists_to_plot = default_gene_lists_to_plot) {
  # Plot gene list membership by decile of LOEUF or pLI
  # df: Dataframe containing gene list membership counts from
  # summarize_gene_lists
  # metric: Metric to use for the plot (either 'loeuf' or 'pli')
  # gene_lists_to_plot: List of gene lists to plot
  # Returns: ggplot object

  # Convert counts to proportions
  props <- df %>%
    group_by(.data$gene_list) %>%
    mutate(prop_in_bin = .data$count / sum(.data$count))
  props <- filter(props, .data$gene_list %in% gene_lists_to_plot)

  ####################################################################
  # Plot gene list per decile
  ####################################################################
  top_legend <- max(props$prop_in_bin * 100)
  bar_plot <- ggplot(
    props,
    aes(x = .data$metric, y = .data$prop_in_bin * 100, fill = .data$gene_list)
  ) +
    geom_bar(position = "dodge", stat = "identity", width = 0.9) +
    theme_classic() +
    theme(
      axis.title = element_text(colour = "black", size = 12),
      axis.text = element_text(colour = "black", size = 10)
    ) +
    scale_fill_manual(values = gene_list_colors, guide = "none")

  top_legend_adj <- 1.0
  for (gene_list in gene_lists_to_plot) {
    bar_plot <- bar_plot +
      annotate(
        "text",
        4.5,
        top_legend * top_legend_adj,
        hjust = 0.5,
        vjust = 1,
        label = gene_list,
        color = gene_list_colors[gene_list]
      )
    top_legend_adj <- top_legend_adj - 0.12
  }

  bar_plot <- bar_plot +
    ylab("Percent of gene list (%)") +
    xlab("LOEUF decile") +
    scale_x_discrete(
      labels = c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th")
    )

  return(bar_plot)
}

plot_roc <- function(df, hi_genes, metric, split_seed = 663) {
  # Plot ROC with a value of AUC displayed on the plot for a given metric,
  # using haploinsufficient genes as the positive class
  # df: Dataframe containing constraint data for genes
  # hi_genes: Dataframe containing haploinsufficient genes
  # metric: Metric to use for ROC plot (either 'loeuf' or 'pli')
  # split_seed: Seed for splitting data into training and testing sets
  # Returns: ggplot object

  # Define haploinsufficient genes
  df <- mutate(df, hi_gene = .data$gene %in% hi_genes$gene)

  # Split data into a training and testing set
  set.seed(split_seed)
  sample <- sample(c(TRUE, FALSE), nrow(df), replace = TRUE, prob = c(0.7, 0.3))
  train <- df[sample, ]
  test <- df[!sample, ]

  # Train generalized linear model
  formula <- as.formula(paste("hi_gene ~", metric))
  m <- glm(formula, data = train, family = "binomial")

  # Apply model predictions
  predicted <- predict(m, test, type = "response")

  # Plot ROC curve
  p <- roc(test$hi_gene, predicted, print.auc = TRUE)

  return(p)
}

combine_roc_plots <- function(
    roc1,
    roc2,
    version1,
    version2,
    title_label,
    color1 = "darkorange1",
    color2 = "darkorchid3") {
  # Combine two ROC plots into a single plot
  # roc1: ROC object from pROC package
  # roc2: ROC object from pROC package
  # version1: Version label for ROC1
  # version2: Version label for ROC2
  # title_label: Title for the plot
  # color1: Color for ROC1
  # color2: Color for ROC2
  # Returns: ggplot object

  # Combine ROC outputs
  roc_data1 <- as.data.frame(cbind(roc1$sensitivities, roc1$specificities))
  roc_data2 <- as.data.frame(cbind(roc2$sensitivities, roc2$specificities))
  roc_data1$version <- version1
  roc_data2$version <- version2
  all_rocs <- rbind(roc_data1, roc_data2)
  colnames(all_rocs) <- c("sensitivity", "specificity", "version")

  # Define AUC labels
  auc1 <- paste("AUC: ", round(roc1$auc[1], 2), sep = "")
  auc2 <- paste("AUC: ", round(roc2$auc[1], 2), sep = "")

  # Plot ROC output
  p <- ggplot(
    all_rocs,
    aes(1 - .data$specificity, .data$sensitivity, color = version)
  ) +
    geom_line(linewidth = 1, alpha = 0.9) +
    theme_classic() +
    theme(
      axis.title = element_text(colour = "black", size = 12),
      axis.text = element_text(colour = "black", size = 10)
    ) +
    scale_color_manual(values = c(color1, color2)) +
    labs(
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Version",
      title = title_label
    ) +
    annotate("text", x = .70, y = .25, label = auc1, color = color1) +
    annotate("text", x = .70, y = .20, label = auc2, color = color2)

  return(p)
}

expected_projections <- function(
    df,
    label = "pLoF",
    sample_size_df = dataset_sample_sizes,
    xlimits = c(100, 1e8)) {
  # Plot displaying the percent of genes that would be expected to have a
  # certain number or variants based on sample size
  # df: Dataframe with columns 'n_variants', 'n_required', and 'rank'
  # label: Text that will be included at the top of the plot
  # sample_size_df: Dataframe with columns 'intercepts', 'linetype', and
  # 'label' for vertical lines on the plot that show the sample size of some
  # key datasets
  # xlimits: Define the limits of the x-axis
  # Returns: ggplot object
  df <- df %>%
    mutate(
      n_variants = forcats::fct_reorder(as.factor(.data$n_variants), .data$n_variants)
    )

  p <- ggplot(df, aes(y = .data$rank, x = .data$n_required, color = .data$n_variants)) +
    geom_line(linewidth = 2) +
    theme_classic() +
    scale_x_log10(labels = comma, limits = xlimits) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("Sample size required") +
    ylab("Percent of human genes") +
    geom_vline(
      data = sample_size_df,
      aes(xintercept = .data$intercepts, linetype = .data$labels),
      key_glyph = "path"
    )

  # Add manual linetype scale for legend
  p <- p +
    scale_linetype_manual(
      name = "Database size",
      values = sample_size_df$linetype,
      labels = sample_size_df$labels
    ) +
    scale_color_discrete(name = ">= N variants\nexpected", h = c(40, 120)) +
    annotate("text", x = xlimits[1], y = 1, hjust = 0, vjust = 1, label = label)

  return(p)
}

plot_projected_sample_size <- function(df) {
  # Plot of the percent of genes that have a variable expected number of
  # variants across sample sizes
  # df: Dataframe consisting of downsampling data per specific genetic ancestry groups
  # (has to include 'global')
  # Returns: ggplot object of the percent of genes that would be expected to have a
  # certain number or variants based on sample size

  # Filter to rows where a respective gene has at least 1 lof, mis, and syn variant
  # within all its respective rows
  df <- df %>%
    group_by(.data$gene) %>%
    filter(min(.data$exp_lof) > 0 & min(.data$exp_mis) > 0 & min(.data$exp_syn) > 0)

  # Convert expected counts and n downsamplings to log scale
  df <- df %>%
    mutate(
      log_exp_lof = log10(.data$exp_lof),
      log_exp_mis = log10(.data$exp_mis),
      log_exp_syn = log10(.data$exp_syn),
      log_n = log10(.data$downsampling)
    )

  # Fit linear models where expected counts are a function of downsampling size for
  # each gene
  df <- df %>%
    summarize(
      lof_fit = list(lm(.data$log_exp_lof ~ .data$log_n)),
      mis_fit = list(lm(.data$log_exp_mis ~ .data$log_n)),
      syn_fit = list(lm(.data$log_exp_syn ~ .data$log_n))
    )

  # Extract slope and intercept for lof variants
  gene_lof_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(.data$lof_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(.data$lof_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(.data$gene) %>%
    summarize(slope = sum(.data$slope), intercept = sum(.data$intercept))

  # Extract slope and intercept for missense variants
  gene_mis_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(.data$mis_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(.data$mis_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(.data$gene) %>%
    summarize(slope = sum(.data$slope), intercept = sum(.data$intercept))

  ####################################################################
  # Process predictions
  ####################################################################
  post_process_predictions <- function(data) {
    # Transform predictions at specific points ('5', '10', '20', '50', '100') using
    # fitted model coefficients
    # data: is the dataframe to use (should contain columns
    # 'gene', 'slope' and 'intercept')
    # Calculates the percentile rank of each value in the n_required column within the
    # respective n_variants group
    data %>%
      mutate(
        `5` = 10^((log10(5) - .data$intercept) / .data$slope),
        `10` = 10^((log10(10) - .data$intercept) / .data$slope),
        `20` = 10^((log10(20) - .data$intercept) / .data$slope),
        `50` = 10^((log10(50) - .data$intercept) / .data$slope),
        `100` = 10^((log10(100) - .data$intercept) / .data$slope)
      ) %>%
      select(-all_of(c("slope", "intercept"))) %>%
      # Select all columns for pivot except 'gene'
      pivot_longer(
        cols = -all_of(c("gene")),
        names_to = "n_variants",
        values_to = "n_required"
      ) %>%
      # Clean and convert to numeric if needed
      mutate(
        n_variants = as.numeric(gsub("X", "", .data$n_variants))
      ) %>%
      group_by(.data$n_variants) %>%
      mutate(
        rank = percent_rank(.data$n_required)
      ) %>%
      ungroup()
  }

  samples_required_lof <- post_process_predictions(gene_lof_fit)
  samples_required_mis <- post_process_predictions(gene_mis_fit)

  ####################################################################
  # Plot projections
  ####################################################################
  lof_projections <- expected_projections(samples_required_lof, "pLoF")
  missense_projections <- expected_projections(samples_required_mis, "Missense")

  return(list(lof = lof_projections, mis = missense_projections))
}


plot_decile_change <- function(df) {
  # Plot of the change in LOEUF deciles between gnomAD versions
  # df: Dataframe consisting of LOEUF deciles, with the v2 values defined by
  # 'oe_lof_upper_bin' and the v4 values defined by 'lof.oe_ci.upper_bin_decile'
  # Returns: ggplot object of the LOEUF decile change from v2 to v4

  # Calculate the decile change by subtracting the v4 value from the v2 value
  df <- df %>% mutate(decile_change = .data$oe_lof_upper_bin - .data$lof.oe_ci.upper_bin_decile)
  # Plot decile changes
  p <- ggplot(
    df,
    aes(x = .data$decile_change)
  ) +
    geom_histogram(stat = "count") +
    theme_classic() +
    theme(
      axis.title = element_text(colour = "black", size = 25),
      axis.text = element_text(colour = "black", size = 20)
    ) +
    scale_x_continuous(breaks = seq(-10, 8, by = 2)) +
    ylab("Count") +
    xlab("Decile change from v2 to v4")

  return(p)
}



plot_metric_comparison <- function(df) {
  # Plot comparison of metrics between gnomAD versions
  # df: Dataframe consisting of 'metric_name' and the corresponding
  # values in 'v2' and 'v4'
  # Returns: ggplot object of metric comparison between v2 to v4

  p <- ggplot(df, aes(x = .data$v2, y = .data$v4)) +
    geom_point(size=0.75, alpha=.25) +
    facet_wrap(.~metric_name, scale="free") +
    theme_classic() +
    theme(axis.text = element_text(size = 15, colour = "black"),
          axis.title = element_text(size = 20, colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15, colour = "black"),
          plot.margin = unit(c(.65, .65, .65, .65), "lines")
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(x = "Value in gnomAD v2", y = "Value in gnomAD v4")

  return(p)
}



plot_observed_vs_expected <- function(df, version) {
    # Plot the observed vs expected values
    # df: Dataframe consisting of metrics to plot in 'metric_name', the
    # observed counts in 'obs', and the expected counts in 'exp', and
    # max values of the two in 'max_limit'
    # Returns: ggplot object of observed vs expected values for the specified version

  # Create list to store plots in
  plot_list <- list()
  # Create dataframe to store correlation results in
  correlation_results <- data.frame(metric = character(), correlation = numeric(), stringsAsFactors = FALSE)
  # Plot observed vs expected coutns for each metric
  for(metric in unique(df$metric_name)) {
    # Filter dataset to the specified metric
    data_subset <- filter(df, .data$metric_name == metric)
    # Pull out the max value of the metric
    max_limit <- max(data_subset$max_limit, na.rm = TRUE)

    # Calculte correlation between observed and expected counts
    correlation <- cor(data_subset$exp, data_subset$obs, method = "pearson")
    correlation_results <- rbind(correlation_results, data.frame(metric = metric, correlation = correlation))

    # Plot results
    p <- ggplot(data_subset, aes(x = exp, y = obs)) +
      geom_point(size = 0.75) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      ggtitle(metric) +  # Add title to identify the metric
      coord_fixed(ratio = 1) +
      xlim(0, max_limit) +
      ylim(0, max_limit) +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", hjust=1, angle=45),
        axis.title = element_blank(),
        plot.margin = unit(c(1,1,1,2.2), "lines"),
        plot.title = element_text(hjust = 0.5)
      )
    # Add plot to the plot list
    plot_list[[metric]] <- p
  }

  # Combine all plots in the plot list
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = 'v')
  combined_plot<- combined_plot + draw_label("Expected Variants", x=0.5, vjust= 7, angle= 0, size=16)
  combined_plot<- combined_plot + draw_label("Observed\n Variants" , y=0.5, vjust= -7.75, angle=90, size=16)
  # Print the correlation results
  print(glue("Correlation results for  observed vs expected counts in {version}:"))
  print(correlation_results)
  return(combined_plot)
}
