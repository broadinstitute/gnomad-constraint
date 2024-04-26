library(conflicted)
library(dplyr)
library(ggplot2)
library(googleCloudStorageR)
library(grid)
library(rlang)
library(tidyr)
library(forcats, include.only = c("fct_recode", "fct_relevel"))
library(purrr, include.only = c("map_dbl", "map_df"))
library(pROC, include.only = c("roc")) # nolint
library(readr, include.only = c("read_tsv", "cols", "col_double"))
library(scales, include.only = c("comma"))
library(stringr, include.only = c("str_sub"))

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

# Create dataframe of ExAC, gnomAD v2, and gnomAD v4 sample sizes
dataset_sample_sizes <- data.frame(
  intercepts = c(60706, 141456, 807162),
  labels = c("ExAC", "gnomADv2", "gnomADv4"),
  linetype = c("dotted", "dashed", "twodash")
)


label_function <- function(x) {
  paste0(x * 10, "-", x * 10 + 10, "%")
}

summarize_gene_lists <- function(df, metric, version) {
  # Output table of gene list membership and
  # plot gene list membership according to LOEUF decile
  # 'df' is dataframe to use
  # 'version' is version of gnomAD to use for the plot (either 'v2' or 'v4')

  # Filter gene data to specified version
  df <- filter(df, !!sym(version))

  # Remove rows where gene_list or metric is not defined and filter to gene lists of
  # interest
  df <- df %>%
    filter(!is.na(.data$gene_list) & !is.na(!!sym(metric))) %>%
    mutate(metric = label_function(!!sym(metric))) %>%
    filter(
      .data$gene_list %in% c(
        "Haploinsufficient",
        "Autosomal Dominant",
        "Autosomal Recessive",
        "Olfactory Genes"
      )
    )

  ####################################################################
  # Summarize counts of gene lists by oe_lof_upper_bin
  ####################################################################
  # Generate counts of gene list membership
  gene_list_sums <- df %>%
    group_by(.data$metric, .data$gene_list, .drop = FALSE) %>%
    summarise(count = n())
  gene_list_sums <- filter(
    gene_list_sums,
    .data$gene_list %in% c(
      "Haploinsufficient",
      "Autosomal Dominant",
      "Autosomal Recessive",
      "Olfactory Genes"
    )
  )

  return(gene_list_sums)
}

plot_gene_lists <- function(df, gene_lists_to_plot, metric) {
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
      axis.title = element_text(colour = "black", size = 12, face = "bold"),
      axis.text = element_text(colour = "black", size = 10)
    ) +
    scale_fill_manual(values = gene_list_colors, guide = "none") +
    annotate(
      "text",
      4.5,
      top_legend,
      hjust = 0.5,
      vjust = 1,
      label = "Haploinsufficient",
      color = gene_list_colors["Haploinsufficient"]
    ) +
    annotate(
      "text",
      4.5,
      top_legend * 0.88,
      hjust = 0.5,
      vjust = 1,
      label = "Autosomal Recessive",
      color = gene_list_colors["Autosomal Recessive"]
    ) +
    annotate(
      "text",
      4.5,
      top_legend * 0.76,
      hjust = 0.5,
      vjust = 1,
      label = "Olfactory Genes",
      color = gene_list_colors["Olfactory Genes"]
    ) +
    ylab("Percent of gene list (%)") +
    xlab("LOEUF decile (%)") +
    scale_x_discrete(
      labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90")
    )

  return(bar_plot)
}

plot_roc <- function(df, hi_genes, metric, split_seed = 663) {
  # Plot ROC and display value for AUC
  # 'df' is dataframe to use
  # 'metric' is which metric to use for ROC plot (either 'loeuf' or 'pli')

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
      axis.title = element_text(colour = "black", size = 12, face = "bold"),
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
  # Generate plot displaying the percent of genes that would be expected to have a
  # certain number or variants based on sample size
  # df: input dataframe with columns 'n_variants', 'n_required',
  # and 'rank'
  # label: text that will be included at the top of the plot
  # xlimits: Define the limits of the x-axis
  df <- df %>%
    mutate(
      n_variants = forcats::fct_reorder(as.factor(.data$n_variants), .data$n_variants)
    )

  p <- ggplot(df, aes(y = .data$rank, x = .data$n_required, color = .data$n_variants)) +
    geom_line(linewidth = 2) +
    theme_classic() +
    scale_x_log10(labels = comma, limits = xlimits) +
    xlab("Sample size required") +
    ylab("Percent of human genes") +
    geom_vline(
      data = sample_size_df,
      aes(xintercept = .data$intercepts, linetype = .data$linetype),
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
  # Get plot of the percent of genes that have a variable expected number of variants
  # across sample sizes
  # df: dataframe consisting of downsampling data per specific genetic ancestry groups
  # (has to include 'global')
  # version: version of gnomAD to use for the plot
  # Returns: plot of the percent of genes that would be expected to have a certain
  # number or variants based on sample size

  # Filter to rows where a respective gene has at least 1 lof, mis, and syn variant
  # within all its respective rows
  # Convert expected counts and n downsamplings to log scale
  # Fit linear models where expected counts are a function of downsampling size for
  # each gene
  df <- df %>%
    group_by(.data$gene) %>%
    filter(min(.data$exp_lof) > 0 & min(.data$exp_mis) > 0 & min(.data$exp_syn) > 0) %>%
    mutate(
      log_exp_lof = log10(.data$exp_lof),
      log_exp_mis = log10(.data$exp_mis),
      log_exp_syn = log10(.data$exp_syn),
      log_n = log10(.data$downsampling)
    ) %>%
    summarize(
      lof_fit = list(lm(.data$log_exp_lof ~ .data$log_n)),
      mis_fit = list(lm(.data$log_exp_mis ~ .data$log_n)),
      syn_fit = list(lm(.data$log_exp_syn ~ .data$log_n))
    )

  # Extract slope and intercept for lof variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
  # why need this step to sum?
  gene_lof_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(.data$lof_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(.data$lof_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(.data$gene) %>%
    summarize(slope = sum(.data$slope), intercept = sum(.data$intercept))

  # Extract slope and intercept for missense variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
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
    # fitted model coefficients 'data' is the dataframe to use (should contain columns
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
        cols = -.data$gene,
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