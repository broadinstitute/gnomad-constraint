label_function <- function(x) {
  paste0(x * 10, "-", x * 10 + 10, "%")
}

plot_gene_lists <- function(df, version, plot_path) {
  # Output table of gene list membership and
  # plot gene list membership according to LOEUF decile
  # 'df' is dataframe to use
  # 'version' is version of gnomAD to use for the plot (either 'v2' or 'v4')

  if (version == "v2") {
    metric <- "oe_lof_upper_bin"
  } else {
    metric <- "lof.oe_ci.upper_bin_decile"
  }

  gene_data$lof.oe_ci.upper_bin_decile
  sort(colnames(gene_data))

  # Filter dataframe to specified version
  df <- filter(df,!!sym(version))

  # Remove rows where gene_list or metric is not defined and filter to gene lists of interest
  df <- df %>%
    filter(!is.na(gene_list) & !is.na(!!sym(metric))) %>%
    mutate(metric = label_function(!!sym(metric))) %>%
    filter(
      gene_list %in% c(
        "Haploinsufficient",
        "Autosomal Dominant",
        "Autosomal Recessive",
        "Olfactory Genes"
      )
    )
  summary_df <- df %>% select(gene, gene_list,!!sym(metric))

  ####################################################################
  # Summarize counts of gene lists by oe_lof_upper_bin
  ####################################################################
  # Generate counts of gene list membership
  gene_list_sums <- df %>%
    group_by(metric, gene_list, .drop = FALSE) %>%
    summarise(count = n())
  gene_list_sums <- dplyr::filter(
      gene_list_sums,
      gene_list %in% c(
        "Haploinsufficient",
        "Autosomal Dominant",
        "Autosomal Recessive",
        "Olfactory Genes"
      )
    )

  summary_gene_list_per_sums <- gene_list_sums %>% spread(gene_list, count)

  # Write out table of gene list membership
  write.table(
    summary_gene_list_per_sums,
    file = paste("gene_list_counts_", version, ".txt", sep = ""),
    quote = FALSE
  )

  # Convert counts to proportions
  props <- gene_list_sums %>%
    group_by(gene_list) %>%
    mutate(prop_in_bin = count / sum(count))
  props <- filter(props, gene_list %in% gene_lists_to_plot)

  ####################################################################
  # Plot gene list per decile
  ####################################################################
  top_legend <- max(props$prop_in_bin * 100)

  bar_plot <- ggplot(props, aes(x = metric, y = prop_in_bin * 100, fill = gene_list)) +
    geom_bar(position = "dodge",
             stat = "identity",
             width = 0.9) +
    theme_classic() +
    theme(
      axis.title = element_text(
        colour = "black",
        size = 12,
        face = "bold"
      ),
      axis.text = element_text(colour = "black", size = 10)
    ) +
    scale_fill_manual(values = gene_list_colors, guide = F) +
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
    scale_x_discrete(labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90"))

  ggsave(
    bar_plot,
    filename = plot_path,
    dpi = 300,
    width = 11,
    height = 6,
    units = "in"
  )
}

plot_roc <- function(df, metric, plot_path) {
  # Plot ROC and display value for AUC
  # 'df' is dataframe to use
  # 'metric' is which metric to use for ROC plot (either 'loeuf' or 'pli')

  # Define haploinsufficient genes
  df <- mutate(df, hi_gene = gene %in% hi_genes$gene)

  # Define column names based on specified metric
  if (metric == "loeuf") {
    v2_metric <- "oe_lof_upper"
    v4_metric <- "lof.oe_ci.upper"
    title_label <- "LOEUF"
  } else if (metric == "pli") {
    v2_metric <- "pLI"
    v4_metric <- "lof.pLI"
    title_label <- "pLI"
  }

  # Filter to where the metric is defined in both v2 and v4
  df <- df %>% dplyr::filter(!is.na(!!sym(v2_metric)) &
                           !is.na(!!sym(v4_metric)))

  # Split data into a training and testing set
  set.seed(663)
  sample <- sample(c(TRUE, FALSE),
           nrow(df),
           replace = TRUE,
           prob = c(0.7, 0.3))
  train <- df[sample,]
  test <- df[!sample,]

  # Train generalized linear model
  v2_formula <- as.formula(paste("hi_gene ~", v2_metric))
  v4_formula <- as.formula(paste("hi_gene ~", v4_metric))

  m2 <- glm(v2_formula, data = train, family = "binomial")
  m4 <- glm(v4_formula, data = train, family = "binomial")

  # Apply model predictions
  predicted_v2 <- predict(m2, test, type = "response")
  predicted_v4 <- predict(m4, test, type = "response")

  # Plot ROC curve and display AUC
  roc_v2 <- roc(test$hi_gene,
        predicted_v2,
        print.auc = T,
        plot = TRUE)
  roc_v4 <- roc(test$hi_gene,
        predicted_v4,
        print.auc = T,
        plot = TRUE)

  # Combine ROC outputs
  roc_data_v2 <- as.data.frame(cbind(roc_v2$sensitivities, roc_v2$specificities))
  roc_data_v4 <- as.data.frame(cbind(roc_v4$sensitivities, roc_v4$specificities))
  roc_data_v2$version <- "v2"
  roc_data_v4$version <- "v4"
  all_rocs <- rbind(roc_data_v2, roc_data_v4)
  colnames(all_rocs) <- c("sensitivity", "specificity", "version")

  # Define AUC labels
  v2_auc <- paste("AUC: ", round(roc_v2$auc[1], 2), sep = "")
  v4_auc <- paste("AUC: ", round(roc_v4$auc[1], 2), sep = "")

  # Define plot colors for gnomAD versions
  v2_color <- "darkorange1"
  v4_color <- "darkorchid3"

  # Plot ROC output
  roc_plot <- ggplot(all_rocs, aes(1 - specificity, sensitivity, color = version)) +
    geom_line(size = 1, alpha = 0.9) +
    theme_classic() +
    theme(
      axis.title = element_text(
        colour = "black",
        size = 12,
        face = "bold"
      ),
      axis.text = element_text(colour = "black", size = 10)
    ) +
    scale_color_manual(values = c(v2_color, v4_color)) +
    labs(
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Version",
      title = title_label
    ) +
    annotate(
      "text",
      x = .70,
      y = .25,
      label = v2_auc,
      color = v2_color
    ) +
    annotate(
      "text",
      x = .70,
      y = .20,
      label = v4_auc,
      color = v4_color
    )

  ggsave(
    roc_plot,
    filename = plot_path,
    dpi = 300,
    width = 6,
    height = 6,
    units = "in"
  )
}

plot_projected_sample_size <- function(df, version) {
  # Get plot of the percent of genes that have a variable expected number of variants
  # across sample sizes
  # df: dataframe consisting of downsampling data per specific genetic ancestry groups
  # (has to include 'global')
  # version: version of gnomAD to use for the plot
  # Returns: plot of the percent of genes that would be expected to have a certain
  # number or variants based on sample size

  # Filter to canoncial/MANE Select transcripts and fit linear models
  if (version == "v2") {
    df <- filter(
        df,
        (canonical == "true") & (pop == "global") & (downsampling >= 100)
      )
  } else {
    df <- filter(
        df,
        (mane_select == "true") &
          (grepl("^ENST", transcript)) &
          (gen_anc == "global") & (downsampling > 100) & (!is.na(gene))
      )
  } # Remove 8 genes with missing gene names

  # Filter to rows where a respective gene has at least 1 lof, mis, and syn variant within all its respective rows
  # Convert expected counts and n downsamplings to log scale
  # Fit linear models where expected counts are a function of downsampling size for each gene
  df <- df %>%
    group_by(gene) %>%
    filter(min(exp_lof) > 0 &
      min(exp_mis) > 0 & min(exp_syn) > 0) %>%
    mutate(
      log_exp_lof = log10(exp_lof),
      log_exp_mis = log10(exp_mis),
      log_exp_syn = log10(exp_syn),
      log_n = log10(downsampling)
    ) %>%
    summarize(
      lof_fit = list(lm(log_exp_lof ~ log_n)),
      mis_fit = list(lm(log_exp_mis ~ log_n)),
      syn_fit = list(lm(log_exp_syn ~ log_n))
    )

  # Extract slope and intercept for lof variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
  # why need this step to sum?
  gene_lof_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(lof_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(lof_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(gene) %>%
    summarize(slope = sum(slope), intercept = sum(intercept))

  # Extract slope and intercept for missense variants
  # Extract slope (coefficient for log_n)
  # Extract intercept
  gene_mis_fit <- df %>%
    mutate(
      slope = purrr::map_dbl(mis_fit, ~ .x$coefficients[2]),
      intercept = purrr::map_dbl(mis_fit, ~ .x$coefficients[1])
    ) %>%
    group_by(gene) %>%
    summarize(slope = sum(slope), intercept = sum(intercept))

  ####################################################################
  # Process predictions
  ####################################################################
  post_process_predictions <- function(data) {
    # Transform predictions at specific points ('5', '10', '20', '50', '100') using fitted model coefficients
    # 'data' is the dataframe to use (should contain columns 'gene', 'slope' and 'intercept')
    # Calculates the percentile rank of each value in the n_required column within the respective n_variants group
    data %>%
      mutate(
        `5` = 10^((log10(5) - intercept) / slope),
        `10` = 10^((log10(10) - intercept) / slope),
        `20` = 10^((log10(20) - intercept) / slope),
        `50` = 10^((log10(50) - intercept) / slope),
        `100` = 10^((log10(100) - intercept) / slope)
      ) %>%
      select(-slope, -intercept) %>%
      gather("n_variants", "n_required", -gene) %>%
      mutate(n_variants = as.numeric(n_variants)) %>%
      group_by(n_variants) %>%
      mutate(rank = percent_rank(n_required)) %>%
      ungroup()
  }

  samples_required_lof <- post_process_predictions(gene_lof_fit)
  samples_required_mis <- post_process_predictions(gene_mis_fit)

  ####################################################################
  # Plot projections
  ####################################################################
  # Define the limits of the x-axis
  xlimits <- c(100, 1e8)

  # Create dataframe of ExAC, gnomAD v2, and gnomAD v4 sample sizes
  lines <- data.frame(
      intercepts = c(60706, 141456, 807162),
      names = c("1", "2", "3")
  )

  expected_projections <- function(projection_df, label = "pLoF") {
    # Generate plot displaying the percent of genes that would be expected to have a certain number or variants based on sample size
    # 'projection_df is the input dataframe with columns 'n_variants', 'n_required', and 'rank'
    # 'label' is the text that will be included at the top of the plot
    p <- projection_df %>%
      mutate(n_variants = forcats::fct_reorder(as.factor(n_variants), n_variants)) %>%
      ggplot() +
      aes(y = rank, x = n_required, color = n_variants) +
      geom_line(size = 2) +
      theme_classic() +
      scale_x_log10(labels = scales::comma, limits = xlimits) +
      xlab("Sample size required") +
      ylab("Percent of human genes") +
      geom_vline(xintercept = 60706, linetype = "dotted") + # Exac size
      geom_vline(xintercept = 141456, linetype = "dashed") + # gnomAD v2 size
      geom_vline(xintercept = 807162, linetype = "twodash") + # gnomAD v4 size

      # Add manual linetype scale for legend
      p <- p + scale_linetype_manual(
        name = "Database size",
        values = c("dotted", "dashed", "twodash"),
        labels = c("ExAC", "gnomADv2", "gnomADv4")
      ) + geom_vline(
        data = lines,
        aes(xintercept = intercepts, linetype = names),
        key_glyph = "path"
      ) + scale_color_discrete(name = ">= N variants\nexpected", h = c(40, 120)) + annotate(
        "text",
        x = xlimits[1],
        y = 1,
        hjust = 0,
        vjust = 1,
        label = label
      )

    return(p)
  }

  lof_projections <- expected_projections(samples_required_lof, "pLoF")
  missense_projections <- expected_projections(samples_required_mis, "Missense")

  ggsave(
    lof_projections,
    filename = paste("lof_ds_projections_", version, ".png", sep = ""),
    dpi = 300,
    width = 8,
    height = 6,
    units = "in"
  )
  ggsave(
    missense_projections,
    filename = paste("mis_ds_projections_", version, ".png", sep = ""),
    dpi = 300,
    width = 8,
    height = 6,
    units = "in"
  )
}
