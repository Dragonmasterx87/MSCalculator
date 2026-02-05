
#' Calculate Module Score from Gene Expression Data
#'
#' This function calculates a module score for a custom gene set from normalized 
#' gene expression data. The module score represents the average expression of 
#' genes in the pathway relative to a control gene set.
#'
#' @param count_matrix A matrix or data frame with genes as rows and samples as columns.
#'                     First column should be gene IDs (gene_id).
#' @param gene_set A character vector of gene IDs (Ensembl IDs or gene symbols) 
#'                 representing the pathway/module of interest.
#' @param sample_columns Optional: numeric vector or character vector specifying which 
#'                       columns to use for analysis. If NULL, uses all columns except gene_id.
#'                       Examples: c(2:5) or c("mCon1", "mCon2", "mRes1", "mRes2")
#' @param method Method for calculating module score: "mean", "median", or "zscore" (default: "mean")
#' @param normalize Logical, whether to normalize counts (log2 transform) before calculation (default: TRUE)
#' @param control_size Number of control genes to use per gene in gene_set (default: 100)
#' @param random_seed Random seed for reproducibility (default: 123)
#'
#' @return A list containing:
#'   - module_scores: Vector of module scores for each sample
#'   - genes_found: Genes from gene_set found in the data
#'   - genes_missing: Genes from gene_set not found in the data
#'   - normalized_data: Normalized expression matrix (if normalize=TRUE)
#'   - samples_used: Names of samples/columns used in analysis
#'
#' @examples
#' # Example 1: Use all samples
#' raw_counts <- read.csv("gene_rawCounts_Res.csv", stringsAsFactors = FALSE)
#' oxphos_genes <- c("ENSMUSG00000029368", "ENSMUSG00000064351")
#' result <- calculate_module_score(raw_counts, oxphos_genes)
#' 
#' # Example 2: Select specific columns by index (columns 2-5 = first 4 samples)
#' result <- calculate_module_score(raw_counts, oxphos_genes, sample_columns = 2:5)
#' 
#' # Example 3: Select specific columns by name
#' result <- calculate_module_score(raw_counts, oxphos_genes, 
#'                                  sample_columns = c("mCon1", "mCon2", "mRes1"))
#' 
#' # Example 4: Mix control and treatment samples
#' result <- calculate_module_score(raw_counts, oxphos_genes, 
#'                                  sample_columns = c(2, 3, 6, 7))
#'
#' @export

calculate_module_score <- function(count_matrix, 
                                   gene_set, 
                                   sample_columns = NULL,
                                   method = "mean",
                                   normalize = TRUE,
                                   control_size = 100,
                                   random_seed = 123) {

  # Load required packages
  if (!require("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.")
  }

  # Set random seed for reproducibility
  set.seed(random_seed)

  # Extract gene IDs and expression matrix
  if (is.data.frame(count_matrix)) {
    gene_ids <- count_matrix[, 1]

    # Handle column selection
    if (is.null(sample_columns)) {
      # Use all columns except the first (gene_id)
      expr_matrix <- as.matrix(count_matrix[, -1])
      sample_names <- colnames(count_matrix)[-1]
    } else if (is.numeric(sample_columns)) {
      # Select by numeric index
      expr_matrix <- as.matrix(count_matrix[, sample_columns])
      sample_names <- colnames(count_matrix)[sample_columns]
    } else if (is.character(sample_columns)) {
      # Select by column name
      if (!all(sample_columns %in% colnames(count_matrix))) {
        missing_cols <- sample_columns[!sample_columns %in% colnames(count_matrix)]
        stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
      }
      expr_matrix <- as.matrix(count_matrix[, sample_columns])
      sample_names <- sample_columns
    } else {
      stop("sample_columns must be NULL, numeric vector, or character vector")
    }

    rownames(expr_matrix) <- gene_ids
    colnames(expr_matrix) <- sample_names

  } else if (is.matrix(count_matrix)) {
    expr_matrix <- count_matrix
    gene_ids <- rownames(expr_matrix)

    # Handle column selection for matrix input
    if (!is.null(sample_columns)) {
      if (is.numeric(sample_columns)) {
        expr_matrix <- expr_matrix[, sample_columns, drop = FALSE]
      } else if (is.character(sample_columns)) {
        expr_matrix <- expr_matrix[, sample_columns, drop = FALSE]
      }
    }
    sample_names <- colnames(expr_matrix)

  } else {
    stop("count_matrix must be a data frame or matrix")
  }

  cat(sprintf("Using %d samples for analysis: %s\n", 
              ncol(expr_matrix), 
              paste(sample_names, collapse = ", ")))

  # Normalize if requested
  if (normalize) {
    # Log2 normalization (add pseudocount to avoid log(0))
    expr_matrix <- log2(expr_matrix + 1)
    cat("Data normalized using log2(count + 1) transformation\n")
  }

  # Check which genes are present
  genes_found <- intersect(gene_set, gene_ids)
  genes_missing <- setdiff(gene_set, gene_ids)

  if (length(genes_found) == 0) {
    stop("None of the genes in gene_set were found in the count_matrix")
  }

  cat(sprintf("Found %d/%d genes from gene_set in the data\n", 
              length(genes_found), length(gene_set)))

  if (length(genes_missing) > 0) {
    cat(sprintf("Missing %d genes: %s\n", 
                length(genes_missing), 
                paste(head(genes_missing, 5), collapse = ", ")))
  }

  # Extract expression for genes in the module
  module_expr <- expr_matrix[genes_found, , drop = FALSE]

  # Calculate module score based on method
  if (method == "mean") {
    # Simple mean of gene expressions
    module_scores <- colMeans(module_expr, na.rm = TRUE)

  } else if (method == "median") {
    # Median of gene expressions
    module_scores <- apply(module_expr, 2, median, na.rm = TRUE)

  } else if (method == "zscore") {
    # Z-score normalization with control genes
    # Calculate mean expression for each gene across all samples
    gene_means <- rowMeans(expr_matrix, na.rm = TRUE)

    # Bin genes by expression level
    gene_bins <- cut(gene_means, breaks = 25, labels = FALSE)
    module_gene_bins <- gene_bins[genes_found]

    # Select control genes from same expression bins
    control_genes <- c()
    for (bin in unique(module_gene_bins)) {
      bin_genes <- gene_ids[gene_bins == bin & !gene_ids %in% genes_found]
      if (length(bin_genes) > 0) {
        n_to_sample <- min(control_size, length(bin_genes))
        control_genes <- c(control_genes, sample(bin_genes, n_to_sample))
      }
    }

    if (length(control_genes) == 0) {
      warning("No control genes found, using simple mean instead")
      module_scores <- colMeans(module_expr, na.rm = TRUE)
    } else {
      # Calculate control mean
      control_expr <- expr_matrix[control_genes, , drop = FALSE]
      control_mean <- colMeans(control_expr, na.rm = TRUE)
      control_sd <- apply(control_expr, 2, sd, na.rm = TRUE)

      # Calculate module mean
      module_mean <- colMeans(module_expr, na.rm = TRUE)

      # Z-score: (module_mean - control_mean) / control_sd
      module_scores <- (module_mean - control_mean) / control_sd
      cat(sprintf("Used %d control genes for z-score calculation\n", length(control_genes)))
    }

  } else {
    stop("method must be one of: 'mean', 'median', 'zscore'")
  }

  # Return results
  result <- list(
    module_scores = module_scores,
    genes_found = genes_found,
    genes_missing = genes_missing,
    normalized_data = expr_matrix,
    samples_used = sample_names,
    method = method
  )

  cat("\nModule score calculation complete!\n")
  return(result)
}


#' Plot Module Scores
#'
#' Visualize module scores across samples with group comparison
#'
#' @param module_result Output from calculate_module_score()
#' @param groups Character vector indicating group membership for each sample.
#'               Must match the length of samples in module_result.
#' @param title Plot title (default: "Module Score")
#' @param colors Vector of colors for groups (default: c("#E69F00", "#56B4E9"))
#' @param show_points Logical, whether to show individual points (default: TRUE)
#'
#' @return A ggplot object
#'
#' @examples
#' # After calculating module scores with specific samples
#' result <- calculate_module_score(counts, genes, sample_columns = c(2:5, 7:9))
#' groups <- c(rep("Control", 4), rep("Treatment", 3))
#' plot_module_scores(result, groups, title = "OXPHOS Module Score")
#'
#' @export

plot_module_scores <- function(module_result, 
                               groups, 
                               title = "Module Score",
                               colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                               show_points = TRUE) {

  # Check if ggplot2 is available
  if (!require("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  # Validate groups length
  if (length(groups) != length(module_result$module_scores)) {
    stop(sprintf("Length of groups (%d) must match number of samples (%d)", 
                length(groups), length(module_result$module_scores)))
  }

  # Create data frame for plotting
  plot_data <- data.frame(
    Sample = names(module_result$module_scores),
    Score = module_result$module_scores,
    Group = factor(groups, levels = unique(groups))
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = colors) +
    labs(title = title, 
         y = paste("Module Score (", module_result$method, ")", sep = ""),
         x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  # Add points if requested
  if (show_points) {
    p <- p + geom_jitter(width = 0.2, size = 3, alpha = 0.8)
  }

  # Add statistical comparison if two groups
  if (length(unique(groups)) == 2) {
    group_levels <- unique(groups)
    score1 <- plot_data$Score[plot_data$Group == group_levels[1]]
    score2 <- plot_data$Score[plot_data$Group == group_levels[2]]

    # T-test
    t_result <- t.test(score1, score2)
    p_value <- t_result$p.value

    # Add p-value annotation
    y_max <- max(plot_data$Score) * 1.1
    p_label <- ifelse(p_value < 0.001, "p < 0.001", 
                     sprintf("p = %.3f", p_value))

    p <- p + annotate("text", x = 1.5, y = y_max, 
                     label = p_label, size = 4, fontface = "italic")
  }

  return(p)
}


#' Statistical Test for Module Scores
#'
#' Perform statistical tests comparing module scores between groups
#'
#' @param module_result Output from calculate_module_score()
#' @param groups Character vector indicating group membership for each sample
#' @param test Type of test: "t.test", "wilcox", or "anova" (default: "t.test")
#'
#' @return A list with test results
#'
#' @export

test_module_scores <- function(module_result, groups, test = "t.test") {

  scores <- module_result$module_scores

  # Validate groups length
  if (length(groups) != length(scores)) {
    stop(sprintf("Length of groups (%d) must match number of samples (%d)", 
                length(groups), length(scores)))
  }

  if (test == "t.test" && length(unique(groups)) == 2) {
    group_levels <- unique(groups)
    score1 <- scores[groups == group_levels[1]]
    score2 <- scores[groups == group_levels[2]]

    result <- t.test(score1, score2)

    cat(sprintf("\nTwo-sample t-test:\n"))
    cat(sprintf("  %s mean: %.4f (n=%d)\n", group_levels[1], mean(score1), length(score1)))
    cat(sprintf("  %s mean: %.4f (n=%d)\n", group_levels[2], mean(score2), length(score2)))
    cat(sprintf("  Difference: %.4f\n", mean(score1) - mean(score2)))
    cat(sprintf("  t-statistic: %.4f\n", result$statistic))
    cat(sprintf("  p-value: %.4e\n", result$p.value))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", result$conf.int[1], result$conf.int[2]))

  } else if (test == "wilcox" && length(unique(groups)) == 2) {
    group_levels <- unique(groups)
    score1 <- scores[groups == group_levels[1]]
    score2 <- scores[groups == group_levels[2]]

    result <- wilcox.test(score1, score2)

    cat(sprintf("\nWilcoxon rank-sum test:\n"))
    cat(sprintf("  %s median: %.4f (n=%d)\n", group_levels[1], median(score1), length(score1)))
    cat(sprintf("  %s median: %.4f (n=%d)\n", group_levels[2], median(score2), length(score2)))
    cat(sprintf("  W-statistic: %.4f\n", result$statistic))
    cat(sprintf("  p-value: %.4e\n", result$p.value))

  } else if (test == "anova") {
    result <- aov(scores ~ groups)
    result_summary <- summary(result)

    cat("\nANOVA results:\n")
    print(result_summary)

    # Post-hoc if significant
    p_val <- result_summary[[1]]$"Pr(>F)"[1]
    if (p_val < 0.05 && length(unique(groups)) > 2) {
      cat("\nPost-hoc pairwise t-tests (with Bonferroni correction):\n")
      posthoc <- pairwise.t.test(scores, groups, p.adjust.method = "bonferroni")
      print(posthoc)
    }

  } else {
    stop("Invalid test or incompatible number of groups")
  }

  return(result)
}


#' Export Module Scores to CSV
#'
#' Export module scores with sample information and group labels
#'
#' @param module_result Output from calculate_module_score()
#' @param groups Character vector indicating group membership for each sample
#' @param pathway_name Name of the pathway/module (for column naming)
#' @param output_file Path to output CSV file
#'
#' @return Invisibly returns the data frame that was exported
#'
#' @examples
#' result <- calculate_module_score(counts, oxphos_genes, sample_columns = 2:9)
#' groups <- c(rep("Control", 4), rep("Treatment", 4))
#' export_scores(result, groups, "OXPHOS", "oxphos_scores.csv")
#'
#' @export

export_scores <- function(module_result, groups, pathway_name, output_file) {

  # Validate groups length
  if (length(groups) != length(module_result$module_scores)) {
    stop(sprintf("Length of groups (%d) must match number of samples (%d)", 
                length(groups), length(module_result$module_scores)))
  }

  # Create output data frame
  output_df <- data.frame(
    Sample = module_result$samples_used,
    Group = groups,
    Score = module_result$module_scores
  )

  colnames(output_df)[3] <- paste0(pathway_name, "_Score")

  # Write to file
  write.csv(output_df, output_file, row.names = FALSE)

  cat(sprintf("\nExported scores to: %s\n", output_file))
  cat(sprintf("Columns: Sample, Group, %s_Score\n", pathway_name))

  return(invisible(output_df))
}
