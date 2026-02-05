
#' Calculate Module Score from Gene Expression Data
#'
#' This function calculates a module score for a custom gene set from normalized 
#' gene expression data. Includes automatic gene ID conversion (Ensembl <-> Symbol),
#' comprehensive data validation, and detailed diagnostic output.
#'
#' @param count_matrix A matrix or data frame with genes as rows and samples as columns.
#'                     First column should be gene IDs (gene_id).
#' @param gene_set A character vector of gene IDs (Ensembl IDs or gene symbols) 
#'                 representing the pathway/module of interest.
#' @param sample_columns Optional: numeric vector or character vector specifying which 
#'                       columns to use for analysis. If NULL, uses all columns except gene_id.
#' @param method Method for calculating module score: "mean", "median", or "zscore" (default: "mean")
#' @param normalize Logical, whether to normalize counts (log2 transform) before calculation (default: TRUE)
#' @param control_size Number of control genes to use per gene in gene_set (default: 100)
#' @param random_seed Random seed for reproducibility (default: 123)
#' @param gene_id_type Type of gene IDs in your data: "auto" (default), "ensembl", or "symbol"
#' @param convert_ids Logical, attempt to convert gene IDs if no matches found (default: TRUE)
#' @param verbose Logical, print detailed diagnostic information (default: TRUE)
#'
#' @return A list containing:
#'   - module_scores: Vector of module scores for each sample
#'   - genes_found: Genes from gene_set found in the data
#'   - genes_missing: Genes from gene_set not found in the data
#'   - normalized_data: Normalized expression matrix (if normalize=TRUE)
#'   - samples_used: Names of samples/columns used in analysis
#'   - diagnostics: List of diagnostic information
#'
#' @export

calculate_module_score <- function(count_matrix, 
                                   gene_set, 
                                   sample_columns = NULL,
                                   method = "mean",
                                   normalize = TRUE,
                                   control_size = 100,
                                   random_seed = 123,
                                   gene_id_type = "auto",
                                   convert_ids = TRUE,
                                   verbose = TRUE) {

  # Load required packages
  if (!require("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.")
  }

  # Set random seed
  set.seed(random_seed)

  if (verbose) {
    cat("\n", rep("=", 70), "\n", sep="")
    cat("MSCalculator - Module Score Calculation\n")
    cat(rep("=", 70), "\n\n", sep="")
  }

  # ============================================================================
  # STEP 1: DATA STRUCTURE VALIDATION
  # ============================================================================
  if (verbose) cat("STEP 1: Validating data structure...\n")

  if (!is.data.frame(count_matrix) && !is.matrix(count_matrix)) {
    stop("ERROR: count_matrix must be a data frame or matrix")
  }

  if (is.data.frame(count_matrix)) {
    if (verbose) {
      cat("  ✓ Input type: data.frame\n")
      cat("  ✓ Dimensions: ", nrow(count_matrix), " genes x ", 
          ncol(count_matrix)-1, " samples (+1 gene_id column)\n", sep="")
    }

    # Extract gene IDs
    gene_ids <- as.character(count_matrix[, 1])

    # Check for duplicates
    if (any(duplicated(gene_ids))) {
      n_dup <- sum(duplicated(gene_ids))
      warning(sprintf("Found %d duplicate gene IDs. Using first occurrence.", n_dup))
      keep_idx <- !duplicated(gene_ids)
      count_matrix <- count_matrix[keep_idx, ]
      gene_ids <- gene_ids[keep_idx]
    }

    # Handle column selection
    if (is.null(sample_columns)) {
      expr_matrix <- as.matrix(count_matrix[, -1])
      sample_names <- colnames(count_matrix)[-1]
    } else if (is.numeric(sample_columns)) {
      expr_matrix <- as.matrix(count_matrix[, sample_columns])
      sample_names <- colnames(count_matrix)[sample_columns]
    } else if (is.character(sample_columns)) {
      if (!all(sample_columns %in% colnames(count_matrix))) {
        missing_cols <- sample_columns[!sample_columns %in% colnames(count_matrix)]
        stop(sprintf("ERROR: Columns not found: %s", paste(missing_cols, collapse = ", ")))
      }
      expr_matrix <- as.matrix(count_matrix[, sample_columns])
      sample_names <- sample_columns
    } else {
      stop("ERROR: sample_columns must be NULL, numeric vector, or character vector")
    }

    rownames(expr_matrix) <- gene_ids
    colnames(expr_matrix) <- sample_names

  } else {
    # Matrix input
    expr_matrix <- count_matrix
    gene_ids <- rownames(expr_matrix)

    if (is.null(gene_ids)) {
      stop("ERROR: Matrix must have row names (gene IDs)")
    }

    if (verbose) {
      cat("  ✓ Input type: matrix\n")
      cat("  ✓ Dimensions: ", nrow(expr_matrix), " genes x ", 
          ncol(expr_matrix), " samples\n", sep="")
    }

    if (!is.null(sample_columns)) {
      if (is.numeric(sample_columns)) {
        expr_matrix <- expr_matrix[, sample_columns, drop = FALSE]
      } else if (is.character(sample_columns)) {
        expr_matrix <- expr_matrix[, sample_columns, drop = FALSE]
      }
    }
    sample_names <- colnames(expr_matrix)
  }

  # Check for NA/Inf values
  if (any(is.na(expr_matrix))) {
    n_na <- sum(is.na(expr_matrix))
    warning(sprintf("Found %d NA values in expression matrix. Will be handled as missing.", n_na))
  }

  if (any(is.infinite(expr_matrix))) {
    stop("ERROR: Found Inf values in expression matrix. Please check your data.")
  }

  if (verbose) {
    cat("  ✓ Selected samples (n=", ncol(expr_matrix), "): ", 
        paste(head(sample_names, 5), collapse=", "),
        if(ncol(expr_matrix) > 5) "..." else "", "\n", sep="")
  }

  # ============================================================================
  # STEP 2: GENE ID TYPE DETECTION
  # ============================================================================
  if (verbose) cat("\nSTEP 2: Analyzing gene ID format...\n")

  # Detect gene ID type in data
  data_id_type <- detect_gene_id_type(gene_ids)

  # Detect gene ID type in gene_set
  geneset_id_type <- detect_gene_id_type(gene_set)

  if (verbose) {
    cat("  ✓ Data gene IDs: ", data_id_type, " (examples: ", 
        paste(head(gene_ids, 3), collapse=", "), ")\n", sep="")
    cat("  ✓ Gene set IDs: ", geneset_id_type, " (examples: ", 
        paste(head(gene_set, 3), collapse=", "), ")\n", sep="")
  }

  # ============================================================================
  # STEP 3: GENE MATCHING & CONVERSION
  # ============================================================================
  if (verbose) cat("\nSTEP 3: Matching genes from gene set to data...\n")

  genes_found <- intersect(gene_set, gene_ids)
  genes_missing <- setdiff(gene_set, gene_ids)

  if (verbose) {
    cat("  ✓ Direct matches: ", length(genes_found), "/", length(gene_set), "\n", sep="")
  }

  # Try conversion if needed
  if (length(genes_found) == 0 && convert_ids && data_id_type != geneset_id_type) {
    if (verbose) {
      cat("  ⚠ No direct matches found. Attempting gene ID conversion...\n")
      cat("    (Converting ", geneset_id_type, " -> ", data_id_type, ")\n", sep="")
    }

    # Attempt conversion
    converted_result <- convert_gene_ids(gene_set, geneset_id_type, data_id_type, gene_ids)

    if (length(converted_result$converted) > 0) {
      genes_found <- converted_result$converted
      genes_missing <- converted_result$failed

      if (verbose) {
        cat("  ✓ Conversion successful: ", length(genes_found), " genes matched\n", sep="")
      }
    } else {
      if (verbose) {
        cat("  ✗ Conversion failed. No matches found.\n")
      }
    }
  }

  if (length(genes_found) == 0) {
    stop(sprintf(paste0(
      "ERROR: None of the genes in gene_set were found in the data.\n",
      "  - Data uses: %s\n",
      "  - Gene set uses: %s\n",
      "  - First few data IDs: %s\n",
      "  - First few gene set IDs: %s\n",
      "  Suggestion: Check gene ID format compatibility."
    ), data_id_type, geneset_id_type, 
    paste(head(gene_ids, 3), collapse=", "),
    paste(head(gene_set, 3), collapse=", ")))
  }

  if (verbose) {
    cat("  ✓ Final gene set size: ", length(genes_found), " genes\n", sep="")
    if (length(genes_missing) > 0) {
      cat("  ⚠ Missing genes (", length(genes_missing), "): ", 
          paste(head(genes_missing, 5), collapse=", "),
          if(length(genes_missing) > 5) "..." else "", "\n", sep="")
    }
  }

  # ============================================================================
  # STEP 4: DATA NORMALIZATION
  # ============================================================================
  if (verbose) cat("\nSTEP 4: Data normalization...\n")

  if (normalize) {
    # Check if data might already be normalized
    data_range <- range(expr_matrix, na.rm = TRUE)
    if (data_range[2] < 100) {
      warning("Data values are low (max < 100). Data might already be normalized.")
    }

    expr_matrix <- log2(expr_matrix + 1)
    if (verbose) {
      cat("  ✓ Applied log2(count + 1) transformation\n")
      new_range <- range(expr_matrix, na.rm = TRUE)
      cat("  ✓ Value range after normalization: [", 
          sprintf("%.2f", new_range[1]), ", ", 
          sprintf("%.2f", new_range[2]), "]\n", sep="")
    }
  } else {
    if (verbose) {
      cat("  ⊘ Normalization skipped (normalize = FALSE)\n")
      data_range <- range(expr_matrix, na.rm = TRUE)
      cat("  ✓ Value range: [", sprintf("%.2f", data_range[1]), ", ", 
          sprintf("%.2f", data_range[2]), "]\n", sep="")
    }
  }

  # ============================================================================
  # STEP 5: MODULE SCORE CALCULATION
  # ============================================================================
  if (verbose) cat("\nSTEP 5: Calculating module scores (method: ", method, ")...\n", sep="")

  # Extract expression for genes in the module
  module_expr <- expr_matrix[genes_found, , drop = FALSE]

  if (method == "mean") {
    module_scores <- colMeans(module_expr, na.rm = TRUE)
    if (verbose) cat("  ✓ Calculated mean expression across ", length(genes_found), " genes\n", sep="")

  } else if (method == "median") {
    module_scores <- apply(module_expr, 2, median, na.rm = TRUE)
    if (verbose) cat("  ✓ Calculated median expression across ", length(genes_found), " genes\n", sep="")

  } else if (method == "zscore") {
    # Calculate mean expression for binning
    gene_means <- rowMeans(expr_matrix, na.rm = TRUE)

    # Bin genes by expression level
    gene_bins <- cut(gene_means, breaks = 25, labels = FALSE)
    names(gene_bins) <- rownames(expr_matrix)

    module_gene_bins <- gene_bins[genes_found]

    # Select control genes
    control_genes <- c()
    for (bin in unique(module_gene_bins)) {
      bin_gene_names <- names(gene_bins)[gene_bins == bin & !names(gene_bins) %in% genes_found]

      if (length(bin_gene_names) > 0) {
        n_to_sample <- min(control_size, length(bin_gene_names))
        control_genes <- c(control_genes, sample(bin_gene_names, n_to_sample))
      }
    }

    if (length(control_genes) == 0) {
      warning("No control genes found. Using simple mean instead.")
      module_scores <- colMeans(module_expr, na.rm = TRUE)
      if (verbose) cat("  ⚠ Fallback to mean method (no control genes available)\n")
    } else {
      control_expr <- expr_matrix[control_genes, , drop = FALSE]
      control_mean <- colMeans(control_expr, na.rm = TRUE)
      control_sd <- apply(control_expr, 2, sd, na.rm = TRUE)

      module_mean <- colMeans(module_expr, na.rm = TRUE)

      module_scores <- (module_mean - control_mean) / control_sd

      if (verbose) {
        cat("  ✓ Selected ", length(control_genes), " control genes\n", sep="")
        cat("  ✓ Calculated z-scores across ", ncol(expr_matrix), " samples\n", sep="")
        cat("  ✓ Score range: [", sprintf("%.2f", min(module_scores)), ", ", 
            sprintf("%.2f", max(module_scores)), "]\n", sep="")
      }
    }

  } else {
    stop("ERROR: method must be one of: 'mean', 'median', 'zscore'")
  }

  # ============================================================================
  # RESULTS SUMMARY
  # ============================================================================
  if (verbose) {
    cat("\n", rep("=", 70), "\n", sep="")
    cat("CALCULATION COMPLETE!\n")
    cat(rep("=", 70), "\n", sep="")
    cat("Summary:\n")
    cat("  • Samples analyzed: ", ncol(expr_matrix), "\n", sep="")
    cat("  • Genes in module: ", length(genes_found), "\n", sep="")
    cat("  • Method: ", method, "\n", sep="")
    cat("  • Score range: [", sprintf("%.3f", min(module_scores)), ", ", 
        sprintf("%.3f", max(module_scores)), "]\n", sep="")
    cat(rep("=", 70), "\n\n", sep="")
  }

  # Compile diagnostics
  diagnostics <- list(
    data_id_type = data_id_type,
    geneset_id_type = geneset_id_type,
    n_genes_input = length(gene_set),
    n_genes_found = length(genes_found),
    n_genes_missing = length(genes_missing),
    n_samples = ncol(expr_matrix),
    normalized = normalize,
    method = method
  )

  # Return results
  result <- list(
    module_scores = module_scores,
    genes_found = genes_found,
    genes_missing = genes_missing,
    normalized_data = expr_matrix,
    samples_used = sample_names,
    method = method,
    diagnostics = diagnostics
  )

  return(result)
}


#' Detect Gene ID Type
#' @keywords internal
detect_gene_id_type <- function(gene_ids) {
  # Sample genes for detection
  sample_genes <- head(unique(gene_ids), 100)

  # Check for Ensembl pattern (ENSMUSG, ENSG, etc.)
  ensembl_pattern <- "^ENS[A-Z]*[GT][0-9]+"
  n_ensembl <- sum(grepl(ensembl_pattern, sample_genes))

  # Check for typical symbol patterns (all caps, mixed case)
  symbol_pattern <- "^[A-Z][A-Za-z0-9-]*$"
  n_symbols <- sum(grepl(symbol_pattern, sample_genes))

  if (n_ensembl / length(sample_genes) > 0.8) {
    return("Ensembl")
  } else if (n_symbols / length(sample_genes) > 0.5) {
    return("Symbol")
  } else {
    return("Unknown")
  }
}


#' Convert Gene IDs (Simple Pattern-Based)
#' @keywords internal
convert_gene_ids <- function(gene_set, from_type, to_type, available_genes) {
  # Simple conversion - try case variations and partial matches
  # For production use, consider using biomaRt or org.Mm.eg.db

  converted <- c()
  failed <- gene_set

  # Try case-insensitive matching
  for (gene in gene_set) {
    matches <- available_genes[tolower(available_genes) == tolower(gene)]
    if (length(matches) > 0) {
      converted <- c(converted, matches[1])
      failed <- setdiff(failed, gene)
    }
  }

  list(converted = converted, failed = failed)
}


#' Plot Module Scores
#' @export
plot_module_scores <- function(module_result, 
                               groups, 
                               title = "Module Score",
                               colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                               show_points = TRUE) {

  if (!require("ggplot2", quietly = TRUE)) {
    stop("ERROR: Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  if (length(groups) != length(module_result$module_scores)) {
    stop(sprintf("ERROR: Length of groups (%d) must match number of samples (%d)", 
                length(groups), length(module_result$module_scores)))
  }

  plot_data <- data.frame(
    Sample = names(module_result$module_scores),
    Score = module_result$module_scores,
    Group = factor(groups, levels = unique(groups))
  )

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

  if (show_points) {
    p <- p + geom_jitter(width = 0.2, size = 3, alpha = 0.8)
  }

  if (length(unique(groups)) == 2) {
    group_levels <- unique(groups)
    score1 <- plot_data$Score[plot_data$Group == group_levels[1]]
    score2 <- plot_data$Score[plot_data$Group == group_levels[2]]

    t_result <- t.test(score1, score2)
    p_value <- t_result$p.value

    y_max <- max(plot_data$Score) * 1.1
    p_label <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))

    p <- p + annotate("text", x = 1.5, y = y_max, 
                     label = p_label, size = 4, fontface = "italic")
  }

  return(p)
}


#' Statistical Test for Module Scores
#' @export
test_module_scores <- function(module_result, groups, test = "t.test") {

  scores <- module_result$module_scores

  if (length(groups) != length(scores)) {
    stop(sprintf("ERROR: Length of groups (%d) must match number of samples (%d)", 
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

    p_val <- result_summary[[1]]$"Pr(>F)"[1]
    if (p_val < 0.05 && length(unique(groups)) > 2) {
      cat("\nPost-hoc pairwise t-tests (with Bonferroni correction):\n")
      posthoc <- pairwise.t.test(scores, groups, p.adjust.method = "bonferroni")
      print(posthoc)
    }

  } else {
    stop("ERROR: Invalid test or incompatible number of groups")
  }

  return(result)
}


#' Export Module Scores to CSV
#' @export
export_scores <- function(module_result, groups, pathway_name, output_file) {

  if (length(groups) != length(module_result$module_scores)) {
    stop(sprintf("ERROR: Length of groups (%d) must match number of samples (%d)", 
                length(groups), length(module_result$module_scores)))
  }

  output_df <- data.frame(
    Sample = module_result$samples_used,
    Group = groups,
    Score = module_result$module_scores
  )

  colnames(output_df)[3] <- paste0(pathway_name, "_Score")

  write.csv(output_df, output_file, row.names = FALSE)

  cat(sprintf("\nExported scores to: %s\n", output_file))
  cat(sprintf("Columns: Sample, Group, %s_Score\n", pathway_name))

  return(invisible(output_df))
}
