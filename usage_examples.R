# =============================================================================
# MODULE SCORE CALCULATOR - USAGE EXAMPLES
# Demonstrating flexible column selection for different experimental designs
# =============================================================================

# Load directly from GitHub
source("https://raw.githubusercontent.com/Dragonmasterx87/MSCalculator/main/module_score_calculator.R")

# Load your data
counts <- read.csv("gene_rawCounts_Res.csv", stringsAsFactors = FALSE)

# Define gene sets
oxphos_genes <- c(
  "ENSMUSG00000029368", 
  "ENSMUSG00000064351", 
  "ENSMUSG00000032554"
)

# =============================================================================
# SCENARIO 1: Use ALL samples (4 controls + 4 treatment)
# =============================================================================
result_all <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = NULL,  # NULL = use all columns
  method = "zscore"
)

groups_all <- c(rep("Control", 4), rep("Treatment", 4))
plot_module_scores(result_all, groups_all, "OXPHOS - All Samples")


# =============================================================================
# SCENARIO 2: Use only CONTROL samples (columns 2-5)
# =============================================================================
result_controls <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = 2:5,  # Columns 2, 3, 4, 5 (mCon1-mCon4)
  method = "mean"
)

print("Control samples only:")
print(result_controls$module_scores)


# =============================================================================
# SCENARIO 3: Use only TREATMENT samples (columns 6-9)
# =============================================================================
result_treatment <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = 6:9,  # Columns 6, 7, 8, 9 (mRes1-mRes4)
  method = "mean"
)

print("Treatment samples only:")
print(result_treatment$module_scores)


# =============================================================================
# SCENARIO 4: Select specific samples by COLUMN NAME
# =============================================================================
result_selected <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = c("mCon1", "mCon3", "mRes1", "mRes3"),  # Specific names
  method = "zscore"
)

groups_selected <- c("Control", "Control", "Treatment", "Treatment")
plot_module_scores(result_selected, groups_selected, "OXPHOS - Selected Samples")


# =============================================================================
# SCENARIO 5: Mixed selection - 2 controls + 3 treatments
# =============================================================================
result_mixed <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = c(2, 4, 6, 7, 8),  # mCon1, mCon3, mRes1, mRes2, mRes3
  method = "zscore"
)

groups_mixed <- c("Control", "Control", "Treatment", "Treatment", "Treatment")
plot_module_scores(result_mixed, groups_mixed, "OXPHOS - Unbalanced Design")
test_module_scores(result_mixed, groups_mixed, test = "wilcox")


# =============================================================================
# SCENARIO 6: Compare multiple pathways with same sample selection
# =============================================================================
inflammation_genes <- c("ENSMUSG00000026103", "ENSMUSG00000020275")
fibrosis_genes <- c("ENSMUSG00000026728", "ENSMUSG00000055435")

# Use same sample selection for all pathways
selected_cols <- c(2, 3, 4, 6, 7, 8)  # 3 controls + 3 treatments
groups_comparison <- c(rep("Control", 3), rep("Treatment", 3))

oxphos_result <- calculate_module_score(counts, oxphos_genes, selected_cols, "zscore")
inflam_result <- calculate_module_score(counts, inflammation_genes, selected_cols, "zscore")
fibro_result <- calculate_module_score(counts, fibrosis_genes, selected_cols, "zscore")

# Combine results
combined_df <- data.frame(
  Sample = oxphos_result$samples_used,
  Group = groups_comparison,
  OXPHOS = oxphos_result$module_scores,
  Inflammation = inflam_result$module_scores,
  Fibrosis = fibro_result$module_scores
)

write.csv(combined_df, "multiple_pathway_scores.csv", row.names = FALSE)

# Statistical tests
cat("\n=== OXPHOS ===")
test_module_scores(oxphos_result, groups_comparison, "t.test")

cat("\n=== Inflammation ===")
test_module_scores(inflam_result, groups_comparison, "t.test")

cat("\n=== Fibrosis ===")
test_module_scores(fibro_result, groups_comparison, "t.test")


# =============================================================================
# SCENARIO 7: Time course or multiple conditions (ANOVA)
# =============================================================================
# If you have multiple conditions, e.g., columns 2-10 representing:
# Control (n=3), Low dose (n=3), High dose (n=3)

result_anova <- calculate_module_score(
  count_matrix = counts,
  gene_set = oxphos_genes,
  sample_columns = 2:9,  # All 8 samples
  method = "zscore"
)

# Define groups for ANOVA
groups_anova <- c(rep("Control", 4), rep("Treatment", 4))
# Or for 3 groups: groups_anova <- c(rep("Control", 3), rep("Low", 3), rep("High", 3))

test_module_scores(result_anova, groups_anova, test = "anova")


# =============================================================================
# SCENARIO 8: Batch-wise analysis (analyze batches separately)
# =============================================================================
# Batch 1: First 4 samples
batch1_result <- calculate_module_score(
  counts, oxphos_genes, 
  sample_columns = 2:5, 
  method = "mean"
)

# Batch 2: Last 4 samples  
batch2_result <- calculate_module_score(
  counts, oxphos_genes,
  sample_columns = 6:9,
  method = "mean"
)

cat("\nBatch 1 scores:", batch1_result$module_scores)
cat("\nBatch 2 scores:", batch2_result$module_scores)


# =============================================================================
# SCENARIO 9: Export individual pathway scores for downstream analysis
# =============================================================================
result_export <- calculate_module_score(
  counts, oxphos_genes,
  sample_columns = 2:9,
  method = "zscore"
)

groups_export <- c(rep("Control", 4), rep("Resmetirom", 4))

# Export using the new function
export_scores(result_export, groups_export, "OXPHOS", "oxphos_module_scores.csv")


# =============================================================================
# SCENARIO 10: Quality control - identify outlier samples
# =============================================================================
# Calculate scores with all samples
qc_result <- calculate_module_score(counts, oxphos_genes, method = "mean")

# Identify potential outliers
scores <- qc_result$module_scores
mean_score <- mean(scores)
sd_score <- sd(scores)
outliers <- abs(scores - mean_score) > 2 * sd_score

if (any(outliers)) {
  cat("\nPotential outlier samples detected:\n")
  print(names(scores)[outliers])

  # Re-run analysis without outliers
  non_outlier_cols <- which(!outliers) + 1  # +1 because gene_id is column 1
  result_clean <- calculate_module_score(
    counts, oxphos_genes,
    sample_columns = non_outlier_cols,
    method = "zscore"
  )
}
