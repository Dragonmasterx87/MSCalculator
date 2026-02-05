# Module Score Calculator

A flexible R function to calculate pathway/gene set module scores from gene expression data with **powerful column selection** for any experimental design.

<p align="center">
  <img src="generated-image.png" width="200"/>
</p>

## 🎯 Key Features

- ✅ **Flexible column selection**: Choose any combination of samples by index or name
- ✅ **Multiple scoring methods**: mean, median, or z-score with control genes
- ✅ **Works with raw counts or normalized data**
- ✅ **Automatic normalization** (log2 transform)
- ✅ **Missing gene reporting**
- ✅ **Statistical testing & visualization**
- ✅ **CSV export functionality**

## Installation

Simply source the R script:

```r
source("module_score_calculator.R")
```

No package installation required! Only base R dependencies: `stats`, and optionally `ggplot2` for plotting.

## Quick Start

```r
# Load your data
raw_counts <- read.csv("gene_rawCounts_Res.csv", stringsAsFactors = FALSE)

# Define your gene set
oxphos_genes <- c("ENSMUSG00000029368", "ENSMUSG00000064351", "ENSMUSG00000032554")

# Calculate module score (use all samples)
result <- calculate_module_score(
  count_matrix = raw_counts,
  gene_set = oxphos_genes,
  method = "zscore"
)

# View results
print(result$module_scores)
```

## 🔥 Column Selection - The Power Feature

The `sample_columns` parameter lets you analyze **any combination** of samples:

### Example 1: Use all samples
```r
result <- calculate_module_score(counts, genes, sample_columns = NULL)
```

### Example 2: Select by column index
```r
# Use only columns 2-5 (e.g., control samples)
result <- calculate_module_score(counts, genes, sample_columns = 2:5)

# Use specific non-consecutive columns
result <- calculate_module_score(counts, genes, sample_columns = c(2, 4, 6, 8))
```

### Example 3: Select by column name
```r
result <- calculate_module_score(
  counts, genes, 
  sample_columns = c("mCon1", "mCon2", "mRes1", "mRes2")
)
```

### Example 4: Flexible experimental designs
```r
# 4 controls + 5 treatments
result <- calculate_module_score(counts, genes, sample_columns = 2:10)
groups <- c(rep("Control", 4), rep("Treatment", 5))

# 3 controls + 3 treatments (skip some samples)
result <- calculate_module_score(counts, genes, sample_columns = c(2:4, 6:8))
groups <- c(rep("Control", 3), rep("Treatment", 3))

# Unbalanced design: 2 controls + 4 treatments
result <- calculate_module_score(counts, genes, sample_columns = c(2, 3, 6, 7, 8, 9))
groups <- c(rep("Control", 2), rep("Treatment", 4))
```

## Complete Workflow Example

```r
source("module_score_calculator.R")

# Load data
counts <- read.csv("gene_rawCounts_Res.csv", stringsAsFactors = FALSE)

# Define gene sets
oxphos_genes <- c("ENSMUSG00000029368", "ENSMUSG00000064351")
inflammation_genes <- c("ENSMUSG00000026103", "ENSMUSG00000020275")

# Select specific samples (3 controls + 3 treatments)
selected_samples <- c(2, 3, 4, 6, 7, 8)

# Calculate scores for multiple pathways using SAME samples
oxphos_result <- calculate_module_score(counts, oxphos_genes, selected_samples, "zscore")
inflam_result <- calculate_module_score(counts, inflammation_genes, selected_samples, "zscore")

# Define groups
groups <- c(rep("Control", 3), rep("Treatment", 3))

# Plot
plot_module_scores(oxphos_result, groups, "OXPHOS Score")
plot_module_scores(inflam_result, groups, "Inflammation Score")

# Statistical tests
test_module_scores(oxphos_result, groups, "t.test")
test_module_scores(inflam_result, groups, "t.test")

# Export results
export_scores(oxphos_result, groups, "OXPHOS", "oxphos_scores.csv")
```

## Main Functions

### `calculate_module_score()`

**Parameters:**
- `count_matrix`: Data frame/matrix with genes as rows, samples as columns (first column = gene IDs)
- `gene_set`: Character vector of gene IDs (Ensembl or symbols)
- `sample_columns`: **[NEW]** Numeric/character vector specifying which columns to analyze (default: NULL = all)
- `method`: Scoring method - `"mean"`, `"median"`, or `"zscore"` (default: `"mean"`)
- `normalize`: Log2 normalize raw counts (default: `TRUE`)
- `control_size`: Number of control genes for zscore method (default: `100`)
- `random_seed`: Reproducibility seed (default: `123`)

**Returns:**
- `module_scores`: Score for each sample
- `genes_found`: Genes found in your data
- `genes_missing`: Genes not found
- `normalized_data`: Normalized matrix
- `samples_used`: Names of analyzed samples

### `plot_module_scores()`

**Parameters:**
- `module_result`: Output from `calculate_module_score()`
- `groups`: Group labels (must match number of samples analyzed)
- `title`: Plot title
- `colors`: Custom color palette
- `show_points`: Show individual points (default: `TRUE`)

**Requires:** `ggplot2`

### `test_module_scores()`

**Parameters:**
- `module_result`: Output from `calculate_module_score()`
- `groups`: Group labels
- `test`: `"t.test"`, `"wilcox"`, or `"anova"`

### `export_scores()`

**Parameters:**
- `module_result`: Output from `calculate_module_score()`
- `groups`: Group labels
- `pathway_name`: Name for the score column
- `output_file`: CSV file path

## Scoring Methods

### 1. Mean (default)
```r
result <- calculate_module_score(counts, genes, method = "mean")
```
- Simple average of gene expressions
- Fast and interpretable
- Best for: Quick comparisons, well-defined pathways

### 2. Median
```r
result <- calculate_module_score(counts, genes, method = "median")
```
- Robust to outliers
- Best for: Data with extreme values

### 3. Z-score (recommended for publication)
```r
result <- calculate_module_score(counts, genes, method = "zscore")
```
- Accounts for background expression using control genes
- Control genes matched by expression level
- Similar to Seurat's `AddModuleScore`
- Best for: Rigorous analysis, publications
- Formula: `(module_mean - control_mean) / control_sd`

## Real-World Use Cases

See `usage_examples.R` for 10 complete scenarios including:

1. **Standard comparison**: 4 controls vs 4 treatments
2. **Subset analysis**: Analyze only controls or only treatments
3. **Dose response**: Multiple treatment groups for ANOVA
4. **Batch analysis**: Separate analysis of different batches
5. **Outlier handling**: QC and re-analysis without outliers
6. **Multiple pathways**: Compare OXPHOS, inflammation, fibrosis, etc.
7. **Unbalanced designs**: 2 controls vs 5 treatments
8. **Time course**: Multiple time points

## Data Format

Your input should have:
- **First column**: Gene IDs (Ensembl IDs or gene symbols)
- **Other columns**: Sample expression values

Example:
```
gene_id       | mCon1  | mCon2  | mCon3  | mCon4  | mRes1  | mRes2  | mRes3  | mRes4
------------- | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------
ENSMUSG000... | 400984 | 499797 | 393084 | 384787 | 459560 | 593235 | 486429 | 461127
ENSMUSG000... | 268658 | 249360 | 267278 | 232033 | 284437 | 231158 | 262000 | 278891
```

## Getting Gene Sets

### From Enrichment Analysis
```r
# After GO/KEGG enrichment with clusterProfiler
oxphos_genes <- enrichment_result@result$geneID[
  enrichment_result@result$Description == "oxidative phosphorylation"
]
oxphos_genes <- unlist(strsplit(oxphos_genes, "/"))
```

### From Databases
- **MSigDB**: http://www.gsea-msigdb.org/gsea/msigdb/
- **GO/KEGG**: Via `clusterProfiler` package
- **Reactome**: Pathway database

### Custom Lists
```r
my_genes <- c("ENSMUSG00000029368", "ENSMUSG00000064351", ...)
```

## Tips & Best Practices

1. **Column Selection Strategy**
   - Use `sample_columns = NULL` to include all samples initially
   - Use specific indices to exclude outliers or focus on subsets
   - Use column names for clarity in complex designs

2. **Choosing a Method**
   - Start with `"mean"` for exploration
   - Use `"zscore"` for final analysis/publications
   - Use `"median"` if you suspect outliers

3. **Normalization**
   - `normalize = TRUE` for raw counts
   - `normalize = FALSE` for TPM, RPKM, FPKM, or pre-normalized data

4. **Sample Size**
   - Minimum 3 samples per group recommended
   - Unbalanced designs are OK (e.g., 2 vs 4 samples)

5. **Gene Set Size**
   - Minimum 5-10 genes recommended
   - Larger sets (20-100 genes) are more robust
   - Very small sets (<5 genes) may be unstable

## Troubleshooting

**"None of the genes found"**
- Check gene ID format matches your data (Ensembl vs symbols)
- Verify spelling and species (mouse vs human)

**"Columns not found"**
- When using column names, check exact spelling
- Use `colnames(counts)` to see available names

**"Length of groups doesn't match"**
- Ensure groups vector length matches selected samples
- If selecting 6 columns, groups must have length 6

**Low p-values but small differences**
- Normal with small sample sizes
- Report effect sizes (mean difference) alongside p-values

## Output Files

```r
# Export single pathway
export_scores(result, groups, "OXPHOS", "oxphos_scores.csv")

# Creates CSV with columns:
# Sample, Group, OXPHOS_Score
```

## Citation

If you use this in publications:
```
Qadir, M.M.F. (2026). Module Score Calculator: Flexible pathway scoring 
for gene expression data. GitHub: https://github.com/[username]/module-score-calculator
```

## License

MIT License - Free to use, modify, and distribute.

## Contact

Questions or suggestions? Open an issue on GitHub!

## Acknowledgments

- Inspired by Seurat's `AddModuleScore`
- GSEA methodology
- Community feedback from collaborative projects
