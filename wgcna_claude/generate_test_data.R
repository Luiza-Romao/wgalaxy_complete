#!/usr/bin/env Rscript
# =============================================================================
# Generate minimal test data for Galaxy WGCNA tool testing
# Creates:
#   test-data/test_vst.tabular    — small VST matrix (300 genes × 12 samples)
#   test-data/test_samples.tabular — sample metadata
#   test-data/test_up_degs.tabular — up-regulated gene list (optional)
#   test-data/test_down_degs.tabular — down-regulated gene list (optional)
# =============================================================================

set.seed(42)
dir.create("test-data", showWarnings = FALSE)

n_genes   <- 300
n_samples <- 12
n_induced <- 6
n_control <- 6

sample_ids <- c(
  paste0("Induced_", seq_len(n_induced)),
  paste0("Control_", seq_len(n_control))
)
gene_ids <- paste0("Gene_", seq_len(n_genes))

# Simulate expression: induced genes up in Induced samples
expr_matrix <- matrix(
  rnorm(n_genes * n_samples, mean = 8, sd = 2),
  nrow   = n_genes,
  ncol   = n_samples,
  dimnames = list(gene_ids, sample_ids)
)

# Introduce 3 artificial co-expression modules
# Module A: genes 1-80 up-regulated in Induced
module_a <- 1:80
expr_matrix[module_a, 1:n_induced] <-
  expr_matrix[module_a, 1:n_induced] + 3 +
  matrix(rnorm(length(module_a) * n_induced, 0, 0.5),
         nrow = length(module_a))

# Module B: genes 81-160 down-regulated in Induced
module_b <- 81:160
expr_matrix[module_b, 1:n_induced] <-
  expr_matrix[module_b, 1:n_induced] - 3 +
  matrix(rnorm(length(module_b) * n_induced, 0, 0.5),
         nrow = length(module_b))

# Module C: genes 161-240 mildly correlated
module_c <- 161:240
base_signal <- rnorm(n_samples, 0, 1.5)
expr_matrix[module_c, ] <- expr_matrix[module_c, ] +
  matrix(rep(base_signal, each = length(module_c)),
         nrow = length(module_c))

# Ensure all values are positive (VST is positive)
expr_matrix <- expr_matrix - min(expr_matrix) + 5

# ---- Write VST matrix ----
vst_df <- as.data.frame(expr_matrix)
vst_df <- cbind(GeneID = rownames(vst_df), vst_df)
write.table(vst_df,
  file      = "test-data/test_vst.tabular",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
cat("Written: test-data/test_vst.tabular\n")

# ---- Write sample metadata ----
sample_info <- data.frame(
  SampleID  = sample_ids,
  Treatment = c(rep("Induced", n_induced), rep("Control", n_control)),
  stringsAsFactors = FALSE
)
write.table(sample_info,
  file      = "test-data/test_samples.tabular",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
cat("Written: test-data/test_samples.tabular\n")

# ---- Write DEG lists ----
write.table(
  data.frame(GeneID = gene_ids[module_a[1:40]]),
  file      = "test-data/test_up_degs.tabular",
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)
cat("Written: test-data/test_up_degs.tabular\n")

write.table(
  data.frame(GeneID = gene_ids[module_b[1:40]]),
  file      = "test-data/test_down_degs.tabular",
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)
cat("Written: test-data/test_down_degs.tabular\n")

cat("\nTest data generation complete.\n")
cat("Files are in: test-data/\n")
