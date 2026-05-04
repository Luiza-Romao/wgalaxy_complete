#!/usr/bin/env Rscript
# =============================================================================
# create_test_subset.R
#
# Generates a small but biologically representative test dataset for the
# WGCNA Galaxy Tool from the full sugarcane VST matrix (SRP309574).
#
# Source data:
#   NCBI SRA accession SRP309574
#   Sugarcane (Saccharum officinarum hybrid), 14 samples:
#     7 x Regular Day (RD — 12h light / 12h dark)
#     7 x Short Day   (SD —  8h light / 16h dark)
#   Expression quantified with Salmon; normalized with DESeq2 VST.
#
# Strategy:
#   - 150 genes with highest variance across samples (enriched for genes
#     responding to photoperiod; ensures the pipeline produces modules with
#     detectable module-trait correlations)
#   - 150 genes sampled randomly from the remaining pool (background signal)
#   - Total: 300 genes x 14 samples
#
# This size is sufficient for WGCNA to detect modules while keeping the test
# files small enough for a GitHub repository (< 1 MB).
#
# Usage:
#   Rscript create_test_subset.R \
#     --vst_matrix   path/to/full_vst_matrix.tsv \
#     --sample_info  path/to/sample_info.tsv \
#     --out_dir      test-data
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--vst_matrix",  type = "character",
              help = "Path to the full VST expression matrix (genes x samples, TSV)"),
  make_option("--sample_info", type = "character",
              help = "Path to the sample metadata file (TSV with SampleID and Treatment columns)"),
  make_option("--out_dir",     type = "character", default = "test-data",
              help = "Output directory for test files [default: test-data]"),
  make_option("--n_top",       type = "integer",   default = 150,
              help = "Number of high-variance genes to include [default: 150]"),
  make_option("--n_random",    type = "integer",   default = 150,
              help = "Number of randomly sampled genes to include [default: 150]"),
  make_option("--seed",        type = "integer",   default = 42,
              help = "Random seed for reproducibility [default: 42]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vst_matrix) || is.null(opt$sample_info)) {
  stop("Both --vst_matrix and --sample_info are required.")
}

set.seed(opt$seed)

# ── Load full VST matrix ──────────────────────────────────────────────────────

cat("Loading VST matrix:", opt$vst_matrix, "\n")
vst <- read.table(
  opt$vst_matrix,
  header           = TRUE,
  sep              = "\t",
  row.names        = 1,
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  quote            = ""
)
cat("  Dimensions:", nrow(vst), "genes x", ncol(vst), "samples\n")

# ── Load sample info ──────────────────────────────────────────────────────────

cat("Loading sample info:", opt$sample_info, "\n")
sample_info <- read.table(
  opt$sample_info,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  quote            = ""
)
cat("  Samples:", nrow(sample_info), "\n")

# Align sample info to columns present in the VST matrix
common_samples <- intersect(sample_info[[1]], colnames(vst))
if (length(common_samples) == 0) {
  stop(
    "No sample names match between the VST matrix columns and the sample info. ",
    "Check that the first column of sample_info contains sample IDs matching ",
    "the column headers of the VST matrix."
  )
}
vst         <- vst[, common_samples, drop = FALSE]
sample_info <- sample_info[sample_info[[1]] %in% common_samples, ]
cat("  Common samples used:", length(common_samples), "\n")

# ── Gene selection ────────────────────────────────────────────────────────────

cat("Selecting genes...\n")

gene_var  <- apply(vst, 1, var)
ranked    <- order(gene_var, decreasing = TRUE)

top_genes  <- rownames(vst)[ranked[seq_len(opt$n_top)]]
remaining  <- rownames(vst)[ranked[(opt$n_top + 1):length(ranked)]]
rand_genes <- sample(remaining, min(opt$n_random, length(remaining)))

selected   <- unique(c(top_genes, rand_genes))

cat("  High-variance genes selected:", length(top_genes), "\n")
cat("  Randomly sampled genes added:", length(rand_genes), "\n")
cat("  Total genes in subset:", length(selected), "\n")

vst_subset <- vst[selected, , drop = FALSE]

# ── Write outputs ─────────────────────────────────────────────────────────────

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

vst_out  <- cbind(GeneID = rownames(vst_subset), vst_subset)
vst_path <- file.path(opt$out_dir, "test_vst.tabular")
write.table(vst_out, vst_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Written:", vst_path, "\n")

sample_path <- file.path(opt$out_dir, "test_samples.tabular")
write.table(sample_info, sample_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Written:", sample_path, "\n")

# ── Sanity check ──────────────────────────────────────────────────────────────

cat("\nSanity check:\n")
cat("  VST subset dimensions:", nrow(vst_subset), "x", ncol(vst_subset), "\n")
cat("  Expression range:", round(min(vst_subset), 2), "to",
    round(max(vst_subset), 2), "\n")
cat("  Treatments in metadata:",
    paste(unique(sample_info[[2]]), collapse = ", "), "\n")

med_var <- round(median(apply(vst_subset, 1, var)), 4)
cat("  Median gene variance:", med_var, "\n")

if (med_var < 0.01) {
  warning(
    "Median gene variance is very low (", med_var, "). ",
    "WGCNA may fail to detect modules. Consider increasing --n_top or ",
    "checking that the VST matrix was not already heavily filtered."
  )
}

cat("\nDone. Commit test-data/ to the repository.\n")
