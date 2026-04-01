#!/usr/bin/env Rscript
# =============================================================================
# WGCNA Galaxy Tool — Main Analysis Script
# Weighted Gene Co-expression Network Analysis (WGCNA)
# Compatible with Galaxy bioinformatics platform
#
# Pipeline:
#   1.  Data loading and sample metadata
#   2.  Expression filtering (presence and variance)
#   3.  Trait matrix preparation
#   4.  Soft-thresholding power selection
#   5.  Blockwise network construction and module detection
#   6.  Module eigengene calculation and module-trait correlation
#   7.  Gene significance and module membership (top module)
#   8.  Hub gene identification
#   9.  kME-based gene filtering
#   10. Cytoscape-ready network export (high-correlation modules)
#   11. Module-trait heatmap (selected modules)
#   12. DEG overlap analysis (optional)
# =============================================================================

# ----------------------------------------------------------------------------- 
# Argument parsing
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  # --- Required inputs ---
  make_option("--vst_matrix",
    type = "character", default = NULL,
    help = "VST-normalized expression matrix (TSV/tabular): genes x samples, first column = gene IDs"),
  make_option("--sample_info",
    type = "character", default = NULL,
    help = "Sample metadata file (TSV): columns SampleID and Treatment (or custom --trait_col)"),

  # --- Optional DEG inputs ---
  make_option("--up_degs",
    type = "character", default = NULL,
    help = "Optional: tabular file of up-regulated gene IDs (first column used)"),
  make_option("--down_degs",
    type = "character", default = NULL,
    help = "Optional: tabular file of down-regulated gene IDs (first column used)"),

  # --- Column names ---
  make_option("--sample_col",
    type = "character", default = "SampleID",
    help = "Column name for sample IDs in metadata [default: SampleID]"),
  make_option("--treatment_col",
    type = "character", default = "Treatment",
    help = "Column name for treatment/group in metadata [default: Treatment]"),
  make_option("--case_label",
    type = "character", default = NULL,
    help = "Label of the 'case' group (e.g. Induced). Auto-detected if NULL"),

  # --- Filtering parameters ---
  make_option("--presence_pct",
    type = "double", default = 0.6,
    help = "Min fraction of samples with non-zero expression [default: 0.60]"),
  make_option("--var_percentile",
    type = "integer", default = 40,
    help = "Variance percentile cutoff for gene filtering [default: 40]"),

  # --- WGCNA network parameters ---
  make_option("--soft_power",
    type = "integer", default = 0,
    help = "Soft-thresholding power (0 = auto-detect) [default: 0]"),
  make_option("--network_type",
    type = "character", default = "signed",
    help = "Network type: signed or unsigned [default: signed]"),
  make_option("--min_module_size",
    type = "integer", default = 50,
    help = "Minimum module size [default: 50]"),
  make_option("--merge_cut_height",
    type = "double", default = 0.25,
    help = "Dendrogram cut height for module merging [default: 0.25]"),
  make_option("--max_block_size",
    type = "integer", default = 5000,
    help = "Max block size for blockwiseModules [default: 5000]"),
  make_option("--n_threads",
    type = "integer", default = 4,
    help = "Number of threads for parallel WGCNA [default: 4]"),

  # --- Post-processing thresholds ---
  make_option("--kme_threshold",
    type = "double", default = 0.8,
    help = "Module membership (|kME|) threshold for gene filtering [default: 0.8]"),
  make_option("--cor_threshold",
    type = "double", default = 0.6,
    help = "|r| threshold for selecting modules for Cytoscape export [default: 0.6]"),
  make_option("--tom_percentile",
    type = "double", default = 0.95,
    help = "TOM percentile cutoff for defining network edges [default: 0.95]"),

  # --- Output paths (supplied by Galaxy) ---
  make_option("--out_plot_clustering",   type = "character", default = "plot_clustering.png"),
  make_option("--out_plot_soft",         type = "character", default = "plot_soft_threshold.png"),
  make_option("--out_plot_dendro",       type = "character", default = "plot_dendrogram.png"),
  make_option("--out_plot_heatmap_all",  type = "character", default = "plot_heatmap_all.png"),
  make_option("--out_plot_mm_gs",        type = "character", default = "plot_mm_gs.png"),
  make_option("--out_plot_hubs",         type = "character", default = "plot_hubs.png"),
  make_option("--out_plot_heatmap_sel",  type = "character", default = "plot_heatmap_selected.png"),
  make_option("--out_plot_deg_bar",      type = "character", default = "plot_deg_bar.png"),
  make_option("--out_table_modules",     type = "character", default = "table_module_assignments.tsv"),
  make_option("--out_table_gs_mm",       type = "character", default = "table_gs_mm.tsv"),
  make_option("--out_table_hubs",        type = "character", default = "table_hub_genes.tsv"),
  make_option("--out_table_sel_modules", type = "character", default = "table_selected_modules.tsv"),
  make_option("--out_table_deg_summary", type = "character", default = "table_deg_network_summary.tsv"),
  make_option("--out_cytoscape_dir",     type = "character", default = "cytoscape_networks"),
  make_option("--out_rdata",             type = "character", default = "wgcna_results.RData")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$vst_matrix))  stop("--vst_matrix is required")
if (is.null(opt$sample_info)) stop("--sample_info is required")
if (!file.exists(opt$vst_matrix))  stop("VST matrix file not found: ", opt$vst_matrix)
if (!file.exists(opt$sample_info)) stop("Sample info file not found: ", opt$sample_info)

# ----------------------------------------------------------------------------- 
# Load libraries
# -----------------------------------------------------------------------------
cat("=== Loading libraries ===\n")
suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(ggrepel)
  library(reshape2)
  library(pheatmap)
  library(igraph)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = opt$n_threads)
cat("WGCNA threads:", opt$n_threads, "\n")

# Create Cytoscape output directory
dir.create(opt$out_cytoscape_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Data loading and sample metadata
# =============================================================================
cat("\n=== STEP 1: Loading data ===\n")

sample_info <- tryCatch(
  read.table(opt$sample_info, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop("Failed to read sample info: ", e$message)
)

# Validate required columns
if (!opt$sample_col %in% colnames(sample_info)) {
  stop("Column '", opt$sample_col, "' not found in sample info. Available: ",
       paste(colnames(sample_info), collapse = ", "))
}
if (!opt$treatment_col %in% colnames(sample_info)) {
  stop("Column '", opt$treatment_col, "' not found in sample info. Available: ",
       paste(colnames(sample_info), collapse = ", "))
}

rownames(sample_info) <- sample_info[[opt$sample_col]]

cat("Sample distribution:\n")
print(table(sample_info[[opt$treatment_col]]))

# Load VST matrix
expr_data <- tryCatch(
  read.table(opt$vst_matrix, header = TRUE, row.names = 1,
             sep = "\t", check.names = FALSE),
  error = function(e) stop("Failed to read VST matrix: ", e$message)
)

cat("\nRaw expression matrix dimensions:", dim(expr_data), "\n")

# Align samples
common_samples <- intersect(colnames(expr_data), rownames(sample_info))
if (length(common_samples) == 0) {
  stop("No common sample IDs between expression matrix and sample metadata.")
}
if (length(common_samples) < ncol(expr_data)) {
  warning("Some samples in expression matrix not in metadata; subsetting to common samples.")
}
expr_data   <- expr_data[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

cat("Samples used in analysis:", length(common_samples), "\n")

# =============================================================================
# STEP 2: Expression filtering (presence + variance)
# =============================================================================
cat("\n=== STEP 2: Filtering genes ===\n")

n_samples      <- ncol(expr_data)
presence_min   <- ceiling(opt$presence_pct * n_samples)

cat("Presence filter: genes must be non-zero in >=", presence_min,
    "samples (", opt$presence_pct * 100, "% )\n")

presence_ok <- rowSums(expr_data > 0) >= presence_min
expr_data   <- expr_data[presence_ok, ]
cat("Genes retained after presence filter:", nrow(expr_data), "\n")

# Check for NA / Inf
na_count  <- sum(is.na(expr_data))
inf_count <- sum(is.infinite(as.matrix(expr_data)))
if (na_count > 0)  warning(na_count, " NA values detected; imputing with row means.")
if (inf_count > 0) warning(inf_count, " Inf values detected; replacing with row means.")

if (na_count > 0 || inf_count > 0) {
  expr_mat <- as.matrix(expr_data)
  for (i in seq_len(nrow(expr_mat))) {
    bad <- is.na(expr_mat[i, ]) | is.infinite(expr_mat[i, ])
    if (any(bad)) expr_mat[i, bad] <- mean(expr_mat[i, !bad], na.rm = TRUE)
  }
  expr_data <- as.data.frame(expr_mat)
}

# Variance filter
gene_vars     <- apply(expr_data, 1, var)
var_cutoff    <- quantile(gene_vars, opt$var_percentile / 100)
cat("Variance threshold (", opt$var_percentile, "th percentile):", round(var_cutoff, 4), "\n")
expr_data <- expr_data[gene_vars > var_cutoff, ]
cat("Genes retained after variance filter:", nrow(expr_data), "\n")

# Transpose: WGCNA requires samples x genes
expr_matrix <- t(as.matrix(expr_data))
cat("Final expression matrix:", nrow(expr_matrix), "samples x", ncol(expr_matrix), "genes\n")

# =============================================================================
# STEP 3: Trait matrix preparation
# =============================================================================
cat("\n=== STEP 3: Building trait matrix ===\n")

groups  <- sample_info[[opt$treatment_col]]
u_groups <- unique(groups)

# Auto-detect case label
case_label <- opt$case_label
if (is.null(case_label) || case_label == "") {
  # Prefer alphabetically last if not specified (or user can specify)
  case_label <- sort(u_groups)[length(sort(u_groups))]
  cat("Auto-detected case label:", case_label, "\n")
} else {
  if (!case_label %in% u_groups) {
    stop("case_label '", case_label, "' not found in treatment column. Available: ",
         paste(u_groups, collapse = ", "))
  }
}

trait_data <- data.frame(
  Case = as.integer(groups == case_label),
  row.names = rownames(sample_info)
)
colnames(trait_data) <- case_label
trait_data <- trait_data[rownames(expr_matrix), , drop = FALSE]

cat("Trait matrix (first rows):\n")
print(head(trait_data))

# =============================================================================
# STEP 4: Sample clustering & outlier check (plot)
# =============================================================================
cat("\n=== STEP 4: Sample clustering ===\n")

sample_tree <- hclust(dist(expr_matrix), method = "average")

png(opt$out_plot_clustering, width = 1200, height = 600, res = 120)
par(mar = c(0, 5, 4, 0))
plot(sample_tree,
     main  = "Sample clustering (check for outliers)",
     sub   = "",
     xlab  = "",
     cex.lab  = 1.2,
     cex.axis = 1.0,
     cex.main = 1.4)
# Draw trait as colored bar below the dendrogram
traitColors <- numbers2colors(trait_data, signed = FALSE)
plotDendroAndColors(sample_tree, traitColors,
                    groupLabels = colnames(trait_data),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
cat("Sample clustering plot saved.\n")

# =============================================================================
# STEP 5: Soft-thresholding power selection
# =============================================================================
cat("\n=== STEP 5: Soft-threshold selection ===\n")

powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(expr_matrix,
                         powerVector  = powers,
                         networkType  = opt$network_type,
                         verbose      = 2)

png(opt$out_plot_soft, width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))

# Scale-free topology fit index
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)",
     ylab = expression(paste("Scale-free topology fit index, ", R^2)),
     type = "n",
     main = "Scale-free topology fit")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "steelblue")
abline(h = 0.85, col = "firebrick", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, col = "firebrick")

par(mfrow = c(1, 1))
dev.off()

soft_power <- if (opt$soft_power > 0) opt$soft_power else sft$powerEstimate
if (is.na(soft_power) || is.null(soft_power)) {
  soft_power <- 6
  warning("Auto power estimation failed. Using power = 6. Check soft-threshold plot.")
}
cat("Soft-thresholding power selected:", soft_power, "\n")

# =============================================================================
# STEP 6: Blockwise network construction and module detection
# =============================================================================
cat("\n=== STEP 6: Network construction ===\n")

cor <- WGCNA::cor   # avoid stats::cor collision

net <- blockwiseModules(
  expr_matrix,
  power             = soft_power,
  networkType       = opt$network_type,
  TOMType           = opt$network_type,
  minModuleSize     = opt$min_module_size,
  reassignThreshold = 0,
  mergeCutHeight    = opt$merge_cut_height,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize      = opt$max_block_size,
  saveTOMs          = FALSE,
  verbose           = 3
)

cor <- stats::cor   # restore

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_matrix)
module_colors_all    <- module_colors

n_modules <- length(unique(module_colors)) - 1  # exclude grey
cat("\nNumber of co-expression modules identified:", n_modules, "\n")
print(table(module_colors))

# Dendrogram plot (block 1)
png(opt$out_plot_dendro, width = 1400, height = 800, res = 120)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module color",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene dendrogram and module assignment (block 1)"
)
dev.off()
cat("Dendrogram plot saved.\n")

# Export full module assignment table
module_assignment <- data.frame(
  Gene   = colnames(expr_matrix),
  Module = module_colors,
  stringsAsFactors = FALSE
)
write.table(module_assignment, file = opt$out_table_modules,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Module assignment table saved.\n")

# =============================================================================
# STEP 7: Module eigengenes and module-trait correlation
# =============================================================================
cat("\n=== STEP 7: Module-trait correlations ===\n")

MEs               <- orderMEs(net$MEs)
module_trait_cor  <- cor(MEs, trait_data, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_matrix))

text_matrix <- paste0(
  signif(module_trait_cor,  2), "\n(",
  signif(module_trait_pval, 1), ")"
)
dim(text_matrix) <- dim(module_trait_cor)

# Full heatmap
png(opt$out_plot_heatmap_all,
    width  = max(600,  150 * ncol(module_trait_cor)),
    height = max(800, 14  * nrow(module_trait_cor)),
    res    = 120)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix        = module_trait_cor,
  xLabels       = colnames(trait_data),
  yLabels       = names(MEs),
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = text_matrix,
  setStdMargins = FALSE,
  cex.text      = 0.5,
  zlim          = c(-1, 1),
  main          = "Module eigengene – trait correlation"
)
dev.off()
cat("Full module-trait heatmap saved.\n")

# =============================================================================
# STEP 8: Focus on the top module
# =============================================================================
cat("\n=== STEP 8: Top module analysis ===\n")

trait_name    <- colnames(trait_data)[1]
cor_abs       <- abs(module_trait_cor[, trait_name])
best_ME_name  <- names(MEs)[which.max(cor_abs)]
best_ME_num   <- as.numeric(gsub("ME", "", best_ME_name))
module_color  <- labels2colors(best_ME_num)

genes_in_module <- colnames(expr_matrix)[net$colors == best_ME_num]

cat("\nTop module correlated with '", trait_name, "':", best_ME_name,
    "( color:", module_color, ")\n")
cat("Genes in module:", length(genes_in_module), "\n")
cat("Pearson r with", trait_name, ":",
    round(module_trait_cor[best_ME_name, trait_name], 4), "\n")

# Gene significance (GS) — correlation with trait
gene_sig      <- as.numeric(cor(expr_matrix, trait_data[[trait_name]], use = "p"))
gene_sig_pval <- corPvalueStudent(gene_sig, nrow(expr_matrix))
names(gene_sig)      <- colnames(expr_matrix)
names(gene_sig_pval) <- colnames(expr_matrix)

# Module membership (MM / kME)
MM      <- as.numeric(cor(expr_matrix[, genes_in_module],
                          MEs[, best_ME_name], use = "p"))
MM_pval <- corPvalueStudent(MM, nrow(expr_matrix))

module_df <- data.frame(
  Gene             = genes_in_module,
  ModuleMembership = MM,
  MM_pvalue        = MM_pval,
  GeneSignificance = gene_sig[genes_in_module],
  GS_pvalue        = gene_sig_pval[genes_in_module],
  stringsAsFactors = FALSE
)
module_df <- module_df[order(-abs(module_df$ModuleMembership)), ]

write.table(module_df, file = opt$out_table_gs_mm,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("GS/MM table saved.\n")

# MM vs GS scatter plot
png(opt$out_plot_mm_gs, width = 900, height = 700, res = 120)
col_safe <- tryCatch(adjustcolor(module_color, alpha.f = 0.5),
                     error = function(e) "steelblue")
plot(
  abs(module_df$ModuleMembership),
  abs(module_df$GeneSignificance),
  xlab = paste("Module membership (kME) —", module_color),
  ylab = paste("Gene significance for", trait_name),
  main = paste("MM vs. GS —", module_color, "module"),
  pch  = 20, col = col_safe
)
abline(lm(abs(GeneSignificance) ~ abs(ModuleMembership), data = module_df),
       col = "firebrick", lwd = 2)
cor_val <- cor(abs(module_df$ModuleMembership), abs(module_df$GeneSignificance))
legend("topleft", bty = "n",
       legend = paste("r =", round(cor_val, 3)))
dev.off()
cat("MM vs GS plot saved.\n")

# =============================================================================
# STEP 9: Hub gene identification
# =============================================================================
cat("\n=== STEP 9: Hub genes ===\n")

adj_matrix   <- adjacency(expr_matrix[, genes_in_module],
                           power = soft_power, type = opt$network_type)
connectivity <- colSums(adj_matrix) - 1

hub_info <- data.frame(
  Gene             = genes_in_module,
  Connectivity     = connectivity,
  ModuleMembership = MM,
  GeneSignificance = gene_sig[genes_in_module],
  stringsAsFactors = FALSE
)
hub_info <- hub_info[order(-hub_info$Connectivity), ]

hub_threshold  <- quantile(hub_info$Connectivity, 0.90)
hub_info$isHub <- hub_info$Connectivity > hub_threshold
top_hubs       <- hub_info[hub_info$isHub, ]

cat("Hub genes identified (top 10% connectivity):", nrow(top_hubs), "\n")
cat("Top 10 hubs:\n")
print(head(top_hubs, 10))

write.table(top_hubs, file = opt$out_table_hubs,
            sep = "\t", row.names = FALSE, quote = FALSE)

png(opt$out_plot_hubs, width = 900, height = 700, res = 120)
plot(
  hub_info$ModuleMembership,
  hub_info$Connectivity,
  xlab = "Module membership (kME)",
  ylab = "Intramodular connectivity",
  main = paste("Hub gene identification —", module_color, "module"),
  pch  = 20,
  col  = ifelse(hub_info$isHub, "firebrick", "grey60")
)
legend("topleft",
       legend = c("Hub (top 10%)", "Other genes"),
       col    = c("firebrick", "grey60"),
       pch    = 20, bty = "n")
dev.off()
cat("Hub gene plot saved.\n")

# =============================================================================
# STEP 10: kME-based gene filtering
# =============================================================================
cat("\n=== STEP 10: kME filtering ===\n")

genes_high_kME       <- module_df[abs(module_df$ModuleMembership) > opt$kme_threshold, ]
genes_in_module_final <- genes_high_kME$Gene

cat("kME filter (|kME| >", opt$kme_threshold, "):\n")
cat("  Genes before filter:", nrow(module_df), "\n")
cat("  Genes after  filter:", length(genes_in_module_final), "\n")

# =============================================================================
# STEP 11: Cytoscape-ready network export (modules with |r| >= threshold)
# =============================================================================
cat("\n=== STEP 11: Cytoscape network export ===\n")

# Recalculate correlations for all modules
module_numbers    <- net$colors
module_colors_all <- labels2colors(module_numbers)
names(module_colors_all) <- colnames(expr_matrix)

module_trait_cor <- cor(MEs, trait_data, use = "p")
cor_with_trait   <- module_trait_cor[, trait_name]

unique_modules <- sort(unique(module_numbers))
unique_modules <- unique_modules[unique_modules != 0]

module_summary <- data.frame(
  ModuleNum   = unique_modules,
  ModuleColor = labels2colors(unique_modules),
  Correlation = cor_with_trait[paste0("ME", unique_modules)],
  stringsAsFactors = FALSE
)

modules_selected <- module_summary[
  abs(module_summary$Correlation) >= opt$cor_threshold, ]

cat("Modules with |r| >=", opt$cor_threshold, ":",
    nrow(modules_selected), "of", nrow(module_summary), "\n")

if (nrow(modules_selected) == 0) {
  warning("No modules passed the correlation threshold (", opt$cor_threshold, ").",
          " Lowering to 0.3 for export.")
  modules_selected <- module_summary[
    abs(module_summary$Correlation) >= 0.3, ]
}

modules_selected <- modules_selected[order(-abs(modules_selected$Correlation)), ]
write.table(modules_selected, file = opt$out_table_sel_modules,
            sep = "\t", row.names = FALSE, quote = FALSE)

selected_modules <- modules_selected$ModuleColor

# Network export loop
for (i in seq_len(nrow(modules_selected))) {

  mod_num   <- modules_selected$ModuleNum[i]
  mod_color <- modules_selected$ModuleColor[i]
  mod_cor   <- round(abs(modules_selected$Correlation[i]), 3)

  cat("\n-------------------------------------------\n")
  cat("Processing module:", mod_color, "(ME", mod_num, ") |r| =", mod_cor, "\n")

  genes_mod <- colnames(expr_matrix)[net$colors == mod_num]
  cat("  Genes in module:", length(genes_mod), "\n")

  # Recalculate kME for this module
  me_col <- paste0("ME", mod_num)
  if (!me_col %in% colnames(MEs)) {
    cat("  ME not found in MEs object, skipping.\n")
    next
  }
  kME_mod <- as.numeric(cor(expr_matrix[, genes_mod], MEs[, me_col], use = "p"))
  names(kME_mod) <- genes_mod

  genes_kME <- names(kME_mod)[abs(kME_mod) > opt$kme_threshold]
  cat("  Genes after |kME| >", opt$kme_threshold, "filter:", length(genes_kME), "\n")

  if (length(genes_kME) < 3) {
    cat("  Too few genes after kME filter, skipping network export.\n")
    next
  }

  # Build adjacency and TOM
  adj_mod <- adjacency(expr_matrix[, genes_kME],
                        power = soft_power, type = opt$network_type)
  TOM_mod <- TOMsimilarity(adj_mod, TOMType = opt$network_type, verbose = 0)
  rownames(TOM_mod) <- colnames(TOM_mod) <- genes_kME

  tom_cutoff <- quantile(TOM_mod[upper.tri(TOM_mod)], opt$tom_percentile)
  cat("  TOM threshold (", opt$tom_percentile * 100, "th pct):", round(tom_cutoff, 4), "\n")

  # Edge table
  idx <- which(upper.tri(TOM_mod), arr.ind = TRUE)
  edge_df <- data.frame(
    fromNode = genes_kME[idx[, 1]],
    toNode   = genes_kME[idx[, 2]],
    weight   = TOM_mod[idx],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[edge_df$weight >= tom_cutoff, ]

  # Node table
  gs_mod  <- gene_sig[genes_kME]
  kme_mod <- kME_mod[genes_kME]
  conn_mod <- colSums(adj_mod[genes_kME, genes_kME]) - 1
  node_df <- data.frame(
    node             = genes_kME,
    module           = mod_color,
    kME              = round(kme_mod, 4),
    GeneSignificance = round(gs_mod, 4),
    Connectivity     = round(conn_mod, 4),
    stringsAsFactors = FALSE
  )

  edge_file <- file.path(opt$out_cytoscape_dir,
                          paste0("edges_", mod_color, "_kMEfilt.txt"))
  node_file <- file.path(opt$out_cytoscape_dir,
                          paste0("nodes_", mod_color, "_kMEfilt.txt"))

  write.table(edge_df, file = edge_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(node_df, file = node_file, sep = "\t", row.names = FALSE, quote = FALSE)

  cat("  Edge file:", edge_file, "\n")
  cat("  Node file:", node_file, "\n")
}

cat("\nNetwork export complete for", nrow(modules_selected), "modules.\n")

# =============================================================================
# STEP 12: Module-trait heatmap for selected modules
# =============================================================================
cat("\n=== STEP 12: Selected module heatmap ===\n")

sel_ME_names <- paste0("ME", modules_selected$ModuleNum)
sel_ME_names <- sel_ME_names[sel_ME_names %in% rownames(module_trait_cor)]

if (length(sel_ME_names) > 0) {
  sel_cor  <- module_trait_cor[sel_ME_names, , drop = FALSE]
  sel_pval <- module_trait_pval[sel_ME_names, , drop = FALSE]

  sel_text <- paste0(signif(sel_cor, 2), "\n(", signif(sel_pval, 1), ")")
  dim(sel_text) <- dim(sel_cor)

  sel_labels <- modules_selected$ModuleColor[
    match(sel_ME_names, paste0("ME", modules_selected$ModuleNum))]

  png(opt$out_plot_heatmap_sel,
      width  = max(500, 200 * ncol(sel_cor)),
      height = max(400, 50  * nrow(sel_cor)),
      res    = 120)
  par(mar = c(6, 12, 3, 3))
  labeledHeatmap(
    Matrix        = sel_cor,
    xLabels       = colnames(trait_data),
    yLabels       = sel_labels,
    ySymbols      = sel_ME_names,
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = sel_text,
    setStdMargins = FALSE,
    cex.text      = 0.8,
    zlim          = c(-1, 1),
    main          = paste("Selected modules (|r| >=", opt$cor_threshold, ")")
  )
  dev.off()
  cat("Selected module heatmap saved.\n")
}

# =============================================================================
# STEP 13 (optional): DEG overlap analysis
# =============================================================================
has_degs <- !is.null(opt$up_degs) && !is.null(opt$down_degs) &&
  opt$up_degs != "None" && opt$down_degs != "None" &&
  file.exists(opt$up_degs) && file.exists(opt$down_degs)

if (has_degs) {
  cat("\n=== STEP 13: DEG overlap analysis ===\n")

  up_df    <- read.table(opt$up_degs,   header = FALSE, sep = "\t")
  down_df  <- read.table(opt$down_degs, header = FALSE, sep = "\t")
  up_genes   <- unique(as.character(up_df[, 1]))
  down_genes <- unique(as.character(down_df[, 1]))
  up_genes   <- up_genes[!is.na(up_genes) & up_genes != ""]
  down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]

  cat("Up-regulated genes:", length(up_genes), "\n")
  cat("Down-regulated genes:", length(down_genes), "\n")

  degs <- data.frame(
    Gene   = c(up_genes, down_genes),
    Status = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
    stringsAsFactors = FALSE
  )
  degs <- degs[!duplicated(degs$Gene), ]
  cat("Unique DEGs:", nrow(degs), "\n")

  node_dir <- opt$out_cytoscape_dir
  network_summary <- do.call(rbind, lapply(selected_modules, function(mod) {
    node_file <- file.path(node_dir, paste0("nodes_", mod, "_kMEfilt.txt"))
    if (!file.exists(node_file)) {
      return(NULL)
    }
    nodes <- read.table(node_file, header = TRUE, sep = "\t")
    gin   <- nodes$node
    data.frame(
      Module       = mod,
      NetworkNodes = length(gin),
      DEGs_total   = length(intersect(gin, degs$Gene)),
      DEGs_up      = length(intersect(gin, degs$Gene[degs$Status == "Up"])),
      DEGs_down    = length(intersect(gin, degs$Gene[degs$Status == "Down"])),
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(network_summary) && nrow(network_summary) > 0) {
    network_summary$Proportion_DEGs <-
      network_summary$DEGs_total / network_summary$NetworkNodes

    write.table(network_summary, file = opt$out_table_deg_summary,
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Update node files with DEG status
    for (mod in selected_modules) {
      node_file <- file.path(node_dir, paste0("nodes_", mod, "_kMEfilt.txt"))
      if (!file.exists(node_file)) next
      nodes <- read.table(node_file, header = TRUE, sep = "\t")
      nodes$DE_status <- "Not_DE"
      nodes$DE_status[nodes$node %in% degs$Gene[degs$Status == "Up"]]   <- "Up"
      nodes$DE_status[nodes$node %in% degs$Gene[degs$Status == "Down"]] <- "Down"
      write.table(nodes,
        file.path(node_dir, paste0("nodes_", mod, "_kMEfilt_with_DEG.txt")),
        sep = "\t", row.names = FALSE, quote = FALSE)
    }

    # Bar plot
    plot_data <- melt(network_summary[, c("Module", "DEGs_up", "DEGs_down")],
                      id.vars      = "Module",
                      variable.name = "Direction",
                      value.name   = "Count")
    plot_data$Direction <- factor(plot_data$Direction,
      levels = c("DEGs_up", "DEGs_down"),
      labels = c("Up-regulated", "Down-regulated"))

    p_deg <- ggplot(plot_data, aes(x = Module, y = Count, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Up-regulated" = "firebrick",
                                    "Down-regulated" = "steelblue")) +
      labs(title = "DEGs per co-expression module (network nodes)",
           x = "Module", y = "Number of DEGs") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(opt$out_plot_deg_bar, plot = p_deg,
           width = max(8, nrow(network_summary) * 0.8), height = 5, dpi = 150)
    cat("DEG barplot saved.\n")
    cat("DEG network summary:\n")
    print(network_summary)
  } else {
    cat("No network node files found for DEG overlap. Skipping.\n")
  }
} else {
  cat("\nNo DEG files provided; skipping overlap analysis.\n")
  # Write empty placeholder for output
  write.table(data.frame(Note = "No DEG files provided"),
              file = opt$out_table_deg_summary, sep = "\t",
              row.names = FALSE, quote = FALSE)
  # Empty DEG plot placeholder
  png(opt$out_plot_deg_bar, width = 400, height = 200, res = 72)
  plot.new()
  text(0.5, 0.5, "No DEG data provided", cex = 1.2)
  dev.off()
}

# =============================================================================
# Save RData for downstream use
# =============================================================================
cat("\n=== Saving RData workspace ===\n")
save(net, expr_matrix, trait_data, MEs, module_colors, module_colors_all,
     module_trait_cor, module_trait_pval, gene_sig, module_df,
     modules_selected, selected_modules, soft_power,
     file = opt$out_rdata)
cat("RData saved to:", opt$out_rdata, "\n")

cat("\n=== WGCNA analysis complete! ===\n")
sessionInfo()
