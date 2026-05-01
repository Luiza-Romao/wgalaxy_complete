#!/usr/bin/env Rscript
# WGCNA Galaxy Tool

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--vst_matrix",   type="character", default=NULL),
  make_option("--sample_info",  type="character", default=NULL),
  make_option("--up_degs",      type="character", default="None"),
  make_option("--down_degs",    type="character", default="None"),
  make_option("--sample_col",   type="character", default="SampleID"),
  make_option("--treatment_col",type="character", default="Treatment"),
  make_option("--case_label",   type="character", default=""),
  make_option("--presence_pct",   type="double",  default=0.6),
  make_option("--var_percentile", type="integer", default=40),
  make_option("--soft_power",       type="integer",   default=0),
  make_option("--network_type",     type="character", default="signed"),
  make_option("--min_module_size",  type="integer",   default=50),
  make_option("--merge_cut_height", type="double",    default=0.25),
  make_option("--max_block_size",   type="integer",   default=5000),
  make_option("--n_threads",        type="integer",   default=4),
  make_option("--kme_threshold",  type="double", default=0.8),
  make_option("--cor_threshold",  type="double", default=0.6),
  make_option("--tom_percentile", type="double", default=0.95),
  make_option("--run_preservation",     type="character", default="no"),
  make_option("--pres_query_vst_matrix",type="character", default=NULL),
  make_option("--ref_vst_matrix",       type="character", default="None"),
  make_option("--ref_sample_info",      type="character", default="None"),
  make_option("--ref_sample_col",       type="character", default="SampleID"),
  make_option("--ref_treatment_col",    type="character", default="Treatment"),
  make_option("--pres_n_perms",         type="integer",   default=200),
  make_option("--pres_max_genes",       type="integer",   default=2000),
  make_option("--pres_random_seed",     type="integer",   default=12345),
  make_option("--out_plot_clustering",   type="character", default="plot_clustering.png"),
  make_option("--out_plot_soft",         type="character", default="plot_soft_threshold.png"),
  make_option("--out_plot_dendro",       type="character", default="plot_dendrogram.png"),
  make_option("--out_plot_heatmap_all",  type="character", default="plot_heatmap_all.png"),
  make_option("--out_plot_mm_gs",        type="character", default="plot_mm_gs.png"),
  make_option("--out_plot_hubs",         type="character", default="plot_hubs.png"),
  make_option("--out_plot_heatmap_sel",  type="character", default="plot_heatmap_selected.png"),
  make_option("--out_plot_deg_bar",      type="character", default="plot_deg_bar.png"),
  make_option("--out_table_modules",     type="character", default="table_module_assignments.tsv"),
  make_option("--out_table_gs_mm",       type="character", default="table_gs_mm.tsv"),
  make_option("--out_table_hubs",        type="character", default="table_hub_genes.tsv"),
  make_option("--out_table_sel_modules", type="character", default="table_selected_modules.tsv"),
  make_option("--out_table_deg_summary", type="character", default="table_deg_summary.tsv"),
  make_option("--out_cytoscape_dir",     type="character", default="cytoscape_output"),
  make_option("--out_rdata",             type="character", default="wgcna_results.RData"),
  make_option("--out_log",               type="character", default="analysis_log.txt"),
  make_option("--out_table_preservation",type="character", default="table_module_preservation.tsv"),
  make_option("--out_plot_zsummary",     type="character", default="plot_zsummary.png"),
  make_option("--out_plot_pres_heatmap", type="character", default="plot_preservation_heatmap.png"),
  make_option("--out_table_pres_class",  type="character", default="table_preservation_classification.tsv")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Logging ------------------------------------------------------------
log_con <- file(opt$out_log, open="wt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
                paste(..., collapse=" "))
  writeLines(msg, log_con); flush(log_con)
  cat(msg, "\n")
}

# ---- Helpers ------------------------------------------------------------
placeholder_png <- function(path, msg="Analysis not requested") {
  png(path, width=500, height=200, res=72)
  plot.new(); text(0.5, 0.5, msg, cex=1.2, col="grey40"); dev.off()
}

placeholder_tsv <- function(path, msg="Analysis not requested") {
  write.table(data.frame(Note=msg), file=path,
              sep="\t", row.names=FALSE, quote=FALSE)
}

load_tabular <- function(path, row1=FALSE, lbl="file") {
  read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE,
             check.names=FALSE, row.names=if(row1) 1 else NULL)
}

suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2); library(ggrepel)
  library(reshape2); library(pheatmap); library(igraph); library(RColorBrewer)
})
options(stringsAsFactors=FALSE)
enableWGCNAThreads(nThreads=opt$n_threads)
dir.create(opt$out_cytoscape_dir, showWarnings=FALSE, recursive=TRUE)

log_msg("WGCNA pipeline started")
log_msg("Threads:", opt$n_threads, "| Network type:", opt$network_type)

# --- STEP 1: Load data ---------------------------------------------------
log_msg("STEP 1: Loading data")
sample_info <- load_tabular(opt$sample_info, lbl="sample info")
rownames(sample_info) <- sample_info[[opt$sample_col]]
expr_data <- load_tabular(opt$vst_matrix, row1=TRUE, lbl="VST matrix")
common <- intersect(colnames(expr_data), rownames(sample_info))
expr_data   <- expr_data[, common, drop=FALSE]
sample_info <- sample_info[common, , drop=FALSE]
log_msg("Samples retained:", length(common), "| Genes loaded:", nrow(expr_data))

# --- STEP 2: Filtering ---------------------------------------------------
log_msg("STEP 2: Gene filtering")
pmin <- ceiling(opt$presence_pct * ncol(expr_data))
expr_data <- expr_data[rowSums(expr_data > 0) >= pmin, ]
gv <- apply(expr_data, 1, var)
vc <- quantile(gv, opt$var_percentile/100)
expr_data <- expr_data[gv > vc, ]
expr_matrix <- t(as.matrix(expr_data))
log_msg("Genes after filtering:", ncol(expr_matrix))

# --- STEP 3: Trait Matrix ------------------------------------------------
log_msg("STEP 3: Building trait matrix")
groups <- sample_info[[opt$treatment_col]]
u_groups <- sort(unique(groups))
trait_data <- as.data.frame(sapply(u_groups, function(g) as.integer(groups == g)))
rownames(trait_data) <- rownames(sample_info)
trait_name <- if (opt$case_label != "" && opt$case_label %in% u_groups) {
  opt$case_label
} else {
  u_groups[length(u_groups)]
}
log_msg("Case (trait of interest):", trait_name)

# --- STEP 4: Sample clustering ------------------------------------------
log_msg("STEP 4: Sample clustering (outlier detection)")
sample_tree <- hclust(dist(expr_matrix), method="average")
png(opt$out_plot_clustering, width=1000, height=600, res=120)
par(mar=c(2, 5, 3, 1))
plot(sample_tree, main="Sample clustering (outlier detection)",
     sub="", xlab="", cex.lab=1.1, cex.axis=1.0, cex.main=1.3)
dev.off()

# --- STEP 5: Soft threshold ---------------------------------------------
log_msg("STEP 5: Soft-thresholding power selection")
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(expr_matrix, powerVector=powers,
                         networkType=opt$network_type, verbose=0)
soft_power <- if (opt$soft_power > 0) {
  opt$soft_power
} else if (!is.na(sft$powerEstimate)) {
  sft$powerEstimate
} else {
  6
}
log_msg("Selected soft power:", soft_power)

png(opt$out_plot_soft, width=1100, height=550, res=120)
par(mfrow=c(1, 2), mar=c(5, 5, 4, 2))
sft_fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
plot(sft$fitIndices[, 1], sft_fit,
     xlab="Soft threshold (power)", ylab="Scale-free topology R^2",
     type="n", main="Scale-free topology fit")
text(sft$fitIndices[, 1], sft_fit, labels=powers, col="red")
abline(h=0.85, col="red", lty=2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft threshold (power)", ylab="Mean connectivity",
     type="n", main="Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, col="red")
dev.off()

# --- STEP 6: Network construction ---------------------------------------
log_msg("STEP 6: Blockwise network construction (power =", soft_power, ")")
net <- blockwiseModules(expr_matrix, power=soft_power,
                        networkType=opt$network_type,
                        TOMType=opt$network_type,
                        minModuleSize=opt$min_module_size,
                        mergeCutHeight=opt$merge_cut_height,
                        numericLabels=FALSE,
                        maxBlockSize=opt$max_block_size,
                        saveTOMs=FALSE, verbose=0)

module_colors <- net$colors
MEs <- orderMEs(net$MEs)
log_msg("Modules detected:", length(unique(module_colors)))

write.table(data.frame(Gene=colnames(expr_matrix), Module=module_colors),
            file=opt$out_table_modules, sep="\t",
            row.names=FALSE, quote=FALSE)

# --- STEP 6.5: Gene dendrogram ------------------------------------------
log_msg("STEP 6.5: Gene dendrogram")
png(opt$out_plot_dendro, width=1100, height=600, res=120)
plotDendroAndColors(net$dendrograms[[1]],
                    module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05,
                    main="Gene clustering dendrogram and module assignment")
dev.off()

# --- STEP 7 & 12: Module-Trait correlation & Heatmaps -------------------
log_msg("STEP 7: Module-trait correlations")
module_trait_cor  <- cor(MEs, trait_data, use="p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_matrix))

# Plot All Heatmap
png(opt$out_plot_heatmap_all, width=1000, height=800, res=120)
labeledHeatmap(Matrix=module_trait_cor, xLabels=colnames(trait_data),
               yLabels=names(MEs),
               colors=blueWhiteRed(50),
               textMatrix=round(module_trait_cor, 2),
               setStdMargins=FALSE, cex.text=0.6,
               main="Module-Trait Correlation")
dev.off()

# Selected Modules Heatmap (|r| >= threshold, excluding grey catch-all module)
mod_cor_case <- module_trait_cor[, trait_name]
keep_mods <- names(mod_cor_case)[abs(mod_cor_case) >= opt$cor_threshold &
                                   names(mod_cor_case) != "MEgrey"]

if (length(keep_mods) > 0) {
  sel_cor  <- module_trait_cor[keep_mods, , drop=FALSE]
  sel_pval <- module_trait_pval[keep_mods, , drop=FALSE]

  clean_mod_names    <- gsub("^ME", "", keep_mods)
  rownames(sel_cor)  <- clean_mod_names
  rownames(sel_pval) <- clean_mod_names

  write.table(
    data.frame(
      Module      = clean_mod_names,
      Correlation = mod_cor_case[keep_mods],
      P_value     = module_trait_pval[keep_mods, trait_name]
    ),
    file=opt$out_table_sel_modules, sep="\t",
    row.names=FALSE, quote=FALSE
  )

  png(opt$out_plot_heatmap_sel, width=800, height=600, res=120)
  labeledHeatmap(Matrix=sel_cor, xLabels=colnames(trait_data),
                 yLabels=rownames(sel_cor),
                 colors=blueWhiteRed(50),
                 textMatrix=round(sel_cor, 2),
                 setStdMargins=FALSE, cex.text=0.6,
                 main=paste("Selected Modules (|r| >=", opt$cor_threshold, ")"))
  dev.off()
} else {
  placeholder_png(opt$out_plot_heatmap_sel,
                  "No modules met correlation threshold")
  placeholder_tsv(opt$out_table_sel_modules,
                  "No modules met correlation threshold")
}

# --- STEP 8: GS + MM table and top-module scatter -----------------------
log_msg("STEP 8: Gene significance and module membership")
gene_module_membership <- as.data.frame(cor(expr_matrix, MEs, use="p"))
mm_pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_module_membership),
                                            nrow(expr_matrix)))
names(gene_module_membership) <- paste0("MM.", gsub("^ME", "", names(MEs)))
names(mm_pvalue)              <- paste0("p.MM.", gsub("^ME", "", names(MEs)))

gene_trait_significance <- as.data.frame(
  cor(expr_matrix, trait_data[, trait_name, drop=FALSE], use="p"))
gs_pvalue <- as.data.frame(
  corPvalueStudent(as.matrix(gene_trait_significance), nrow(expr_matrix)))
names(gene_trait_significance) <- paste0("GS.", trait_name)
names(gs_pvalue)               <- paste0("p.GS.", trait_name)

# Per-gene self MM (membership in its own module)
mm_self <- vapply(seq_along(module_colors), function(i) {
  col_idx <- match(paste0("MM.", module_colors[i]), names(gene_module_membership))
  if (is.na(col_idx)) NA_real_ else gene_module_membership[i, col_idx]
}, numeric(1))
mm_p_self <- vapply(seq_along(module_colors), function(i) {
  col_idx <- match(paste0("p.MM.", module_colors[i]), names(mm_pvalue))
  if (is.na(col_idx)) NA_real_ else mm_pvalue[i, col_idx]
}, numeric(1))

gs_mm_table <- data.frame(
  Gene       = colnames(expr_matrix),
  Module     = module_colors,
  GS         = gene_trait_significance[, 1],
  GS_pvalue  = gs_pvalue[, 1],
  MM         = mm_self,
  MM_pvalue  = mm_p_self
)
write.table(gs_mm_table, file=opt$out_table_gs_mm,
            sep="\t", row.names=FALSE, quote=FALSE)

# Identify top module (highest |cor| with case trait, excluding grey)
non_grey <- mod_cor_case[names(mod_cor_case) != "MEgrey"]
if (length(non_grey) > 0) {
  top_module <- gsub("^ME", "", names(which.max(abs(non_grey))))
  log_msg("Top module:", top_module)
  in_top <- module_colors == top_module

  png(opt$out_plot_mm_gs, width=800, height=700, res=120)
  par(mar=c(5, 5, 4, 2))
  mm_col <- paste0("MM.", top_module)
  scatter_col <- if (top_module %in% colors()) top_module else "darkblue"
  verboseScatterplot(abs(gene_module_membership[in_top, mm_col]),
                     abs(gene_trait_significance[in_top, 1]),
                     xlab=paste("Module Membership (kME) in", top_module, "module"),
                     ylab=paste("Gene Significance for", trait_name),
                     main=paste("MM vs GS \u2014", top_module, "module"),
                     cex.main=1.2, cex.lab=1.1, cex.axis=1.0,
                     col=scatter_col, pch=19)
  dev.off()
} else {
  top_module <- NA_character_
  placeholder_png(opt$out_plot_mm_gs,
                  "No non-grey module passed correlation threshold")
}

# --- STEP 9: Hub gene identification ------------------------------------
log_msg("STEP 9: Hub gene identification")
if (!is.na(top_module)) {
  in_top <- module_colors == top_module
  top_genes <- colnames(expr_matrix)[in_top]

  if (length(top_genes) >= 5) {
    adj_top <- adjacency(expr_matrix[, in_top],
                         power=soft_power, type=opt$network_type)
    kIM <- rowSums(adj_top) - 1     # subtract self-connection
    mm_col <- paste0("MM.", top_module)
    mm_top <- gene_module_membership[in_top, mm_col]
    gs_top <- gene_trait_significance[in_top, 1]

    hub_table <- data.frame(
      Gene   = top_genes,
      Module = top_module,
      kIM    = kIM,
      MM     = mm_top,
      GS     = gs_top
    )
    hub_table <- hub_table[order(-hub_table$kIM), ]
    n_hubs <- max(5, ceiling(0.10 * nrow(hub_table)))
    hub_table$Is_hub <- seq_len(nrow(hub_table)) <= n_hubs
    write.table(hub_table, file=opt$out_table_hubs,
                sep="\t", row.names=FALSE, quote=FALSE)

    hub_col <- if (top_module %in% colors()) top_module else "darkblue"
    png(opt$out_plot_hubs, width=900, height=700, res=120)
    par(mar=c(5, 5, 4, 2))
    plot(hub_table$MM, hub_table$kIM,
         xlab="Module Membership (kME)",
         ylab="Intramodular connectivity (kIM)",
         main=paste("Hub gene identification \u2014", top_module, "module"),
         pch=21,
         bg=ifelse(hub_table$Is_hub, hub_col, "grey80"),
         cex=ifelse(hub_table$Is_hub, 1.6, 1.0),
         cex.lab=1.1, cex.main=1.2)
    if (any(hub_table$Is_hub)) {
      hubs <- hub_table[hub_table$Is_hub, ]
      lbl_n <- min(15, nrow(hubs))
      text(hubs$MM[seq_len(lbl_n)], hubs$kIM[seq_len(lbl_n)],
           labels=hubs$Gene[seq_len(lbl_n)],
           pos=4, cex=0.7, col="black")
    }
    legend("topleft", legend=c("Hub", "Non-hub"),
           pch=21, pt.bg=c(hub_col, "grey80"), bty="n")
    dev.off()
  } else {
    placeholder_png(opt$out_plot_hubs,
                    "Top module too small for hub analysis")
    placeholder_tsv(opt$out_table_hubs,
                    "Top module too small for hub analysis")
  }
} else {
  placeholder_png(opt$out_plot_hubs, "No top module identified")
  placeholder_tsv(opt$out_table_hubs, "No top module identified")
}

# --- STEP 10 & 11: kME filtering + Cytoscape export ---------------------
log_msg("STEP 11: Cytoscape network export")
if (length(keep_mods) > 0) {
  sel_modules_clean <- gsub("^ME", "", keep_mods)
  log_msg("Exporting Cytoscape networks for", length(sel_modules_clean),
          "module(s):", paste(sel_modules_clean, collapse=", "))

  exported <- 0
  for (mod in sel_modules_clean) {
    in_mod <- module_colors == mod
    if (sum(in_mod) < 5) {
      log_msg("  skipping", mod, "(fewer than 5 genes)")
      next
    }

    mm_col <- paste0("MM.", mod)
    if (!(mm_col %in% names(gene_module_membership))) {
      log_msg("  skipping", mod, "(no MM column)")
      next
    }

    gene_mm    <- gene_module_membership[in_mod, mm_col]
    keep_genes <- abs(gene_mm) >= opt$kme_threshold
    if (sum(keep_genes) < 5) {
      # fallback: top 50 by |kME| if threshold is too strict for this module
      ord <- order(-abs(gene_mm))
      keep_genes <- rep(FALSE, length(gene_mm))
      keep_genes[ord[seq_len(min(50, length(ord)))]] <- TRUE
      log_msg("  ", mod, ": kME threshold too strict, using top",
              sum(keep_genes), "genes by |kME|")
    }

    mod_genes_idx <- which(in_mod)[keep_genes]
    mod_genes     <- colnames(expr_matrix)[mod_genes_idx]
    if (length(mod_genes) < 5) next

    adj_mod <- adjacency(expr_matrix[, mod_genes],
                         power=soft_power, type=opt$network_type)
    TOM_mod <- TOMsimilarity(adj_mod, TOMType=opt$network_type, verbose=0)
    dimnames(TOM_mod) <- list(mod_genes, mod_genes)

    # TOM percentile threshold (off-diagonal only)
    tom_off <- TOM_mod
    diag(tom_off) <- NA
    tom_cut <- quantile(tom_off, opt$tom_percentile, na.rm=TRUE)

    edge_file <- file.path(opt$out_cytoscape_dir,
                           paste0("edges_", mod, ".txt"))
    node_file <- file.path(opt$out_cytoscape_dir,
                           paste0("nodes_", mod, ".txt"))

    exportNetworkToCytoscape(TOM_mod,
                             edgeFile=edge_file, nodeFile=node_file,
                             weighted=TRUE, threshold=tom_cut,
                             nodeNames=mod_genes,
                             nodeAttr=rep(mod, length(mod_genes)))
    exported <- exported + 1
  }
  log_msg("Cytoscape networks written:", exported)
  if (exported == 0) {
    write.table(data.frame(Note="No module produced an exportable network"),
                file=file.path(opt$out_cytoscape_dir, "no_networks_exported.txt"),
                sep="\t", row.names=FALSE, quote=FALSE)
  }
} else {
  log_msg("No modules selected \u2014 skipping Cytoscape export")
  write.table(data.frame(Note="No modules met correlation threshold"),
              file=file.path(opt$out_cytoscape_dir,
                             "no_modules_selected.txt"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

# --- STEP 13: DEG overlap (optional) ------------------------------------
log_msg("STEP 13: DEG overlap analysis")
has_degs <- (opt$up_degs   != "None" && file.exists(opt$up_degs)) ||
            (opt$down_degs != "None" && file.exists(opt$down_degs))

if (has_degs) {
  read_deg <- function(p) {
    if (p == "None" || !file.exists(p)) return(character(0))
    x <- read.table(p, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                    fill=TRUE, comment.char="")
    # Drop a header row if it looks like a non-gene label
    if (nrow(x) > 0 && tolower(x[1, 1]) %in%
        c("gene", "geneid", "id", "gene_id", "name", "symbol")) {
      x <- x[-1, , drop=FALSE]
    }
    unique(as.character(x[, 1]))
  }
  up <- read_deg(opt$up_degs)
  dn <- read_deg(opt$down_degs)
  log_msg("Up-DEGs:", length(up), "| Down-DEGs:", length(dn))

  modules_unique <- unique(module_colors)
  deg_summary <- data.frame(
    Module      = modules_unique,
    Module_size = vapply(modules_unique,
                         function(m) sum(module_colors == m), integer(1)),
    Up_DEGs     = vapply(modules_unique,
                         function(m) sum(colnames(expr_matrix)[module_colors == m] %in% up),
                         integer(1)),
    Down_DEGs   = vapply(modules_unique,
                         function(m) sum(colnames(expr_matrix)[module_colors == m] %in% dn),
                         integer(1))
  )
  deg_summary$Total_DEGs <- deg_summary$Up_DEGs + deg_summary$Down_DEGs
  deg_summary$Pct_DEG <- round(100 * deg_summary$Total_DEGs / deg_summary$Module_size, 2)
  deg_summary <- deg_summary[order(-deg_summary$Total_DEGs), ]
  write.table(deg_summary, file=opt$out_table_deg_summary,
              sep="\t", row.names=FALSE, quote=FALSE)

  m_long <- reshape2::melt(
    deg_summary[, c("Module", "Up_DEGs", "Down_DEGs")],
    id.vars="Module", variable.name="Direction", value.name="Count")
  m_long$Module <- factor(m_long$Module, levels=deg_summary$Module)

  png(opt$out_plot_deg_bar, width=1000, height=700, res=120)
  p <- ggplot(m_long, aes(x=Module, y=Count, fill=Direction)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=c(Up_DEGs="#d73027", Down_DEGs="#4575b4")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title="DEG distribution per module",
         x="Module", y="DEG count")
  print(p)
  dev.off()
} else {
  placeholder_png(opt$out_plot_deg_bar,    "DEG analysis not requested")
  placeholder_tsv(opt$out_table_deg_summary, "DEG analysis not requested")
}

# --- STEP 14: Preservation Analysis -------------------------------------
log_msg("STEP 14: Module preservation")
if (opt$run_preservation == "yes") {
  ref_expr     <- load_tabular(opt$ref_vst_matrix, row1=TRUE)
  shared_genes <- intersect(colnames(expr_matrix), rownames(ref_expr))
  log_msg("Shared genes for preservation:", length(shared_genes))

  if (length(shared_genes) > 50) {
    multiData  <- list(Query = list(data = expr_matrix[, shared_genes]),
                       Ref   = list(data = t(ref_expr[shared_genes, ])))
    multiColor <- list(Query = module_colors[shared_genes])

    mp <- modulePreservation(multiData, multiColor, referenceNetworks = 1,
                             nPermutations = opt$pres_n_perms,
                             randomSeed    = opt$pres_random_seed,
                             verbose       = 3)

    stats <- mp$preservation$Z[[1]][[2]]
    stats <- stats[rownames(stats) != "gold", , drop=FALSE]
    write.table(stats, file=opt$out_table_preservation,
                sep="\t", quote=FALSE)

    # Z-summary scatter plot
    png(opt$out_plot_zsummary, width=800, height=600, res=120)
    bg_colors <- ifelse(rownames(stats) %in% colors(), rownames(stats), "grey80")
    plot(stats$moduleSize, stats$Zsummary.pres, pch=21, bg=bg_colors, cex=2,
         main="Module Preservation Z-summary",
         xlab="Module Size", ylab="Zsummary")
    abline(h=c(2, 10), col=c("blue", "red"), lty=2)
    legend("topright",
           legend=c("Z > 10: highly preserved",
                    "Z > 2: moderately preserved"),
           lty=2, col=c("red", "blue"), bty="n", cex=0.8)
    dev.off()

    # Preservation statistics heatmap (all Z-columns)
    z_cols <- grep("^Z", colnames(stats), value=TRUE)
    if (length(z_cols) > 1) {
      z_mat <- as.matrix(stats[, z_cols, drop=FALSE])
      z_mat[z_mat > 30] <- 30                # cap extreme values for display
      png(opt$out_plot_pres_heatmap, width=900, height=700, res=120)
      pheatmap(z_mat, main="Module Preservation Z-statistics",
               cluster_rows=FALSE, cluster_cols=FALSE,
               color=colorRampPalette(c("white", "yellow", "red"))(50),
               border_color=NA, fontsize=9)
      dev.off()
    } else {
      placeholder_png(opt$out_plot_pres_heatmap,
                      "Insufficient Z-statistics for heatmap")
    }

    # Module preservation classification table
    pres_class <- data.frame(
      Module         = rownames(stats),
      Module_size    = stats$moduleSize,
      Zsummary       = stats$Zsummary.pres,
      Classification = ifelse(stats$Zsummary.pres >= 10, "Highly_preserved",
                       ifelse(stats$Zsummary.pres >=  2, "Moderately_preserved",
                                                          "Not_preserved"))
    )
    write.table(pres_class, file=opt$out_table_pres_class,
                sep="\t", row.names=FALSE, quote=FALSE)

  } else {
    placeholder_tsv(opt$out_table_preservation,
                    "Too few shared genes for preservation analysis")
    placeholder_png(opt$out_plot_zsummary,
                    "Too few shared genes for preservation analysis")
    placeholder_png(opt$out_plot_pres_heatmap,
                    "Too few shared genes for preservation analysis")
    placeholder_tsv(opt$out_table_pres_class,
                    "Too few shared genes for preservation analysis")
  }
} else {
  placeholder_tsv(opt$out_table_preservation,
                  "Preservation analysis not requested")
  placeholder_png(opt$out_plot_zsummary,
                  "Preservation analysis not requested")
  placeholder_png(opt$out_plot_pres_heatmap,
                  "Preservation analysis not requested")
  placeholder_tsv(opt$out_table_pres_class,
                  "Preservation analysis not requested")
}

log_msg("Pipeline finished successfully")
save.image(file=opt$out_rdata)
close(log_con)
