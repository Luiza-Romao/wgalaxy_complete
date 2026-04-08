#!/usr/bin/env Rscript
# =============================================================================
# WGCNA Galaxy Tool v2 — Main Analysis Script
# Weighted Gene Co-expression Network Analysis (WGCNA)
# Compatible with Galaxy bioinformatics platform
#
# Pipeline:
#   1.  Data loading and sample metadata
#   2.  Expression filtering (presence and variance)
#   3.  Trait matrix preparation
#   4.  Sample clustering / outlier detection
#   5.  Soft-thresholding power selection
#   6.  Blockwise network construction and module detection
#   7.  Module eigengenes and module-trait correlation
#   8.  Gene significance and module membership (top module)
#   9.  Hub gene identification
#   10. kME-based gene filtering
#   11. Cytoscape-ready network export (high-correlation modules)
#   12. Module-trait heatmap (selected modules)
#   13. DEG overlap analysis              [OPTIONAL]
#   14. Gene Ontology / Enrichment        [OPTIONAL — NEW]
#   15. Module Preservation / Z-summary   [OPTIONAL — NEW]
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  # Required inputs
  make_option("--vst_matrix",   type="character", default=NULL),
  make_option("--sample_info",  type="character", default=NULL),
  # Optional DEG inputs
  make_option("--up_degs",      type="character", default="None"),
  make_option("--down_degs",    type="character", default="None"),
  # Column names
  make_option("--sample_col",   type="character", default="SampleID"),
  make_option("--treatment_col",type="character", default="Treatment"),
  make_option("--case_label",   type="character", default=""),
  # Filtering
  make_option("--presence_pct",   type="double",  default=0.6),
  make_option("--var_percentile", type="integer", default=40),
  # Network
  make_option("--soft_power",       type="integer",   default=0),
  make_option("--network_type",     type="character", default="signed"),
  make_option("--min_module_size",  type="integer",   default=50),
  make_option("--merge_cut_height", type="double",    default=0.25),
  make_option("--max_block_size",   type="integer",   default=5000),
  make_option("--n_threads",        type="integer",   default=4),
  # Thresholds
  make_option("--kme_threshold",  type="double", default=0.8),
  make_option("--cor_threshold",  type="double", default=0.6),
  make_option("--tom_percentile", type="double", default=0.95),
  # ── NEW: GO / Enrichment ──────────────────────────────────────────────────
  make_option("--run_enrichment",  type="character", default="no"),
  make_option("--orgdb_package",   type="character", default="none"),
  make_option("--gene_id_type",    type="character", default="SYMBOL"),
  make_option("--run_kegg",        type="character", default="no"),
  make_option("--kegg_organism",   type="character", default=""),
  make_option("--enrich_pval_cut", type="double",    default=0.05),
  make_option("--enrich_qval_cut", type="double",    default=0.2),
  make_option("--enrich_min_gs",   type="integer",   default=10),
  make_option("--enrich_max_gs",   type="integer",   default=500),
  # ── NEW: Module Preservation ─────────────────────────────────────────────
  make_option("--run_preservation",  type="character", default="no"),
  make_option("--ref_vst_matrix",    type="character", default="None"),
  make_option("--ref_sample_info",   type="character", default="None"),
  make_option("--ref_sample_col",    type="character", default="SampleID"),
  make_option("--ref_treatment_col", type="character", default="Treatment"),
  make_option("--pres_n_perms",      type="integer",   default=200),
  make_option("--pres_max_genes",    type="integer",   default=2000),
  make_option("--pres_random_seed",  type="integer",   default=12345),
  # Output paths — core
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
  make_option("--out_cytoscape_dir",     type="character", default="cytoscape_networks"),
  make_option("--out_rdata",             type="character", default="wgcna_results.RData"),
  # Output paths — enrichment
  make_option("--out_enrich_dir",       type="character", default="enrichment_results"),
  make_option("--out_table_go_summary", type="character", default="table_go_summary.tsv"),
  make_option("--out_plot_go_dot",      type="character", default="plot_go_dotplot.png"),
  make_option("--out_plot_go_bar",      type="character", default="plot_go_barplot.png"),
  # Output paths — preservation
  make_option("--out_table_preservation",type="character", default="table_module_preservation.tsv"),
  make_option("--out_plot_zsummary",     type="character", default="plot_zsummary.png"),
  make_option("--out_plot_pres_heatmap", type="character", default="plot_preservation_heatmap.png"),
  make_option("--out_table_pres_class",  type="character", default="table_preservation_classification.tsv")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$vst_matrix) || opt$vst_matrix == "None")
  stop("--vst_matrix is required")
if (is.null(opt$sample_info) || opt$sample_info == "None")
  stop("--sample_info is required")
if (!file.exists(opt$vst_matrix))
  stop("VST matrix not found: ", opt$vst_matrix)
if (!file.exists(opt$sample_info))
  stop("Sample info not found: ", opt$sample_info)

# ----------------------------------------------------------------------------- 
# Helpers
# -----------------------------------------------------------------------------
placeholder_png <- function(path, msg="Analysis not requested") {
  png(path, width=500, height=200, res=72)
  plot.new(); text(0.5, 0.5, msg, cex=1.2, col="grey40")
  dev.off()
}
placeholder_tsv <- function(path, msg="Analysis not requested") {
  write.table(data.frame(Note=msg), file=path, sep="\t",
              row.names=FALSE, quote=FALSE)
}
load_tabular <- function(path, row1=FALSE, lbl="file") {
  tryCatch(
    read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE,
               check.names=FALSE, row.names=if(row1) 1 else NULL),
    error=function(e) stop("Cannot read ", lbl, ": ", e$message)
  )
}

# -----------------------------------------------------------------------------
# Load core libraries
# -----------------------------------------------------------------------------
cat("=== Loading libraries ===\n")
suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2); library(ggrepel)
  library(reshape2); library(pheatmap); library(igraph); library(RColorBrewer)
})
options(stringsAsFactors=FALSE)
enableWGCNAThreads(nThreads=opt$n_threads)
dir.create(opt$out_cytoscape_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(opt$out_enrich_dir,    showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# STEP 1: Load data
# =============================================================================
cat("\n=== STEP 1: Loading data ===\n")
sample_info <- load_tabular(opt$sample_info, lbl="sample info")
for (col in c(opt$sample_col, opt$treatment_col))
  if (!col %in% colnames(sample_info))
    stop("Column '", col, "' not found in sample info. Available: ",
         paste(colnames(sample_info), collapse=", "))
rownames(sample_info) <- sample_info[[opt$sample_col]]
cat("Sample distribution:\n"); print(table(sample_info[[opt$treatment_col]]))

expr_data <- load_tabular(opt$vst_matrix, row1=TRUE, lbl="VST matrix")
cat("Raw dimensions:", dim(expr_data), "\n")
common <- intersect(colnames(expr_data), rownames(sample_info))
if (!length(common)) stop("No common sample IDs between VST matrix and metadata.")
expr_data   <- expr_data[, common, drop=FALSE]
sample_info <- sample_info[common, , drop=FALSE]
cat("Samples used:", length(common), "\n")

# =============================================================================
# STEP 2: Gene filtering
# =============================================================================
cat("\n=== STEP 2: Filtering genes ===\n")
pmin <- ceiling(opt$presence_pct * ncol(expr_data))
expr_data <- expr_data[rowSums(expr_data > 0) >= pmin, ]
cat("After presence filter:", nrow(expr_data), "\n")
m <- as.matrix(expr_data)
bad <- is.na(m) | is.infinite(m)
if (any(bad)) {
  for (i in which(rowSums(bad)>0)) { ok <- !bad[i,]; m[i,!ok] <- mean(m[i,ok]) }
  expr_data <- as.data.frame(m)
}
gv <- apply(expr_data, 1, var)
vc <- quantile(gv, opt$var_percentile/100)
expr_data <- expr_data[gv > vc, ]
cat("After variance filter:", nrow(expr_data), "genes\n")
expr_matrix <- t(as.matrix(expr_data))
cat("Final matrix:", nrow(expr_matrix), "samples x", ncol(expr_matrix), "genes\n")

# =============================================================================
# STEP 3: Trait matrix
# =============================================================================
cat("\n=== STEP 3: Trait matrix ===\n")
groups   <- sample_info[[opt$treatment_col]]
u_groups <- unique(groups)
case_label <- if (nchar(opt$case_label)==0) {
  lbl <- sort(u_groups)[length(sort(u_groups))]
  cat("Auto case label:", lbl, "\n"); lbl
} else {
  if (!opt$case_label %in% u_groups)
    stop("case_label '", opt$case_label, "' not found in Treatment column.")
  opt$case_label
}
trait_data <- setNames(
  data.frame(as.integer(groups==case_label), row.names=rownames(sample_info)),
  case_label
)
trait_data <- trait_data[rownames(expr_matrix),,drop=FALSE]
trait_name <- colnames(trait_data)[1]
cat("Trait:", trait_name, "\n"); print(head(trait_data))

# =============================================================================
# STEP 4: Sample clustering
# =============================================================================
cat("\n=== STEP 4: Sample clustering ===\n")
stree  <- hclust(dist(expr_matrix), method="average")
trcols <- numbers2colors(trait_data, signed=FALSE)
png(opt$out_plot_clustering, width=1200, height=700, res=120)
plotDendroAndColors(stree, trcols, groupLabels=colnames(trait_data),
                    main="Sample dendrogram and trait heatmap")
dev.off()

# =============================================================================
# STEP 5: Soft-threshold
# =============================================================================
cat("\n=== STEP 5: Soft-threshold ===\n")
powers <- c(1:10, seq(12,30,2))
sft <- pickSoftThreshold(expr_matrix, powerVector=powers,
                          networkType=opt$network_type, verbose=2)
png(opt$out_plot_soft, width=1200, height=600, res=120)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Power", ylab=expression(R^2), type="n", main="Scale-free topology")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="steelblue")
abline(h=0.85, col="firebrick", lty=2)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Power", ylab="Mean connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="firebrick")
par(mfrow=c(1,1)); dev.off()
soft_power <- if (opt$soft_power>0) opt$soft_power else sft$powerEstimate
if (is.na(soft_power)||is.null(soft_power)) { soft_power <- 6; warning("Using power=6") }
cat("Power selected:", soft_power, "\n")

# =============================================================================
# STEP 6: Network construction
# =============================================================================
cat("\n=== STEP 6: Network construction ===\n")
cor <- WGCNA::cor
net <- blockwiseModules(
  expr_matrix, power=soft_power,
  networkType=opt$network_type, TOMType=opt$network_type,
  minModuleSize=opt$min_module_size, reassignThreshold=0,
  mergeCutHeight=opt$merge_cut_height, numericLabels=TRUE,
  pamRespectsDendro=FALSE, maxBlockSize=opt$max_block_size,
  saveTOMs=FALSE, verbose=3
)
cor <- stats::cor
module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_matrix)
module_colors_all    <- module_colors
cat("Modules:", length(unique(module_colors))-1, "\n")
print(table(module_colors))
png(opt$out_plot_dendro, width=1400, height=800, res=120)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                    "Module", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    main="Gene dendrogram (block 1)")
dev.off()
write.table(data.frame(Gene=colnames(expr_matrix), Module=module_colors),
            file=opt$out_table_modules, sep="\t", row.names=FALSE, quote=FALSE)

# =============================================================================
# STEP 7: Module-trait correlation
# =============================================================================
cat("\n=== STEP 7: Module-trait correlations ===\n")
MEs               <- orderMEs(net$MEs)
module_trait_cor  <- cor(MEs, trait_data, use="p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_matrix))
txt <- paste0(signif(module_trait_cor,2),"\n(",signif(module_trait_pval,1),")")
dim(txt) <- dim(module_trait_cor)
png(opt$out_plot_heatmap_all,
    width=max(600,200*ncol(module_trait_cor)),
    height=max(800,15*nrow(module_trait_cor)), res=120)
par(mar=c(6,10,3,3))
labeledHeatmap(Matrix=module_trait_cor, xLabels=colnames(trait_data),
               yLabels=names(MEs), ySymbols=names(MEs),
               colorLabels=FALSE, colors=blueWhiteRed(50),
               textMatrix=txt, setStdMargins=FALSE,
               cex.text=0.5, zlim=c(-1,1),
               main="Module eigengene – trait correlation")
dev.off()

# =============================================================================
# STEP 8: Top module — GS and kME
# =============================================================================
cat("\n=== STEP 8: Top module ===\n")
best_ME_name  <- names(MEs)[which.max(abs(module_trait_cor[,trait_name]))]
best_ME_num   <- as.numeric(gsub("ME","",best_ME_name))
module_color  <- labels2colors(best_ME_num)
genes_in_mod  <- colnames(expr_matrix)[net$colors==best_ME_num]
cat("Top module:", best_ME_name, "(", module_color, ")  r =",
    round(module_trait_cor[best_ME_name,trait_name],4), "\n")

gene_sig      <- setNames(as.numeric(cor(expr_matrix,trait_data[[trait_name]],use="p")),
                           colnames(expr_matrix))
gene_sig_pval <- corPvalueStudent(gene_sig, nrow(expr_matrix))
MM      <- as.numeric(cor(expr_matrix[,genes_in_mod], MEs[,best_ME_name], use="p"))
MM_pval <- corPvalueStudent(MM, nrow(expr_matrix))
module_df <- data.frame(Gene=genes_in_mod, ModuleMembership=MM, MM_pvalue=MM_pval,
                         GeneSignificance=gene_sig[genes_in_mod],
                         GS_pvalue=gene_sig_pval[genes_in_mod])
module_df <- module_df[order(-abs(module_df$ModuleMembership)),]
write.table(module_df, file=opt$out_table_gs_mm, sep="\t", row.names=FALSE, quote=FALSE)
col_safe <- tryCatch(adjustcolor(module_color,.5), error=function(e)"steelblue")
png(opt$out_plot_mm_gs, width=900, height=700, res=120)
plot(abs(module_df$ModuleMembership), abs(module_df$GeneSignificance),
     xlab=paste("kME —",module_color), ylab=paste("GS for",trait_name),
     main=paste("MM vs GS —",module_color), pch=20, col=col_safe)
abline(lm(abs(GeneSignificance)~abs(ModuleMembership),data=module_df),
       col="firebrick",lwd=2)
legend("topleft",bty="n",
       legend=paste("r =",round(cor(abs(module_df$ModuleMembership),
                                     abs(module_df$GeneSignificance)),3)))
dev.off()

# =============================================================================
# STEP 9: Hub genes
# =============================================================================
cat("\n=== STEP 9: Hub genes ===\n")
adj_m        <- adjacency(expr_matrix[,genes_in_mod], power=soft_power, type=opt$network_type)
conn         <- colSums(adj_m)-1
hub_info     <- data.frame(Gene=genes_in_mod, Connectivity=conn,
                            ModuleMembership=MM,
                            GeneSignificance=gene_sig[genes_in_mod])
hub_info     <- hub_info[order(-hub_info$Connectivity),]
hub_info$isHub <- hub_info$Connectivity > quantile(hub_info$Connectivity,.9)
write.table(hub_info[hub_info$isHub,], file=opt$out_table_hubs,
            sep="\t",row.names=FALSE,quote=FALSE)
png(opt$out_plot_hubs, width=900, height=700, res=120)
plot(hub_info$ModuleMembership, hub_info$Connectivity,
     xlab="kME", ylab="Connectivity", main=paste("Hub genes —",module_color),
     pch=20, col=ifelse(hub_info$isHub,"firebrick","grey60"))
legend("topleft",legend=c("Hub (top 10%)","Other"),
       col=c("firebrick","grey60"),pch=20,bty="n")
dev.off()

# =============================================================================
# STEP 10: kME filtering
# =============================================================================
genes_high_kME <- module_df[abs(module_df$ModuleMembership)>opt$kme_threshold,]
cat("|kME| >", opt$kme_threshold, ": before=", nrow(module_df),
    " after=", nrow(genes_high_kME), "\n")

# =============================================================================
# STEP 11: Cytoscape export
# =============================================================================
cat("\n=== STEP 11: Cytoscape export ===\n")
module_colors_all <- labels2colors(net$colors)
names(module_colors_all) <- colnames(expr_matrix)
module_trait_cor  <- cor(MEs, trait_data, use="p")
unique_mods <- sort(unique(net$colors)); unique_mods <- unique_mods[unique_mods!=0]
mod_summary <- data.frame(
  ModuleNum   = unique_mods,
  ModuleColor = labels2colors(unique_mods),
  Correlation = module_trait_cor[paste0("ME",unique_mods), trait_name]
)
modules_selected <- mod_summary[abs(mod_summary$Correlation)>=opt$cor_threshold,]
if (!nrow(modules_selected)) {
  warning("No modules at cor_threshold=",opt$cor_threshold,"; lowering to 0.3")
  modules_selected <- mod_summary[abs(mod_summary$Correlation)>=0.3,]
}
modules_selected <- modules_selected[order(-abs(modules_selected$Correlation)),]
write.table(modules_selected,file=opt$out_table_sel_modules,
            sep="\t",row.names=FALSE,quote=FALSE)
selected_modules <- modules_selected$ModuleColor

for (i in seq_len(nrow(modules_selected))) {
  mn  <- modules_selected$ModuleNum[i]; mc <- modules_selected$ModuleColor[i]
  mec <- paste0("ME",mn)
  cat("  Module:", mc, "\n")
  gm <- colnames(expr_matrix)[net$colors==mn]
  if (!mec %in% colnames(MEs)) next
  kME_i <- setNames(as.numeric(cor(expr_matrix[,gm],MEs[,mec],use="p")), gm)
  gk    <- names(kME_i)[abs(kME_i)>opt$kme_threshold]
  if (length(gk)<3) next
  adj_i <- adjacency(expr_matrix[,gk], power=soft_power, type=opt$network_type)
  TOM_i <- TOMsimilarity(adj_i, TOMType=opt$network_type, verbose=0)
  rownames(TOM_i) <- colnames(TOM_i) <- gk
  tc   <- quantile(TOM_i[upper.tri(TOM_i)], opt$tom_percentile)
  idx  <- which(upper.tri(TOM_i), arr.ind=TRUE)
  edf  <- data.frame(fromNode=gk[idx[,1]], toNode=gk[idx[,2]], weight=TOM_i[idx])
  edf  <- edf[edf$weight>=tc,]
  ndf  <- data.frame(node=gk, module=mc, kME=round(kME_i[gk],4),
                      GeneSignificance=round(gene_sig[gk],4),
                      Connectivity=round(colSums(adj_i)-1,4))
  write.table(edf, file.path(opt$out_cytoscape_dir,paste0("edges_",mc,"_kMEfilt.txt")),
              sep="\t",row.names=FALSE,quote=FALSE)
  write.table(ndf, file.path(opt$out_cytoscape_dir,paste0("nodes_",mc,"_kMEfilt.txt")),
              sep="\t",row.names=FALSE,quote=FALSE)
}
cat("Cytoscape export complete.\n")

# =============================================================================
# STEP 12: Selected module heatmap
# =============================================================================
cat("\n=== STEP 12: Selected module heatmap ===\n")
sel_ME <- paste0("ME",modules_selected$ModuleNum)
sel_ME <- sel_ME[sel_ME %in% rownames(module_trait_cor)]
if (length(sel_ME)>0) {
  sc <- module_trait_cor[sel_ME,,drop=FALSE]
  sp <- module_trait_pval[sel_ME,,drop=FALSE]
  st <- paste0(signif(sc,2),"\n(",signif(sp,1),")")
  dim(st) <- dim(sc)
  sl <- modules_selected$ModuleColor[match(sel_ME,paste0("ME",modules_selected$ModuleNum))]
  png(opt$out_plot_heatmap_sel, width=max(500,220*ncol(sc)),
      height=max(400,55*nrow(sc)), res=120)
  par(mar=c(6,12,3,3))
  labeledHeatmap(Matrix=sc, xLabels=colnames(trait_data),
                 yLabels=sl, ySymbols=sel_ME,
                 colorLabels=FALSE, colors=blueWhiteRed(50),
                 textMatrix=st, setStdMargins=FALSE, cex.text=0.8, zlim=c(-1,1),
                 main=paste("Selected modules (|r| >=",opt$cor_threshold,")"))
  dev.off()
} else placeholder_png(opt$out_plot_heatmap_sel, "No modules above threshold")

# =============================================================================
# STEP 13 (optional): DEG overlap
# =============================================================================
has_degs <- opt$up_degs!="None" && opt$down_degs!="None" &&
  file.exists(opt$up_degs) && file.exists(opt$down_degs)
if (has_degs) {
  cat("\n=== STEP 13: DEG overlap ===\n")
  ug <- unique(na.omit(as.character(read.table(opt$up_degs,  header=FALSE,sep="\t")[,1])))
  dg <- unique(na.omit(as.character(read.table(opt$down_degs,header=FALSE,sep="\t")[,1])))
  ug <- ug[ug!=""]; dg <- dg[dg!=""]
  degs <- data.frame(Gene=c(ug,dg),
                     Status=c(rep("Up",length(ug)),rep("Down",length(dg))))
  degs <- degs[!duplicated(degs$Gene),]
  ns <- do.call(rbind, lapply(selected_modules, function(mod) {
    nf <- file.path(opt$out_cytoscape_dir, paste0("nodes_",mod,"_kMEfilt.txt"))
    if (!file.exists(nf)) return(NULL)
    gin <- read.table(nf,header=TRUE,sep="\t")$node
    data.frame(Module=mod, NetworkNodes=length(gin),
               DEGs_total=length(intersect(gin,degs$Gene)),
               DEGs_up   =length(intersect(gin,degs$Gene[degs$Status=="Up"])),
               DEGs_down =length(intersect(gin,degs$Gene[degs$Status=="Down"])))
  }))
  if (!is.null(ns) && nrow(ns)>0) {
    ns$Proportion_DEGs <- ns$DEGs_total/ns$NetworkNodes
    write.table(ns, file=opt$out_table_deg_summary,
                sep="\t",row.names=FALSE,quote=FALSE)
    for (mod in selected_modules) {
      nf <- file.path(opt$out_cytoscape_dir, paste0("nodes_",mod,"_kMEfilt.txt"))
      if (!file.exists(nf)) next
      nd <- read.table(nf,header=TRUE,sep="\t")
      nd$DE_status <- "Not_DE"
      nd$DE_status[nd$node %in% degs$Gene[degs$Status=="Up"]]   <- "Up"
      nd$DE_status[nd$node %in% degs$Gene[degs$Status=="Down"]] <- "Down"
      write.table(nd, file.path(opt$out_cytoscape_dir,
                                 paste0("nodes_",mod,"_kMEfilt_with_DEG.txt")),
                  sep="\t",row.names=FALSE,quote=FALSE)
    }
    pd <- melt(ns[,c("Module","DEGs_up","DEGs_down")],
               id.vars="Module",variable.name="Direction",value.name="Count")
    pd$Direction <- factor(pd$Direction, c("DEGs_up","DEGs_down"),
                            c("Up-regulated","Down-regulated"))
    ggsave(opt$out_plot_deg_bar,
           ggplot(pd,aes(x=Module,y=Count,fill=Direction))+
             geom_bar(stat="identity",position="dodge")+
             scale_fill_manual(values=c("Up-regulated"="firebrick",
                                        "Down-regulated"="steelblue"))+
             labs(title="DEGs per module",x="Module",y="Count")+
             theme_minimal(base_size=13)+
             theme(axis.text.x=element_text(angle=45,hjust=1)),
           width=max(8,nrow(ns)*0.8),height=5,dpi=150)
  } else {
    placeholder_tsv(opt$out_table_deg_summary,"No node files found")
    placeholder_png(opt$out_plot_deg_bar,"No DEG data matched network")
  }
} else {
  placeholder_tsv(opt$out_table_deg_summary,"DEG analysis not requested")
  placeholder_png(opt$out_plot_deg_bar,"DEG analysis not requested")
}

# =============================================================================
# STEP 14 (optional): GO / Functional Enrichment Analysis
# =============================================================================
run_enrichment <- tolower(trimws(opt$run_enrichment)) == "yes"

if (run_enrichment) {
  cat("\n=== STEP 14: Gene Ontology / Enrichment Analysis ===\n")

  # Conditionally load enrichment packages
  for (pkg in c("clusterProfiler","enrichplot")) {
    if (!requireNamespace(pkg, quietly=TRUE))
      stop("Package '", pkg, "' required. Install: BiocManager::install('", pkg, "')")
    suppressPackageStartupMessages(library(pkg, character.only=TRUE))
  }

  org_pkg <- trimws(opt$orgdb_package)
  if (org_pkg == "none" || org_pkg == "" || org_pkg == "None")
    stop("--orgdb_package required when --run_enrichment yes")
  if (!requireNamespace(org_pkg, quietly=TRUE))
    stop("OrgDb package '", org_pkg, "' not installed. ",
         "Install: BiocManager::install('", org_pkg, "')")
  suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
  orgdb <- get(org_pkg)
  cat("OrgDb:", org_pkg, "  keyType:", opt$gene_id_type, "\n")

  universe_genes <- colnames(expr_matrix)
  all_go  <- list()
  all_keg <- list()
  run_kegg_flag <- tolower(trimws(opt$run_kegg)) == "yes"

  for (mod in selected_modules) {
    cat("\n  Processing module:", mod, "\n")
    mn_i    <- modules_selected$ModuleNum[modules_selected$ModuleColor == mod]
    genes_i <- colnames(expr_matrix)[net$colors == mn_i]
    if (length(genes_i) < opt$enrich_min_gs) {
      cat("  Too few genes, skip.\n"); next
    }

    # ── GO: BP, MF, CC ─────────────────────────────────────────────────────
    for (ont in c("BP","MF","CC")) {
      cat("    GO", ont, "... ")
      ego <- tryCatch(
        enrichGO(gene=genes_i, universe=universe_genes, OrgDb=orgdb,
                 keyType=opt$gene_id_type, ont=ont,
                 pAdjustMethod="BH",
                 pvalueCutoff=opt$enrich_pval_cut,
                 qvalueCutoff=opt$enrich_qval_cut,
                 minGSSize=opt$enrich_min_gs, maxGSSize=opt$enrich_max_gs,
                 readable=FALSE),
        error=function(e) { cat("ERROR:", e$message, "\n"); NULL }
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        df <- as.data.frame(ego); df$Module <- mod; df$Ontology <- ont
        all_go[[paste(mod,ont,sep="_")]] <- df
        write.table(df,
          file.path(opt$out_enrich_dir, paste0("GO_",ont,"_",mod,".tsv")),
          sep="\t",row.names=FALSE,quote=FALSE)
        cat(nrow(df), "terms\n")
      } else cat("0 terms\n")
    }

    # ── KEGG ───────────────────────────────────────────────────────────────
    if (run_kegg_flag) {
      if (nchar(trimws(opt$kegg_organism)) == 0)
        stop("--kegg_organism required when --run_kegg yes")
      cat("    KEGG", opt$kegg_organism, "... ")
      ek <- tryCatch(
        enrichKEGG(gene=genes_i, organism=opt$kegg_organism,
                   pAdjustMethod="BH",
                   pvalueCutoff=opt$enrich_pval_cut,
                   qvalueCutoff=opt$enrich_qval_cut,
                   minGSSize=opt$enrich_min_gs, maxGSSize=opt$enrich_max_gs),
        error=function(e) { cat("ERROR:", e$message,"\n"); NULL }
      )
      if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
        dk <- as.data.frame(ek); dk$Module <- mod
        all_keg[[mod]] <- dk
        write.table(dk,
          file.path(opt$out_enrich_dir, paste0("KEGG_",mod,".tsv")),
          sep="\t",row.names=FALSE,quote=FALSE)
        cat(nrow(dk),"pathways\n")
      } else cat("0 pathways\n")
    }
  }

  # ── Combined summary ───────────────────────────────────────────────────────
  if (length(all_go) > 0) {
    go_all <- do.call(rbind, all_go)
    write.table(go_all, file=opt$out_table_go_summary,
                sep="\t",row.names=FALSE,quote=FALSE)
    cat("GO summary:", nrow(go_all), "terms across", length(all_go), "module-ont pairs\n")

    # Figures based on GO-BP only (most informative)
    bp_list <- all_go[grepl("_BP$",names(all_go))]

    if (length(bp_list) > 0) {
      top <- do.call(rbind, lapply(bp_list, function(df) {
        df <- df[order(df$p.adjust),]; head(df, 8)
      }))
      top$GR_num <- sapply(top$GeneRatio, function(x) {
        v <- as.numeric(strsplit(x,"/")[[1]]); v[1]/v[2]
      })
      top$neg_log_padj <- -log10(top$p.adjust + 1e-300)
      top$Description  <- factor(top$Description,
                                  levels=rev(unique(top$Description[order(top$p.adjust)])))
      nmods <- length(unique(top$Module))

      # Dot plot
      p_dot <- ggplot(top, aes(x=GR_num, y=Description,
                                size=Count, color=neg_log_padj)) +
        geom_point(alpha=0.9) +
        scale_color_gradient(low="steelblue", high="firebrick",
                              name=expression(-log[10](p.adj))) +
        scale_size_continuous(range=c(2,9), name="Gene count") +
        facet_wrap(~Module, scales="free_y", ncol=min(nmods,3)) +
        labs(title="GO Biological Process (top 8 per module)",
             x="Gene ratio", y=NULL) +
        theme_bw(base_size=11) +
        theme(strip.background=element_rect(fill="grey90"),
              axis.text.y=element_text(size=8),
              plot.title=element_text(face="bold"))
      ggsave(opt$out_plot_go_dot, p_dot,
             width=min(22,8+nmods*3), height=min(26,5+nrow(top)*0.35),
             dpi=150, limitsize=FALSE)

      # Bar plot
      p_bar <- ggplot(top, aes(x=reorder(Description,neg_log_padj),
                                y=neg_log_padj, fill=Module)) +
        geom_col(show.legend=FALSE) + coord_flip() +
        facet_wrap(~Module, scales="free_y", ncol=min(nmods,3)) +
        labs(title="GO-BP enrichment (top 8 per module)",
             x=NULL, y=expression(-log[10](p.adj))) +
        theme_bw(base_size=11) +
        theme(strip.background=element_rect(fill="grey90"),
              axis.text.y=element_text(size=8),
              plot.title=element_text(face="bold"))
      ggsave(opt$out_plot_go_bar, p_bar,
             width=min(22,8+nmods*3), height=min(26,5+nrow(top)*0.35),
             dpi=150, limitsize=FALSE)
      cat("GO plots saved.\n")
    } else {
      placeholder_png(opt$out_plot_go_dot, "No GO-BP terms found")
      placeholder_png(opt$out_plot_go_bar, "No GO-BP terms found")
    }
  } else {
    cat("No significant GO terms found.\n")
    placeholder_tsv(opt$out_table_go_summary, "No significant GO terms")
    placeholder_png(opt$out_plot_go_dot, "No significant GO terms")
    placeholder_png(opt$out_plot_go_bar, "No significant GO terms")
  }

} else {
  cat("\nEnrichment not requested — skipping STEP 14.\n")
  placeholder_tsv(opt$out_table_go_summary, "GO enrichment not requested")
  placeholder_png(opt$out_plot_go_dot, "GO enrichment not requested")
  placeholder_png(opt$out_plot_go_bar, "GO enrichment not requested")
}

# =============================================================================
# STEP 15 (optional): Module Preservation / Z-summary
#
# Compares co-expression structure of modules defined in the QUERY (main/test)
# network against a REFERENCE network (e.g. healthy vs pathogen-infected).
#
# Uses WGCNA::modulePreservation() — Langfelder & Horvath 2011.
# Zsummary interpretation:
#   < 2    : not preserved (structure differs between conditions)
#   2 – 10 : moderately preserved
#   > 10   : highly preserved (same module exists in both organisms/conditions)
# =============================================================================
run_preservation <- tolower(trimws(opt$run_preservation)) == "yes"

if (run_preservation) {
  cat("\n=== STEP 15: Module Preservation / Z-summary ===\n")

  for (v in c("ref_vst_matrix","ref_sample_info")) {
    val <- opt[[v]]
    if (is.null(val) || val=="None")
      stop("--", v, " is required for preservation analysis.")
    if (!file.exists(val)) stop(v, " file not found: ", val)
  }

  # Load reference data
  cat("Loading reference dataset...\n")
  ref_si <- load_tabular(opt$ref_sample_info, lbl="reference sample info")
  if (!opt$ref_sample_col %in% colnames(ref_si))
    stop("Column '", opt$ref_sample_col, "' not in reference metadata.")
  rownames(ref_si) <- ref_si[[opt$ref_sample_col]]

  ref_expr <- load_tabular(opt$ref_vst_matrix, row1=TRUE, lbl="reference VST")
  cat("Reference dimensions:", dim(ref_expr), "\n")

  ref_common <- intersect(colnames(ref_expr), rownames(ref_si))
  if (!length(ref_common)) stop("No common samples in reference.")
  ref_expr <- ref_expr[, ref_common, drop=FALSE]
  cat("Reference samples:", length(ref_common), "\n")

  # Shared genes between test and reference
  shared <- intersect(colnames(expr_matrix), rownames(ref_expr))
  cat("Shared genes:", length(shared), "\n")
  if (length(shared) < 50)
    stop("Too few shared genes (", length(shared), "). Check gene ID format.")

  # Build matrices restricted to shared genes
  test_sub <- expr_matrix[, shared, drop=FALSE]
  ref_mat  <- t(as.matrix(ref_expr[shared, , drop=FALSE]))
  cols_sub <- module_colors[shared]

  # Optional stratified subsampling for large datasets
  max_g <- opt$pres_max_genes
  if (ncol(test_sub) > max_g) {
    cat("Subsetting to", max_g, "genes (stratified by module).\n")
    set.seed(opt$pres_random_seed)
    mod_list <- split(names(cols_sub), cols_sub)
    per_mod  <- max(5L, max_g %/% length(mod_list))
    keep     <- unique(unlist(lapply(mod_list, function(g) {
      if (length(g) > per_mod) sample(g, per_mod) else g
    })))
    test_sub <- test_sub[, keep, drop=FALSE]
    ref_mat  <- ref_mat[,  keep, drop=FALSE]
    cols_sub <- cols_sub[keep]
  }
  cat("Test matrix for preservation:", dim(test_sub), "\n")
  cat("Reference matrix:", dim(ref_mat), "\n")
  cat("Running modulePreservation (", opt$pres_n_perms, "permutations) — may take minutes...\n")

  set.seed(opt$pres_random_seed)
  mp <- tryCatch(
    modulePreservation(
      multiData  = list(Test=list(data=test_sub), Reference=list(data=ref_mat)),
      multiColor = list(Test=cols_sub),
      dataIsExpr        = TRUE,
      referenceNetworks = 1,
      nPermutations     = opt$pres_n_perms,
      randomSeed        = opt$pres_random_seed,
      quickCor          = 0,
      verbose           = 2
    ),
    error=function(e) { cat("modulePreservation error:", e$message, "\n"); NULL }
  )

  if (is.null(mp)) {
    for (p in c(opt$out_table_preservation, opt$out_table_pres_class))
      placeholder_tsv(p, "modulePreservation failed")
    for (p in c(opt$out_plot_zsummary, opt$out_plot_pres_heatmap))
      placeholder_png(p, "modulePreservation failed")
  } else {

    # Extract preservation statistics
    # Z-statistics: mp$preservation$Z[[ref]][[testInRef]]
    Z_stats  <- mp$preservation$Z[[1]][[2]]
    obs_stat <- mp$preservation$observed[[1]][[2]]

    mod_sizes <- table(cols_sub)

    pres_df <- data.frame(
      Module        = rownames(Z_stats),
      ModuleSize    = as.integer(mod_sizes[rownames(Z_stats)]),
      Zsummary      = round(Z_stats$Zsummary.pres, 3),
      Zdensity      = round(Z_stats$Zdensity.pres, 3),
      Zconnectivity = round(Z_stats$Zconnectivity.pres, 3),
      MedianRank    = round(obs_stat$medianRank.pres, 3),
      stringsAsFactors=FALSE
    )
    pres_df <- pres_df[pres_df$Module != "grey", ]
    pres_df <- pres_df[order(-pres_df$Zsummary), ]
    pres_df$Preservation <- cut(pres_df$Zsummary,
      breaks=c(-Inf, 2, 10, Inf),
      labels=c("Not preserved","Moderately preserved","Highly preserved"))

    write.table(pres_df, file=opt$out_table_preservation,
                sep="\t",row.names=FALSE,quote=FALSE)
    write.table(pres_df[,c("Module","ModuleSize","Zsummary","Preservation")],
                file=opt$out_table_pres_class,
                sep="\t",row.names=FALSE,quote=FALSE)
    cat("Preservation table saved.\n")
    cat("Summary:\n"); print(table(pres_df$Preservation))

    # ── Z-summary scatter plot ────────────────────────────────────────────
    pres_plot <- pres_df[is.finite(pres_df$Zsummary) & !is.na(pres_df$ModuleSize),]

    pcolors <- c("Not preserved"       ="firebrick",
                  "Moderately preserved"="#F39C12",
                  "Highly preserved"    ="#27AE60")

    p_z <- ggplot(pres_plot,
                  aes(x=log10(ModuleSize+1), y=Zsummary,
                      color=Preservation, label=Module)) +
      geom_hline(yintercept=c(2,10), linetype="dashed",
                 color=c("firebrick","#27AE60"), linewidth=0.8) +
      geom_point(aes(size=ModuleSize), alpha=0.85) +
      geom_text_repel(size=3.2, max.overlaps=25,
                      segment.color="grey60", box.padding=0.35) +
      scale_color_manual(values=pcolors, name="Preservation") +
      scale_size_continuous(range=c(3,13), name="Module size") +
      annotate("text", y=c(1.0,5.5,11.5),
               x=rep(min(log10(pres_plot$ModuleSize+1)), 3),
               label=c("Not preserved (Z < 2)",
                       "Moderately preserved (2 ≤ Z ≤ 10)",
                       "Highly preserved (Z > 10)"),
               hjust=0, size=3.2,
               color=c("firebrick","#E67E22","#27AE60")) +
      labs(
        title    = "Module Preservation: Z-summary",
        subtitle = "Query network modules tested in the Reference network",
        x        = expression(log[10]("module size")),
        y        = "Zsummary"
      ) +
      theme_bw(base_size=13) +
      theme(legend.position="bottom",
            plot.title=element_text(face="bold"),
            plot.subtitle=element_text(color="grey40"))
    ggsave(opt$out_plot_zsummary, p_z, width=10, height=7, dpi=150)
    cat("Z-summary plot saved.\n")

    # ── Multi-statistic heatmap ───────────────────────────────────────────
    hvars <- c("Zsummary","Zdensity","Zconnectivity","MedianRank")
    hvars <- hvars[hvars %in% colnames(pres_df)]
    hmat  <- as.matrix(pres_df[, hvars]); rownames(hmat) <- pres_df$Module
    hmat_cap <- hmat; hmat_cap[hmat_cap > 15] <- 15
    hmat_cap[is.na(hmat_cap)] <- 0

    row_ann <- data.frame(Preservation=pres_df$Preservation,
                           row.names=pres_df$Module)
    ann_colors <- list(Preservation=pcolors)

    png(opt$out_plot_pres_heatmap,
        width=900, height=max(400, 32*nrow(hmat_cap)), res=120)
    pheatmap(hmat_cap,
      annotation_row=row_ann, annotation_colors=ann_colors,
      color=colorRampPalette(c("white","#3498DB","#1A5276"))(100),
      cluster_cols=FALSE, fontsize_row=9, fontsize_col=11,
      main="Module preservation statistics\n(Z-scores, capped at 15)",
      display_numbers=round(hmat,2), number_format="%.1f", fontsize_number=7)
    dev.off()
    cat("Preservation heatmap saved.\n")
    cat("\nTop 10 modules by preservation:\n")
    print(head(pres_df[,c("Module","ModuleSize","Zsummary","Preservation")],10))
  }

} else {
  cat("\nPreservation not requested — skipping STEP 15.\n")
  placeholder_tsv(opt$out_table_preservation, "Preservation not requested")
  placeholder_tsv(opt$out_table_pres_class,   "Preservation not requested")
  placeholder_png(opt$out_plot_zsummary,       "Preservation not requested")
  placeholder_png(opt$out_plot_pres_heatmap,   "Preservation not requested")
}

# =============================================================================
# Save workspace
# =============================================================================
cat("\n=== Saving RData ===\n")
objs <- c("net","expr_matrix","trait_data","MEs","module_colors","module_colors_all",
          "module_trait_cor","module_trait_pval","gene_sig","module_df",
          "modules_selected","selected_modules","soft_power")
save(list=objs[sapply(objs,exists)], file=opt$out_rdata)
cat("Saved:", opt$out_rdata, "\n")
cat("\n=== WGCNA v2 complete! ===\n")
sessionInfo()
