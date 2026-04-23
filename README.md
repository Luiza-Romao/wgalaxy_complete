# WGCNA Galaxy Tool 

**Weighted Gene Co-expression Network Analysis** for the Galaxy bioinformatics platform.

**Recent improvements:**
- ✅ Gene Ontology / functional enrichment analysis on WGCNA modules (clusterProfiler)
- ✅ Module preservation / Z-summary to compare networks between conditions or organisms

---

## File Structure

| File | Description |
|------|-------------|
| `wgcna.xml` | Galaxy tool wrapper — inputs, outputs, UI, help text |
| `macros.xml` | Reusable XML macros (requirements, parameter sections) |
| `wgcna.R` | Main R analysis script |
| `conda_environment.yml` | Conda env with all dependencies |
| `generate_test_data.R` | Generates minimal test datasets |

---

## Installation

### 1. Copy tool files
```bash
cp -r wgcna_galaxy/ /galaxy/tools/wgcna_galaxy/
```

### 2. Register in Galaxy config
Add to `tool_conf.xml`:
```xml
<section id="transcriptomics" name="Transcriptomics">
    <tool file="wgcna_galaxy/wgcna.xml" />
</section>
```

### 3. Create conda environment
```bash
conda env create -f wgcna_galaxy/conda_environment.yml
```
For additional OrgDb packages (e.g. rat, yeast, zebrafish):
```bash
conda activate wgcna_galaxy_v2
Rscript -e "BiocManager::install(c('org.Rn.eg.db','org.Sc.sgd.db','org.Dr.eg.db'))"
```

### 4. Restart Galaxy

---

## Input Formats

### VST Expression Matrix (tabular)
```
GeneID       Sample_A1   Sample_A2   Sample_B1   Sample_B2
AT1G01010    8.23        7.95        5.12        5.44
AT1G01020    6.10        6.35        9.87        10.02
```
- Tab-separated; first column = gene IDs; header = sample IDs

### Sample Metadata (tabular)
```
SampleID     Treatment
Sample_A1    Infected
Sample_A2    Infected
Sample_B1    Healthy
Sample_B2    Healthy
```
- Tab-separated
- Header required
- Column names are configurable in the tool
  
### DEG Lists (tabular, optional)
```
AT1G01010
AT2G05230
```
- Tab-separated, no header required
- First column = gene IDs (must match expression matrix)
- Provide separately for up- and down-regulated genes

### Aditional VST Matrix (for Module preservation / Z-summary to compare networks between conditions or organisms, tabular)
Same format as the main VST matrix but from a different condition or organism.

---

## Pipeline steps

```
INPUT: VST matrix + sample metadata
   │
   ├─  1. Load & align samples
   ├─  2. Filter genes (presence + variance)
   ├─  3. Build binary trait matrix
   ├─  4. Sample clustering (outlier detection)
   ├─  5. Soft-threshold power selection
   ├─  6. Blockwise TOM network construction
   ├─  7. Module eigengenes + module-trait Pearson r
   ├─  8. Gene significance (GS) + kME — top module
   ├─  9. Hub gene identification
   ├─ 10. kME filtering (core module genes)
   ├─ 11. Cytoscape export (edge + node tables)
   ├─ 12. Selected module heatmap
   ├─ 13. DEG overlap analysis             [OPTIONAL]
   ├─ 14. GO / KEGG enrichment             [OPTIONAL]
   └─ 15. Module preservation / Z-summary  [OPTIONAL]
```

---

## WGCNA tool Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Presence fraction | 0.60 | Min fraction of samples with non-zero expression |
| Variance percentile | 40 | Remove genes below this variance percentile |
| Soft power (β) | 0 (auto) | Network adjacency power; 0 = auto-detect |
| Network type | signed | `signed` preserves co-expression direction |
| Min module size | 50 | Minimum genes per co-expression module |
| Merge cut height | 0.25 | Modules with eigengene r ≥ 0.75 are merged |
| kME threshold | 0.80 | Module membership filter for core genes |
| Correlation threshold | 0.60 | Pearson’s correlation coefficient cutoff for Cytoscape module export |
| TOM percentile | 0.95 | Top 5% TOM edges exported to Cytoscape |

---
## All Outputs Reference

### Plots (PNG)
| File | Description |
|------|-------------|
| Sample clustering | Hierarchical dendrogram + trait color bar |
| Soft-threshold | R² and mean connectivity vs β |
| Gene dendrogram | Module color bar |
| Full heatmap | All modules × traits |
| MM vs GS | Module membership vs gene significance |
| Hub genes | Connectivity vs kME |
| Selected heatmap | High-correlation modules only |
| DEG barplot | Up/down DEGs per network module |
| **GO dot plot** | Top GO-BP terms, dot-sized by gene count |
| **GO bar plot** | Enrichment by –log10(p.adj) |
| **Z-summary scatter** | Preservation score per module |
| **Preservation heatmap** | Zsummary, Zdensity, Zconnectivity |

### Tables (TSV)
| File | Description |
|------|-------------|
| Module assignments | Gene → module color |
| GS/MM | Gene significance + kME (top module) |
| Hub genes | Top 10% connectivity genes |
| Selected modules | High-correlation module summary |
| DEG summary | DEG counts per network module |
| **GO summary** | All enriched GO terms (combined) |
| **GO/KEGG per module** | Per-module enrichment (collection) |
| **Preservation full** | All Z-statistics |
| **Preservation class** | Not/Moderate/Highly preserved |

### Cytoscape Networks (collection of tabular files)
- `edges_{module}_kMEfilt.txt` — Edge table: `fromNode`, `toNode`, `weight` (TOM)
- `nodes_{module}_kMEfilt.txt` — Node table: `node`, `module`, `kME`, `GeneSignificance`, `Connectivity`
- `nodes_{module}_kMEfilt_with_DEG.txt` *(optional)* — Nodes with `DE_status` column

### Other
- `wgcna_results.RData` — All R objects for downstream custom analysis
- Analysis log — stdout/stderr

---
### Gene Ontology (GO) clusterProfiler

### What it does
Runs `enrichGO()` (clusterProfiler) across all three GO ontologies
(Biological Process, Molecular Function, Cellular Component) and optionally
`enrichKEGG()` for each module that passes the correlation threshold.


| Parameter | Default | Description |
|-----------|---------|-------------|
| OrgDb package | — | Bioconductor annotation DB for your organism |
| Gene ID type | SYMBOL | Must match your expression matrix row names |
| KEGG enrichment | No | Run KEGG in addition to GO |
| KEGG organism code | — | 3-letter code e.g. `hsa`, `mmu`, `ath`, `sce` |
| Adjusted p-value cutoff | 0.05 | BH-corrected significance threshold |
| q-value cutoff | 0.20 | Additional q-value filter |
| Min gene-set size | 10 | Discard very small GO terms |
| Max gene-set size | 500 | Discard very large GO terms |

### Supported organisms (OrgDb)

| Organism | Package |
|----------|---------|
| Homo sapiens | org.Hs.eg.db |
| Mus musculus | org.Mm.eg.db |
| Rattus norvegicus | org.Rn.eg.db |
| Arabidopsis thaliana | org.At.tair.db |
| Saccharomyces cerevisiae | org.Sc.sgd.db |
| Caenorhabditis elegans | org.Ce.eg.db |
| Drosophila melanogaster | org.Dm.eg.db |
| Danio rerio | org.Dr.eg.db |
| Sus scrofa | org.Ss.eg.db |
| Bos taurus | org.Bt.eg.db |
| Gallus gallus | org.Gg.eg.db |
| Custom | Any OrgDb in the "custom" field |

### Outputs
- `table_go_summary.tsv` — all significant GO terms across all modules and ontologies
- `enrichment_output/GO_BP_<module>.tsv` — per-module GO-BP results
- `enrichment_output/GO_MF_<module>.tsv` — per-module GO-MF results
- `enrichment_output/GO_CC_<module>.tsv` — per-module GO-CC results
- `enrichment_output/KEGG_<module>.tsv` — per-module KEGG results (if enabled)
- `plot_go_dotplot.png` — dot plot (top 8 BP terms per module, faceted)
- `plot_go_barplot.png` — bar plot (–log10 p.adj per module, faceted)

---

## Module Preservation / Z-summary 

### What it does
Uses WGCNA's `modulePreservation()` to test whether modules identified in the
**query** (test) network also exist as coherent co-expression units in a
**reference** network.

**Typical use cases:**
- Pathogen-infected vs. healthy plants/animals
- Same tissue, different developmental stages
- Conserved modules across species (e.g. human vs. mouse)

### Z-summary interpretation

| Zsummary | Level | Interpretation |
|----------|-------|----------------|
| < 2 | Not preserved | Module structure lost in reference |
| 2–10 | Moderately preserved | Partial conservation |
| > 10 | Highly preserved | Same module exists in both conditions |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Reference VST matrix | — | Required: baseline/healthy dataset |
| Reference sample metadata | — | Required: metadata for reference |
| Number of permutations | 200 | More → accurate but slower |
| Max genes for test | 2000 | Stratified subsampling limit |
| Random seed | 12345 | For reproducibility |

### How it works internally
1. Finds genes shared between the test and reference datasets
2. Builds expression matrices restricted to shared genes
3. Applies stratified subsampling if `max_genes` is exceeded
4. Runs `modulePreservation()` with the module colors from the main analysis
5. Extracts Zsummary, Zdensity, Zconnectivity statistics

### Outputs
- `table_module_preservation.tsv` — full Z-statistics per module
- `table_preservation_classification.tsv` — simplified Not/Moderate/Highly table
- `plot_zsummary.png` — scatter plot: Zsummary vs log10(module size), colored by level
- `plot_preservation_heatmap.png` — heatmap of all preservation statistics

---

## Cytoscape Import Guide

1. Open Cytoscape
2. **File → Import → Network from File** → select `edges_{module}_kMEfilt.txt`
   - Source: `fromNode`, Target: `toNode`, Interaction type: `weight`
3. **File → Import → Table from File** → select `nodes_{module}_kMEfilt.txt`
   - Select "Node Table" and map `node` to node names
4. Apply a layout (e.g. **Layout → Prefuse Force Directed**)
5. Style nodes by `kME` or `GeneSignificance` columns


---

## Generate Test Data

```bash
Rscript generate_test_data.R
# test-data/test_vst.tabular
# test-data/test_samples.tabular
# test-data/test_up_degs.tabular
# test-data/test_down_degs.tabular
```

For preservation testing, use the same script and treat one dataset as reference.

---

## Troubleshooting

**Auto soft-power fails → uses 6**
→ Inspect the soft-threshold plot; set `--soft_power` manually to where R² ≥ 0.85.

**No modules at correlation threshold**
→ Lower `cor_threshold` (e.g. to 0.4). Biological signal may be weaker.

**GO enrichment: 0 terms found**
→ Verify `gene_id_type` matches your row names exactly (try SYMBOL vs ENSEMBL).
→ Check that the OrgDb package covers your organism.

**KEGG enrichment fails**
→ Requires internet access. Check KEGG organism code at https://www.genome.jp/kegg/catalog/org_list.html

**Preservation very slow**
→ Reduce `pres_n_perms` to 100 or `pres_max_genes` to 1000.

**Few shared genes in preservation**
→ Ensure both datasets use the same gene identifier format in row names.

---

## Citation

If you use this tool, please cite:

> Langfelder P, Horvath S. **WGCNA: an R package for weighted gene co-expression network analysis**. *BMC Bioinformatics* 2008, **9**:559.

> Langfelder P, Luo R, Oldham MC, Horvath S. **Is my network module preserved and reproducible?** *PLoS Computational Biology* 2011, **7**(1):e1001057.

> Yu G et al. **clusterProfiler: a universal enrichment tool for interpreting omics data**. *The Innovation* 2021, **2**(3):100141.
