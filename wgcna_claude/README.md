# WGCNA Galaxy Tool

**Weighted Gene Co-expression Network Analysis** integrated into the Galaxy bioinformatics platform.

---

## Files

| File | Description |
|------|-------------|
| `wgcna.xml` | Galaxy tool wrapper (inputs, outputs, UI) |
| `macros.xml` | Reusable XML macros (requirements, parameters) |
| `wgcna.R` | Main analysis R script |
| `conda_environment.yml` | Conda environment for reproducibility |
| `generate_test_data.R` | Script to create minimal test datasets |

---

## Installation in Galaxy

### Option 1 — Via Galaxy ToolShed (recommended for admins)
Upload the tool to a local or public ToolShed and install from Galaxy's Admin panel.

### Option 2 — Manual installation

1. Copy all files to a directory inside your Galaxy tool root, e.g.:
   ```bash
   cp -r wgcna_galaxy/ /galaxy/tools/wgcna_galaxy/
   ```

2. Register the tool in `tool_conf.xml` (or `local_tool_conf.xml`):
   ```xml
   <section id="transcriptomics" name="Transcriptomics">
       <tool file="wgcna_galaxy/wgcna.xml" />
   </section>
   ```

3. Create the conda environment (Galaxy will do this automatically if
   `conda_auto_install = True` in `galaxy.yml`):
   ```bash
   conda env create -f wgcna_galaxy/conda_environment.yml
   ```

4. Restart Galaxy and the tool will appear in the tool panel.

---

## Required Input Formats

### VST Expression Matrix (`tabular`)
```
GeneID      Sample_A_1  Sample_A_2  Sample_B_1  Sample_B_2
Gene_001    8.23        7.95        5.12        5.44
Gene_002    6.10        6.35        9.87        10.02
...
```
- Tab-separated
- First column: gene identifiers (row names)
- Remaining columns: one per sample (header = sample IDs)

### Sample Metadata (`tabular`)
```
SampleID    Treatment
Sample_A_1  Induced
Sample_A_2  Induced
Sample_B_1  Control
Sample_B_2  Control
```
- Tab-separated
- Header required
- Column names are configurable in the tool

### DEG Gene Lists (`tabular`, optional)
```
Gene_001
Gene_015
Gene_089
```
- Tab-separated, no header required
- First column = gene IDs (must match expression matrix)
- Provide separately for up- and down-regulated genes

---

## Pipeline Steps

```
INPUT: VST matrix + sample metadata
   │
   ├─ 1. Load & align data
   ├─ 2. Filter: presence (≥60% samples) + variance (40th percentile)
   ├─ 3. Build trait matrix (binary: case vs control)
   ├─ 4. Sample clustering (outlier detection)
   ├─ 5. Soft-threshold selection (scale-free topology)
   ├─ 6. Blockwise network construction (signed, TOM-based)
   ├─ 7. Module eigengenes + module-trait Pearson correlation
   ├─ 8. Gene significance (GS) + module membership (kME) — top module
   ├─ 9. Hub gene identification (top 10% intramodular connectivity)
   ├─ 10. kME filtering (|kME| > 0.8)
   ├─ 11. Cytoscape export — modules with |r| ≥ 0.6
   └─ 12. DEG overlap analysis (optional)
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Presence fraction | 0.60 | Min fraction of samples with non-zero expression |
| Variance percentile | 40 | Remove genes below this variance percentile |
| Soft power (β) | 0 (auto) | Network adjacency power; 0 = auto-detect |
| Network type | signed | `signed` preserves co-expression direction |
| Min module size | 50 | Minimum genes per co-expression module |
| Merge cut height | 0.25 | Modules with eigengene r ≥ 0.75 are merged |
| kME threshold | 0.80 | Module membership filter for core genes |
| Correlation threshold | 0.60 | |r| cutoff for Cytoscape module export |
| TOM percentile | 0.95 | Top 5% TOM edges exported to Cytoscape |

---

## Outputs

### Plots (PNG)
| Output | Description |
|--------|-------------|
| Sample clustering | Hierarchical clustering + trait color bar (outlier check) |
| Soft-threshold | R² and mean connectivity vs β |
| Gene dendrogram | Cluster dendrogram with module color bar |
| Full module-trait heatmap | All modules × all traits |
| MM vs GS | Module membership vs gene significance (top module) |
| Hub gene plot | Connectivity vs kME, hubs highlighted in red |
| Selected modules heatmap | High-correlation modules only |
| DEG barplot *(optional)* | Up/down DEGs per network module |

### Tables (TSV)
| Output | Description |
|--------|-------------|
| Module assignments | Gene → module color for all genes |
| GS/MM table | Gene significance + module membership (top module) |
| Hub genes | Hub gene list with connectivity and kME |
| Selected modules | Summary of high-correlation modules |
| DEG network summary *(optional)* | DEG counts per network module |

### Cytoscape Networks (collection of tabular files)
- `edges_{module}_kMEfilt.txt` — Edge table: `fromNode`, `toNode`, `weight` (TOM)
- `nodes_{module}_kMEfilt.txt` — Node table: `node`, `module`, `kME`, `GeneSignificance`, `Connectivity`
- `nodes_{module}_kMEfilt_with_DEG.txt` *(optional)* — Nodes with `DE_status` column

### Other
- `wgcna_results.RData` — All R objects for downstream custom analysis
- Analysis log — stdout/stderr

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

## Citation

If you use this tool, please cite:

> Langfelder P, Horvath S. **WGCNA: an R package for weighted gene co-expression network analysis**.
> *BMC Bioinformatics* 2008, **9**:559. https://doi.org/10.1186/1471-2105-9-559

> Zhang B, Horvath S. **A general framework for weighted gene co-expression network analysis**.
> *Statistical Applications in Genetics and Molecular Biology* 2005, **4**(1).

---

## Generating Test Data

```bash
Rscript generate_test_data.R
# Creates: test-data/test_vst.tabular
#          test-data/test_samples.tabular
#          test-data/test_up_degs.tabular
#          test-data/test_down_degs.tabular
```

---

## Troubleshooting

**Auto soft-power estimation failed**
→ Inspect the soft-threshold plot. Manually set `--soft_power` to the value where
the scale-free fit index (R²) crosses 0.85.

**No modules passed correlation threshold**
→ Lower `cor_threshold` (e.g. to 0.4). This can happen with few samples or weak biology.

**Out of memory during blockwiseModules**
→ Reduce `max_block_size` (e.g. to 3000). This splits large gene sets into smaller blocks.

**Sample IDs mismatch**
→ Ensure column headers in the VST matrix exactly match the `SampleID` column in metadata.
