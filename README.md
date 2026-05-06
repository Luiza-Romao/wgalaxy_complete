# WGCNA Galaxy Tool

Weighted Gene Co-expression Network Analysis for the Galaxy bioinformatics platform.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Table of contents

- [Overview](#overview)
- [File structure](#file-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input formats](#input-formats)
- [Pipeline steps](#pipeline-steps)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Module preservation / Z-summary](#module-preservation--z-summary)
- [Cytoscape import guide](#cytoscape-import-guide)
- [Generate test data](#generate-test-data)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Authors and funding](#authors-and-funding)
- [License](#license)

---

## Overview

This tool implements a complete WGCNA pipeline on VST-normalized RNA-seq data
inside the Galaxy platform, producing co-expression modules, module-trait
correlation statistics, hub gene tables, and Cytoscape-ready network files.
An optional module preservation analysis (Z-summary) allows cross-condition
or cross-species comparison of co-expression structure.

GO enrichment analysis on WGCNA output modules is available as a separate
companion tool at: https://github.com/Luiza-Romao/go-enrichment-galaxy

---

## File structure

| File | Description |
|---|---|
| `wgcna.xml` | Galaxy tool wrapper — inputs, outputs, UI, help text |
| `macros.xml` | Reusable XML macros (requirements, parameter sections) |
| `wgcna.R` | Main R analysis script |
| `conda_environment.yml` | Reproducible conda environment with pinned dependencies |
| `create_test_subset.R` | R script that generates minimal test datasets for local validation | Utils directory
| `BUGS.md` | Development bug log with root causes and fixes |
| `CHANGELOG.md` | Version history |

---

## Requirements

All Bioconductor packages are pinned to release 3.20 (R 4.4). Mixing packages
from different Bioconductor releases causes dependency resolution failures at
environment build time.

| Package | Version | Source |
|---|---|---|
| r-base | 4.4 | conda-forge |
| r-wgcna | 1.74 | bioconda |
| bioconductor-limma | 3.62.1 | Bioconductor 3.20 |
| bioconductor-impute | 1.80.0 | Bioconductor 3.20 |
| bioconductor-preprocesscore | 1.68.0 | Bioconductor 3.20 |
| r-ggplot2 | 3.5.1 | conda-forge |
| r-ggrepel | 0.9.6 | conda-forge |
| r-reshape2 | 1.4.4 | conda-forge |
| r-pheatmap | 1.0.12 | conda-forge |
| r-igraph | 2.1.1 | conda-forge |
| r-rcolorbrewer | 1.1_3 | conda-forge |
| r-optparse | 1.7.5 | conda-forge |
| r-fastcluster | 1.2.6 | conda-forge |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Luiza-Romao/wgalaxy_complete.git
cd wgalaxy_complete
```

### 2. Create the conda environment

```bash
conda env create -n wgcna_galaxy -f conda_environment.yml
```

### 3. Test locally with Planemo

Before deploying to a production Galaxy instance, validate the tool in a local
simulated Galaxy environment using Planemo:

```bash
conda activate wgcna_galaxy
planemo serve --port 8081
```

Open `http://localhost:8081` in a browser and submit a job using the test data
generated in the step below. This is the recommended validation method for
development and testing.

### 4. Validate syntax and run automated tests

```bash
Rscript -e "parse('wgcna.R'); cat('Syntax OK\n')"
Rscript create_test_subset.R
planemo test wgcna.xml
```

### 5. Deploy to a Galaxy instance

Copy the tool files to the Galaxy tools directory:

```bash
cp -r wgalaxy_complete/ /path/to/galaxy/tools/wgcna_galaxy/
```

Register the tool in `tool_conf.xml`:

```xml
<section id="transcriptomics" name="Transcriptomics">
    <tool file="wgcna_galaxy/wgcna.xml" />
</section>
```

Restart Galaxy after registering the tool.

---

## Generating test data

Follow the steps to generate a VST matrix with 300 genes from *Saccharum spp.* to test the tool.
Open the directory test-data to have more information.

```bash
Rscript create_test_subset.R # Find this file at the utils directory
```

This script writes two files to `test-data/`:

```
test-data/test_vst.tabular
test-data/test_samples.tabular
```

For preservation testing, use the main VST matrix as the query and a second
dataset or a subsampled version as the reference.

---

## Input formats

### VST expression matrix (required)

Tab-separated; first column contains gene IDs; remaining columns contain
sample expression values; first row is the header with sample IDs.

```
GeneID       Sample_A1   Sample_A2   Sample_B1   Sample_B2
AT1G01010    8.23        7.95        5.12        5.44
AT1G01020    6.10        6.35        9.87        10.02
```

### Sample metadata (required)

Tab-separated with a header row. At minimum two columns are required: one for
sample IDs and one for treatment group labels. Column names are configurable
in the tool parameters.

```
SampleID     Treatment
Sample_A1    Infected
Sample_A2    Infected
Sample_B1    Healthy
Sample_B2    Healthy
```

### DEG lists (optional)

Tab-separated, no header required. First column contains gene IDs that must
match the identifiers in the expression matrix. Provide separate files for
up-regulated and down-regulated genes.

```
AT1G01010
AT2G05230
```

### Reference VST matrix (optional, for module preservation)

Same format as the main VST matrix but from a different condition, tissue, or
organism. Required only when the Z-summary module preservation analysis is
enabled.

---

## Pipeline steps

```
INPUT: VST matrix + sample metadata
   |
   ├─  1. Load and align samples
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
   ├─ 13. DEG overlap analysis             [optional]
   └─ 14. Module preservation / Z-summary  [optional]
```

---

## Parameters

### Gene filtering

| Parameter | Default | Description |
|---|---|---|
| Presence fraction | 0.60 | Minimum fraction of samples with non-zero expression for a gene to be retained |
| Variance percentile | 40 | Genes below this variance percentile are removed before network construction |

### Network construction

| Parameter | Default | Description |
|---|---|---|
| Soft power (beta) | 0 (auto) | Adjacency matrix power; 0 triggers automatic selection based on scale-free topology fit |
| Network type | signed | `signed` preserves co-expression direction; `unsigned` treats positive and negative correlations equally |
| Min module size | 50 | Minimum number of genes required to form a co-expression module |
| Merge cut height | 0.25 | Modules whose eigengene Pearson r is equal to or greater than 0.75 are merged |
| Max block size | 5000 | Maximum genes per block in blockwise network construction |
| Threads | 1 | Number of parallel threads for TOM computation |

### Thresholds

| Parameter | Default | Description |
|---|---|---|
| kME threshold | 0.80 | Module membership cutoff for defining core module genes |
| Correlation threshold | 0.60 | Minimum module-trait Pearson r for Cytoscape export and the selected heatmap |
| TOM percentile | 0.95 | Top fraction of TOM edge weights exported to Cytoscape; 0.95 retains the top 5% |

### Module preservation

| Parameter | Default | Description |
|---|---|---|
| Reference VST matrix | — | Required when preservation analysis is enabled |
| Reference sample metadata | — | Required when preservation analysis is enabled |
| Number of permutations | 200 | More permutations increase statistical accuracy but increase runtime proportionally |
| Max genes for test | 2000 | Stratified subsampling limit per module |
| Random seed | 12345 | Ensures reproducibility of permutation results |

---

## Outputs

### Plots (PNG)

| Output | Description |
|---|---|
| Sample clustering | Hierarchical dendrogram with trait color bar for outlier detection |
| Soft-threshold | R-squared and mean connectivity plotted against beta values |
| Gene dendrogram | Hierarchical gene dendrogram with module color bar |
| Full heatmap | Pearson correlation of all modules against all traits |
| MM vs GS | Module membership vs gene significance scatter for the top module |
| Hub genes | Intramodular connectivity vs kME for hub gene identification |
| Selected heatmap | Module-trait heatmap restricted to high-correlation modules |
| DEG barplot | Up- and down-regulated DEG counts per network module (optional) |
| Z-summary scatter | Zsummary vs log10(module size), colored by preservation level (optional) |
| Preservation heatmap | Zsummary, Zdensity, and Zconnectivity statistics per module (optional) |

### Tables (TSV)

| Output | Description |
|---|---|
| Module assignments | Gene-to-module color assignment for all genes in the network |
| GS and MM | Gene significance and kME values for the top correlated module |
| Hub genes | Top 10% intramodular connectivity genes with kME and connectivity scores |
| Selected modules | Summary of modules passing the correlation threshold |
| DEG summary | DEG overlap counts per network module (optional) |
| Preservation full statistics | Zsummary, Zdensity, Zconnectivity per module (optional) |
| Preservation classification | Simplified Not / Moderately / Highly preserved per module (optional) |

### Cytoscape networks (collection)

One edge file and one node file are generated per module passing the
correlation threshold.

- `edges_{module}_kMEfilt.txt` — columns: `fromNode`, `toNode`, `weight` (TOM value)
- `nodes_{module}_kMEfilt.txt` — columns: `node`, `module`, `kME`, `GeneSignificance`, `Connectivity`
- `nodes_{module}_kMEfilt_with_DEG.txt` — as above with an additional `DE_status` column, generated when DEG lists are provided

### Other

- `wgcna_results.RData` — complete R workspace with all objects for downstream custom analysis
- Analysis log — timestamped execution log with step-level diagnostics

---

## Module preservation / Z-summary

### What it does

Uses `WGCNA::modulePreservation()` to test whether modules identified in a
query network also exist as coherent co-expression units in a reference
network. This answers the question: does a module found under one condition
or in one organism also exist under another?

Typical use cases include comparing pathogen-infected versus healthy tissue,
the same tissue at different developmental stages, or conserved modules across
related species such as sugarcane and sorghum.

### Z-summary interpretation

| Zsummary | Level | Interpretation |
|---|---|---|
| Less than 2 | Not preserved | Module co-expression structure is absent in the reference |
| 2 to 10 | Moderately preserved | Partial conservation of co-expression structure |
| Greater than 10 | Highly preserved | The same module exists coherently in both datasets |

### How it works

1. Identifies genes shared between the query and reference expression matrices
2. Restricts shared genes to those that passed the variance and presence filter
   in Step 2 and received a module assignment in the main analysis
3. Applies stratified subsampling per module if the gene count exceeds `max_genes`
4. Runs `modulePreservation()` with the module color vector from the main analysis
5. Filters out the `grey` (unassigned genes) and `gold` (internal null reference)
   pseudo-modules from user-facing results
6. Classifies each module by Zsummary threshold and writes the classification table

---

## Cytoscape import guide

1. Open Cytoscape.
2. Select **File > Import > Network from File** and choose `edges_{module}_kMEfilt.txt`.
   Set Source column to `fromNode`, Target column to `toNode`, and Interaction
   type to `weight`.
3. Select **File > Import > Table from File** and choose `nodes_{module}_kMEfilt.txt`.
   Select Node Table and map the `node` column to node names.
4. Apply a layout such as **Layout > Prefuse Force Directed**.
5. Style nodes using the `kME` or `GeneSignificance` columns for visual encoding.

---

## Troubleshooting

**Conda dependency fails with "seemingly installed but failed to build job environment"**

The most common cause is a version string containing a hyphen in `macros.xml`.
Conda requires underscore notation for CRAN package versions. For example,
`r-rcolorbrewer=1.1-3` must be written as `r-rcolorbrewer=1.1_3`. Check all
version strings in the `<requirements>` block against the conda package index.
See BUG-01 in `BUGS.md` for the full root cause analysis.

**R packages fail to load with version mismatch errors at runtime**

This indicates that packages from different Bioconductor releases are declared
in the same environment. All Bioconductor packages must belong to the same
release. The `conda_environment.yml` in this repository pins all packages to
Bioconductor 3.20 (R 4.4). If you modify the environment, verify version
compatibility at https://bioconductor.org/about/release-announcements/. See
BUG-02 in `BUGS.md` for details.

**Z-summary returns all-NA values without raising an error**

This occurs when the shared gene list includes genes that were removed during
the filtering step in Step 2 and therefore carry no module assignment. Ensure
that the shared gene vector is intersected with `names(module_colors)` before
being passed to `modulePreservation()`. The current version implements this
fix. If the issue persists, verify that both expression matrices use the same
gene identifier format in their row names. See BUG-11 in `BUGS.md`.

**Auto soft-power selection falls back to 6**

Inspect the soft-threshold plot. If the R-squared curve does not reach 0.85
within the tested power range, set `soft_power` manually to the value where
the curve first approaches a plateau. Common values for plant transcriptomes
range from 8 to 14.

**No modules pass the correlation threshold**

Lower `cor_threshold` to 0.4 or 0.5. If the biological signal between
conditions is moderate, high thresholds may exclude all modules. Review the
full module-trait heatmap to assess the correlation distribution before
adjusting this value.

**Preservation analysis is very slow**

Reduce `pres_n_perms` to 100 for exploratory runs. For publication-quality
results, 200 to 500 permutations are recommended. Reducing `pres_max_genes`
to 1000 also shortens runtime when modules are large, at a minor cost to
statistical precision.

---

## Citation

If you use this tool in your research, please cite:

**This tool:**
> Romão, L. O. & Vicentini, R. (2026). WGCNA Galaxy Tool.
> GitHub: https://github.com/Luiza-Romao/wgalaxy_complete
> [![DOI](https://zenodo.org/badge/1198684285.svg)](https://doi.org/10.5281/zenodo.20057034)

**WGCNA R package:**
> Langfelder P, Horvath S. WGCNA: an R package for weighted gene co-expression
> network analysis. *BMC Bioinformatics* 2008, **9**:559.
> https://doi.org/10.1186/1471-2105-9-559

**Module preservation method:**
> Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved
> and reproducible? *PLoS Computational Biology* 2011, **7**(1):e1001057.
> https://doi.org/10.1371/journal.pcbi.1001057

**Galaxy platform:**
> The Galaxy Community. The Galaxy platform for accessible, reproducible, and
> collaborative data analyses: 2024 update. *Nucleic Acids Research* 2024,
> **52**(W1):W83–W94. https://doi.org/10.1093/nar/gkac995

---

## Authors and funding

Developed by **Luiza Oliveira Romão** under FAPESP doctoral fellowship
2025/08740-0, supervised by **Prof. Dr. Renato Vicentini**.

Laboratório de Bioinformática e Biologia de Sistemas  
Instituto de Biologia — Universidade Estadual de Campinas (UNICAMP)

Developed within the **Centro de Melhoramento Molecular de Plantas (CeM²P)**,
funded by FAPESP grant 2022/04006-2.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
