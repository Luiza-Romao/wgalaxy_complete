# test-data

This directory contains the test dataset used by `planemo test wgcna.xml`.

---

## Files

| File | Description |
|---|---|
| `test_vst.tabular` | VST-normalized expression matrix — 300 genes x 14 samples |
| `test_samples.tabular` | Sample metadata — 14 samples, 2 treatments (RD and SD) |

---

## Data origin

The test files are a subset of a publicly available sugarcane RNA-seq dataset
deposited at NCBI SRA under accession **SRP309574**.

The 14 samples represent two photoperiodic conditions applied to
*Saccharum officinarum* hybrid plants:

- **RD** (Regular Day — 12h light / 12h dark) — 7 replicates
- **SD** (Short Day — 8h light / 16h dark) — 7 replicates

Expression was quantified with Salmon against the sugarcane reference
transcriptome and normalized using DESeq2 variance-stabilizing transformation
(VST). The 300 genes in the test subset were selected to include 150
high-variance genes (enriched for transcripts responding to photoperiod
treatment) and 150 randomly sampled genes, ensuring that the WGCNA pipeline
produces detectable co-expression modules and module-trait correlations during
automated testing.

---

## Regenerating the subset

If you have access to the full VST matrix from the pipeline, run:

```bash
conda activate wgcna_galaxy
Rscript create_test_subset.R \
  --vst_matrix   path/to/full_vst_matrix.tsv \
  --sample_info  path/to/sample_info.tsv \
  --out_dir      test-data
```

The script selects genes deterministically given the fixed random seed (42)
and produces identical output on any machine with the same input data.
