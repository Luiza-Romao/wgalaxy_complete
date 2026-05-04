# WGCNA Galaxy Tool — Bug Log

Reconstructed from development sessions (April–May 2026).  
Each entry follows the format: **symptom → root cause → fix applied → status**.

---

## Environment & Installation

### BUG-01 · Conda dependency fails to build job environment
- **Symptom:** `Conda dependency seemingly installed but failed to build job environment` immediately after submitting any job.
- **Root cause:** `r-rcolorbrewer` was declared with version `1.1-3` (CRAN notation using a hyphen). Conda requires underscore notation: `1.1_3`. The Galaxy dependency resolver located the package in the index but failed to install it due to the invalid version string.
- **Fix:** Changed version string to `1.1_3` in `macros.xml` `<requirements>` block.
- **Status:** Closed

---

### BUG-02 · Bioconductor package version conflicts
- **Symptom:** Conda environment resolves but R packages fail to load with version mismatch errors at runtime.
- **Root cause:** Packages from different Bioconductor releases (3.18 and 3.19, corresponding to R 4.3 and R 4.4 respectively) were mixed in the same `<requirements>` block. Bioconductor packages within a release are co-versioned and not interchangeable across releases.
- **Fix:** Pinned all Bioconductor packages to a single release (3.19 / R 4.4) and verified versions against the Bioconductor package matrix before declaring them in `macros.xml`.
- **Status:** Closed

---

## XML & Galaxy Interface

### BUG-03 · "Error occurred while building command line for tool"
- **Symptom:** Job fails within seconds with the error `Error occurred while building command line for tool 'wgcna_analysis'`. No R code is executed.
- **Root cause:** Parameter names in the `<command>` block of `wgcna.xml` did not match the option names parsed by `optparse` in `wgcna.R`. Galaxy builds the command line by substituting `$param_name` from the XML; when the name is wrong, it injects an empty string or fails to construct a valid argument.
- **Fix:** Audited all `--option-name` flags in the `<command>` block against the `make_option()` calls in the R script and corrected every mismatch.
- **Status:** Closed

---

### BUG-04 · Help section renders as escaped HTML
- **Symptom:** The Galaxy tool help page displays raw HTML tags as literal text (e.g., `<p class="infomark"><strong>What it does</strong></p>`).
- **Root cause (a):** `.. code-block::` directive was used without a language specifier. Galaxy uses `docutils` without Sphinx/Pygments, which does not support this directive. When the RST parser aborts on the directive, it falls back to rendering the entire section as escaped text.
- **Root cause (b):** Emoji characters (❌ ⚠️ ✅) inside RST grid tables corrupt the column width alignment. A single misaligned byte causes the docutils table parser to fail and propagate the error to subsequent blocks.
- **Root cause (c):** `**bold** *(NEW)*` constructs inside table cells created unclosed inline markup nodes, which cascaded into broken HTML for everything following.
- **Fix:** Replaced `.. code-block::` with standard RST `::` paragraph ending. Removed all emoji from tables. Replaced all RST grid tables with RST bullet lists, which are immune to alignment errors. Removed `*(NEW)*` and similar AI-generated markers.
- **Status:** Closed

---

## R Script — Core Pipeline

### BUG-05 · Syntax error: `unexpected 'else'` on `soft_power` assignment
- **Symptom:** R script fails immediately with `Error: unexpected 'else' in "else"`. Job aborts before any analysis runs.
- **Root cause:** R closes an `if` expression at end-of-line when no continuation is visible. The multi-line `if / else if / else` block for `soft_power` was written without enclosing braces, so R parsed the first line as a complete expression and treated the subsequent `else` as a stray token.
- **Fix:** Wrapped all branches in explicit curly braces:
  ```r
  soft_power <- if (opt$soft_power > 0) {
    opt$soft_power
  } else if (!is.na(sft$powerEstimate)) {
    sft$powerEstimate
  } else {
    6
  }
  ```
- **Status:** Closed

---

### BUG-06 · `--out_log` argument not declared in optparse
- **Symptom:** R script fails with `Error in getopt(spec, opt = args) : Unknown option: --out_log`.
- **Root cause:** The `--out_log` parameter was added to the XML `<command>` block but never registered in the `make_option()` list in the R script.
- **Fix:** Added `make_option("--out_log", type="character", default="analysis_log.txt")` to the options list.
- **Status:** Closed

---

## R Script — Step 7 (Module–Trait Correlation)

### BUG-07 · MEgrey module included in filtered correlation heatmap
- **Symptom:** The heatmap of modules passing the correlation threshold included the grey module, which is the WGCNA catch-all for unassigned genes and has no biological meaning.
- **Root cause:** The filter `pres_df[pres_df$Module != "grey", ]` was not applied before generating the Step 7 heatmap.
- **Fix:** Added `module_colors` filter to exclude `MEgrey` before computing the filtered module list for Step 7 plots.
- **Status:** Closed

---

### BUG-08 · "ME" prefix retained in module labels on plots
- **Symptom:** Module names in heatmaps and tables showed as `MEturquoise`, `MEblue`, etc. instead of `turquoise`, `blue`.
- **Root cause:** `colnames(MEs)` returns names with the `ME` prefix (eigengene convention). The `gsub("^ME", "", ...)` strip was missing from the plotting and table-writing blocks in Step 7.
- **Fix:** Applied `gsub("^ME", "", module_name)` to all labels before plotting and writing to TSV.
- **Status:** Closed

---

### BUG-09 · `cex.text` missing from filtered module heatmap
- **Symptom:** Text labels overflowed or were cut off in the filtered-module `labeledHeatmap`, while the full heatmap rendered correctly.
- **Root cause:** `cex.text=0.6` was set in the full heatmap call but omitted in the filtered-module version, causing inconsistent text sizing and occasional overflow when module names were long.
- **Fix:** Added `cex.text=0.6` to the filtered `labeledHeatmap()` call to match the full heatmap.
- **Status:** Closed

---

## R Script — Step 14/15 (Module Preservation / Z-summary)

### BUG-10 · Z-summary outputs not delivered; entire Step 14 job fails
- **Symptom:** Step 14 job is marked as failed. None of the Z-summary outputs (scatter plot, heatmap, classification table) appear in the Galaxy history.
- **Root cause:** `out_plot_pres_heatmap` and `out_table_pres_class` were declared as Galaxy `<data>` outputs in `wgcna.xml` but were never generated by the R script. Galaxy waits for all declared outputs; when they are absent, it marks the job as failed and suppresses all other outputs from that step, including the Z-summary scatter plot that was actually being generated correctly.
- **Fix:** Implemented the missing heatmap and classification table generation in the R script. Added `placeholder_png()` calls for all output paths in the `length(shared) < 50` early-exit branch so Galaxy always receives a file.
- **Status:** Closed

---

### BUG-11 · `modulePreservation()` returns all-NA Z-scores silently
- **Symptom:** Z-summary scatter plot and table contain only NA values for all modules. No error is thrown.
- **Root cause:** `shared` genes were computed as `intersect(colnames(query_expr_mat), rownames(ref_expr))` but not intersected with `names(module_colors)`. Genes present in the raw VST matrices but removed during Step 2 filtering (presence/variance thresholds) were included in `cols_sub` with no module assignment, producing NA entries. `modulePreservation()` silently propagates NAs through all statistics when the color vector contains NAs.
- **Fix:**
  ```r
  shared <- intersect(colnames(query_expr_mat), rownames(ref_expr))
  shared <- intersect(shared, names(module_colors))  # restrict to network genes
  ```
- **Status:** Closed

---

### BUG-12 · "gold" pseudo-module included in preservation results
- **Symptom:** Z-summary table contains a row for a module named `gold` which does not correspond to any biological module. Its statistics are inflated and misleading.
- **Root cause:** `WGCNA::modulePreservation()` internally adds a `gold` pseudo-module (a random gene sample used as a null reference for the permutation test). This module must be excluded from user-facing results. The filter only excluded `grey` but not `gold`.
- **Fix:** Changed the filter from:
  ```r
  pres_df[pres_df$Module != "grey", ]
  ```
  to:
  ```r
  pres_df[!pres_df$Module %in% c("grey", "gold"), ]
  ```
- **Status:** Closed

---

## Summary Table

| ID | Component | Severity | Status |
|---|---|---|---|
| BUG-01 | Conda / macros.xml | Critical | Closed |
| BUG-02 | Conda / macros.xml | Critical | Closed |
| BUG-03 | wgcna.xml command block | Critical | Closed |
| BUG-04 | wgcna.xml help RST | Moderate | Closed |
| BUG-05 | wgcna.R parser | Critical | Closed |
| BUG-06 | wgcna.R optparse | Critical | Closed |
| BUG-07 | wgcna.R Step 7 | Moderate | Closed |
| BUG-08 | wgcna.R Step 7 | Minor | Closed |
| BUG-09 | wgcna.R Step 7 | Minor | Closed |
| BUG-10 | wgcna.R Step 14 | Critical | Closed |
| BUG-11 | wgcna.R Step 15 | Critical | Closed |
| BUG-12 | wgcna.R Step 15 | Moderate | Closed |

All 12 bugs closed as of May 2026. Tool validated end-to-end on sugarcane RD/SD dataset (SRA accessions SRP309574).
