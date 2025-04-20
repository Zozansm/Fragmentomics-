
# Fragmentomics: Analysis of Fragment End Motifs in cfDNA

Welcome to the official repository for the **Fragmentomics** Bachelor’s thesis project. This study explores the diagnostic potential of **fragment end motifs (FEMs)** and **fragment size distributions** in circulating cell-free DNA (cfDNA) using data from patients with cancer and no cancer. The repository includes the full analysis pipeline implemented in **Bash** and **R**, with scripts for filtering sequencing data, computing motif frequencies, and visualizing results.

---

## 📁 Repository Contents

| File/Folder | Description |
|-------------|-------------|
| `01_filter_bam.sh` | Filters BAM files to retain high-quality, properly paired reads from autosomes only. |
| `02_extract_motifs.sh` | Extracts 5′ end 4-mer sequences from cfDNA fragments. |
| `03_calculate_mds.R` | Calculates Motif Diversity Score (MDS) per sample using normalized Shannon entropy. |
| `04_build_motif_matrix.R` | Constructs a motif frequency matrix across samples for downstream analysis. |
| `05_heatmap_fem.R` | Generates a heatmap to visualize global FEM patterns across all samples. |
| `compare_mds_by_group.R` | Compares MDS scores between cancer and non-cancer groups with Wilcoxon test. |
| `compare_motifs_wilcox.R` | Performs motif-wise Wilcoxon testing and ranks significant FEMs. |
| `selectedMotifs_boxplot.R` | Plots boxplots of selected motifs to visualize frequency shifts by group. |
| `fem_output/` | Folder containing result files like heatmaps, summary stats, and motif tables. |

---

## 📖 Project Overview

### Introduction

Fragmentation patterns in cfDNA reflect biological processes like apoptosis and chromatin structure, and may carry diagnostic information. This study investigates whether **5′ end 4-mer motifs** and **fragment size** characteristics differ between cancer and control patients.

### Objectives

- Identify statistically significant differences in FEMs between cancer and control samples.
- Explore global diversity of FEMs using MDS.
- Visualize motif enrichment via heatmaps and boxplots.
- Evaluate motifs for biomarker potential.

---

## 🧪 Methods Summary

1. **BAM Filtering** (`01_filter_bam.sh`)
   - Removes low-quality reads.
   - Restricts to autosomal chromosomes.
   - Prepares data for FEM extraction.

2. **Motif Extraction** (`02_extract_motifs.sh`)
   - Extracts 4-mer sequences from the 5′ ends of cfDNA fragments.
   - Uses reference genome via `bedtools getfasta`.

3. **Motif Frequency Analysis**
   - Build motif frequency matrix (`04_build_motif_matrix.R`).
   - Visualize FEMs in a clustered heatmap (`05_heatmap_fem.R`).
   - Test FEM-level differences (`compare_motifs_wilcox.R`).
   - Plot top motifs (`selectedMotifs_boxplot.R`).

4. **Statistical Testing**
   - Wilcoxon rank-sum tests to compare cancer vs. non-cancer.
   - Adjust p-values using Benjamini-Hochberg correction.

---

## 🧬 Results Summary

- Certain motifs (e.g., `TGTA`, `TATT`, `TAAA`) showed differential abundance between cancer and non-cancer groups.
- Global MDS did not significantly differentiate cancer status.
- Stratified heatmaps revealed clusters enriched in hematopoietic and metastatic cancers.

See figures and tables in the `plots/` or `fem_output/` directories for full visual results.

---

```r
install.packages(c("ggplot2", "pheatmap", "readr", "readxl", "dplyr", "tibble", "stringr", "ggpubr", "RColorBrewer"))
## 🚀 Usage

Run the pipeline step by step or selectively:

```bash
# 1. Filter BAM files
bash 01_filter_bam.sh

# 2. Extract 4-mer FEMs
bash 02_extract_motifs.sh

# 3. Compute MDS
Rscript 03_calculate_mds.R

# 4. Build motif matrix
Rscript 04_build_motif_matrix.R

# 5. Visualize results
Rscript 05_heatmap_fem.R
Rscript compare_motifs_wilcox.R
