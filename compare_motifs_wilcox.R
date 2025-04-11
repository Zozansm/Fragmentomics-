#!/usr/bin/env Rscript

# --------------------------------------------------------------
# Script: compare_motifs_wilcox.R
# Purpose: Compare 4-mer motif frequencies between Cancer and Control
# Author: Zozan Ismail, 2025
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tibble)
  library(stringr)
})

# -------------------------------
# Load motif frequency matrix
# -------------------------------
motif_mat <- read_tsv("fem_output/motif_frequency_matrix.tsv", show_col_types = FALSE) %>%
  column_to_rownames("motif") %>%
  as.data.frame()

# -------------------------------
# Load metadata
# -------------------------------
meta <- read_excel("Cohortfil_Zozan_FEM_withsubclasses250408.xlsx") %>%
  rename(SampleID = StudieID, Status = Utfall) %>%
  mutate(
    SampleID = as.character(SampleID),
    Status = case_when(
      str_to_lower(Status) == "cancer" ~ "Cancer",
      str_to_lower(Status) == "not cancer" ~ "Control",
      TRUE ~ Status
    )
  )

# Keep only samples present in both data sets
common_samples <- intersect(meta$SampleID, colnames(motif_mat))
motif_mat <- motif_mat[, common_samples]
meta <- meta %>% filter(SampleID %in% common_samples)

# -------------------------------
# Wilcoxon test for each motif
# -------------------------------
results <- lapply(rownames(motif_mat), function(motif) {
  vec <- as.numeric(motif_mat[motif, ])
  group <- meta$Status[match(colnames(motif_mat), meta$SampleID)]

  if (length(unique(group)) == 2) {
    test <- wilcox.test(vec[group == "Cancer"], vec[group == "Control"])
    median_diff <- median(vec[group == "Cancer"]) - median(vec[group == "Control"])
    return(data.frame(
      motif = motif,
      p_value = test$p.value,
      median_diff = median_diff
    ))
  } else {
    return(NULL)
  }
})

# Combine results and adjust p-values
df <- bind_rows(results) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

# -------------------------------
# Output
# -------------------------------
write_tsv(df, "fem_output/wilcox_motif_comparison.tsv")
cat("✅ Wilcoxon test results saved to fem_output/wilcox_motif_comparison.tsv\n")
