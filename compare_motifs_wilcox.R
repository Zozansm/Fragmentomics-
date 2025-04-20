#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# Script: compare_motifs_wilcox.R
# Purpose: Perform motif-wise statistical comparison between cancer
#          and non-cancer samples using Wilcoxon rank-sum test.
# Output: Table with p-values, adjusted p-values, group medians,
#         and median differences (Δ median)
# Author: Zozan Ismail, 2025
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

# -------------------------------
# 1. Load motif frequency matrix
# -------------------------------
motif_file <- "fem_output/motif_frequency_matrix.tsv"
motif_mat <- read_tsv(motif_file, show_col_types = FALSE) %>%
  column_to_rownames("motif") %>%
  as.data.frame()

# -------------------------------
# 2. Load sample metadata
# -------------------------------
metadata_file <- "Cohortfil_Zozan_FEM_withsubclasses250408.xlsx"
meta <- read_excel(metadata_file) %>%
  rename(SampleID = StudieID, Status = Utfall) %>%
  mutate(
    SampleID = as.character(SampleID),
    Status = case_when(
      str_to_lower(Status) == "cancer" ~ "Cancer",
      str_to_lower(Status) == "not cancer" ~ "Non-cancer",
      TRUE ~ Status
    )
  )

# -------------------------------
# 3. Match and filter samples
# -------------------------------
common_samples <- intersect(meta$SampleID, colnames(motif_mat))
motif_mat <- motif_mat[, common_samples]
meta <- meta %>% filter(SampleID %in% common_samples)

# -------------------------------
# 4. Perform Wilcoxon test + Δ median for each motif
# -------------------------------
results <- lapply(rownames(motif_mat), function(motif) {
  vec <- as.numeric(motif_mat[motif, ])
  group <- meta$Status[match(colnames(motif_mat), meta$SampleID)]
  
  values_cancer <- vec[group == "Cancer"]
  values_control <- vec[group == "Non-cancer"]
  
  if (length(values_cancer) > 0 && length(values_control) > 0) {
    test <- wilcox.test(values_cancer, values_control)
    
    median_cancer <- median(values_cancer, na.rm = TRUE)
    median_control <- median(values_control, na.rm = TRUE)
    delta_median <- median_cancer - median_control
    
    return(data.frame(
      Motif = motif,
      P_value = test$p.value,
      Median_Cancer = round(median_cancer, 2),
      Median_Non_Cancer = round(median_control, 2),
      Delta_Median = round(delta_median, 2)
    ))
  } else {
    return(NULL)
  }
})

# -------------------------------
# 5. Combine and adjust p-values
# -------------------------------
results_df <- bind_rows(results) %>%
  mutate(
    Adj_P_value = round(p.adjust(P_value, method = "fdr"), 3),
    P_value = round(P_value, 3)
  ) %>%
  arrange(Adj_P_value)

# -------------------------------
# 6. Save result table
# -------------------------------
output_file <- "fem_output/wilcox_motif_comparison_with_effectsize.tsv"
write_tsv(results_df, output_file)
cat("✅ Motif-wise comparison saved to:", output_file, "\n")
