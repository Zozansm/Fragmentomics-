#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script: 05_heatmap_fem.R
# Purpose: Generate a heatmap of normalized FEM motif frequencies
#          with annotations for cancer status and subtype.
# Author: Zozan S., 2025
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(pheatmap)
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(RColorBrewer)
})

# ------------------------------
# 1. Load motif frequency matrix
# ------------------------------
motif_matrix <- read_tsv("fem_output/motif_frequency_matrix.tsv")
rownames(motif_matrix) <- motif_matrix$motif
motif_matrix <- motif_matrix[, -1]  # Remove 'motif' column

# ------------------------------
# 2. Load and clean metadata
# ------------------------------
meta_raw <- read_excel("Cohortfil_Zozan_FEM_withsubclasses250408.xlsx")

metadata <- meta_raw %>%
  rename(
    SampleID = StudieID,
    Cancer_Status = Utfall,
    Subtype = Disease_Stage
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    Cancer_Status = case_when(
      tolower(Cancer_Status) == "cancer" ~ "Cancer",
      tolower(Cancer_Status) == "not cancer" ~ "Control",
      TRUE ~ as.character(Cancer_Status)
    )
  ) %>%
  filter(SampleID %in% colnames(motif_matrix))

# Reorder motif matrix to match metadata
motif_matrix <- motif_matrix[, metadata$SampleID]

# ------------------------------
# 3. Filter top 100 variable motifs with fallback
# ------------------------------
top_n <- 100
motif_var <- apply(motif_matrix, 1, var)

cat("📊 Motif variance summary:\n")
print(summary(motif_var))

top_motifs <- names(sort(motif_var, decreasing = TRUE))[1:top_n]

if (length(top_motifs) == 0 || any(is.na(top_motifs))) {
  warning("⚠️ No top motifs found — using full motif matrix.")
  motif_matrix_subset <- motif_matrix
} else {
  cat("✅ Top motifs selected:\n")
  print(head(top_motifs))
  motif_matrix_subset <- motif_matrix[top_motifs, ]
}

# ------------------------------
# 4. Prepare annotations
# ------------------------------
annotation_df <- metadata %>%
  select(SampleID, Cancer_Status, Subtype) %>%
  column_to_rownames("SampleID")

ann_colors <- list(
  Cancer_Status = c("Cancer" = "#e41a1c", "Control" = "#377eb8"),
  Subtype = c(
    "Localized"     = "#4daf4a",
    "Metastatic"    = "#984ea3",
    "Hematopoetic"  = "#ff7f00",
    "Autoimmune"    = "#a6cee3",
    "Inflammatory"  = "#1f78b4",
    "Infectious"    = "#b2df8a",
    "No diagnosis"  = "#33a02c",
    "Other"         = "#fb9a99"
  )
)

# Custom color palette
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# ------------------------------
# 5. Plot and save heatmap
# ------------------------------
pheatmap(
  as.matrix(motif_matrix_subset),
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = color_palette,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",  # Normalize each motif across samples
  clustering_method = "ward.D2",
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  fontsize = 10,
  border_color = NA,
  filename = "fem_output/heatmap_fem.pdf",
  width = 12,
  height = 8
)

cat("✅ FEM heatmap saved to fem_output/heatmap_fem.pdf\n")
