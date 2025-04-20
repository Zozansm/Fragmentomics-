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
      tolower(Cancer_Status) == "not cancer" ~ "Non-cancer",
      TRUE ~ as.character(Cancer_Status)
    )
  ) %>%
  filter(SampleID %in% colnames(motif_matrix))

# Reorder motif matrix to match metadata
motif_matrix <- motif_matrix[, metadata$SampleID]

# ------------------------------
# 3. Use all 256 motifs
# ------------------------------
cat("📊 Using all 256 motifs for heatmap.\n")
motif_matrix_subset <- motif_matrix

# ------------------------------
# 4. Prepare annotations
# ------------------------------
annotation_df <- metadata %>%
  select(SampleID, Cancer_Status, Subtype) %>%
  column_to_rownames("SampleID")

# Improved and distinct color palette
ann_colors <- list(
  Cancer_Status = c("Cancer" = "#e41a1c", "Non-cancer" = "#377eb8"),
  Subtype = c(
    "Localized"     = "#1b9e77",  # teal green
    "Metastatic"    = "#d95f02",  # orange
    "Hematopoetic"  = "#7570b3",  # purple
    "Autoimmune"    = "#e7298a",  # pink/red
    "Inflammatory"  = "#66a61e",  # olive green
    "Infectious"    = "#e6ab02",  # mustard yellow
    "No diagnosis"  = "#a6761d",  # brown
    "Other"         = "#666666"   # grey
  )
)

# Custom color gradient for motif values
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
