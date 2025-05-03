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
motif_matrix <- read_tsv(
  "fem_output/motif_frequency_matrix.tsv",
  col_types = cols(.default = col_double(), motif = col_character())
)

# Convert numeric column names to character
colnames(motif_matrix) <- as.character(colnames(motif_matrix))
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
    Cancer_Status = str_to_title(str_trim(Cancer_Status)),
    Cancer_Status = case_when(
      Cancer_Status == "Cancer" ~ "Cancer",
      Cancer_Status == "Not Cancer" ~ "Non-cancer",
      TRUE ~ Cancer_Status
    )
  )

# ------------------------------
# 3. Match and align metadata and matrix
# ------------------------------
matched_order <- match(colnames(motif_matrix), metadata$SampleID)

# Handle missing metadata samples
if (any(is.na(matched_order))) {
  missing_samples <- colnames(motif_matrix)[is.na(matched_order)]
  warning("⚠️ The following samples are missing in metadata and will be excluded:\n",
          paste(missing_samples, collapse = ", "))
  
  # Remove unmatched columns from motif matrix
  motif_matrix <- motif_matrix[, !is.na(matched_order)]
  matched_order <- matched_order[!is.na(matched_order)]
}

# Subset and align metadata
metadata <- metadata[matched_order, ]

# ------------------------------
# 4. Prepare annotations
# ------------------------------
annotation_df <- metadata %>%
  select(SampleID, Cancer_Status, Subtype) %>%
  column_to_rownames("SampleID")

ann_colors <- list(
  Cancer_Status = c("Cancer" = "#e41a1c", "Non-cancer" = "#377eb8"),
  Subtype = c(
    "Localized"     = "#1b9e77",
    "Metastatic"    = "#d95f02",
    "Hematopoetic"  = "#7570b3",
    "Autoimmune"    = "#e7298a",
    "Inflammatory"  = "#66a61e",
    "Infectious"    = "#e6ab02",
    "No diagnosis"  = "#a6761d",
    "Other"         = "#666666"
  )
)

color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# ------------------------------
# 5. Generate the heatmap
# ------------------------------
pheatmap(
  as.matrix(motif_matrix),
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = color_palette,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",  # Z-score normalization by motif
  clustering_method = "ward.D2",
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  fontsize = 10,
  border_color = NA,
  width = 14,
  height = 20,
  filename = "fem_output/heatmap_fem.pdf"
)

cat("✅ FEM heatmap saved to fem_output/heatmap_fem.pdf\n")
