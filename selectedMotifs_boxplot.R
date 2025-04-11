#!/usr/bin/env Rscript

# --------------------------------------------------------------
# Script: plot_selected_motifs_boxplot.R
# Purpose: Plot selected FEMs (Cancer vs Not cancer)
# Author: Zozan Ismail, 2025
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# -------------------------------
# Define motifs of interest
# -------------------------------
selected_motifs <- c("TAAA", "TATT", "TACA", "TGTA", "CTCC", "CTCG", "GAGG", "AGGG")

# -------------------------------
# Load motif matrix and metadata
# -------------------------------
motif_mat <- read_tsv("fem_output/motif_frequency_matrix.tsv", show_col_types = FALSE) %>%
  pivot_longer(-motif, names_to = "SampleID", values_to = "Frequency") %>%
  filter(motif %in% selected_motifs)

meta <- read_excel("Cohortfil_Zozan_FEM_withsubclasses250408.xlsx") %>%
  rename(SampleID = StudieID, Status = Utfall) %>%
  mutate(
    SampleID = as.character(SampleID),
    Status = case_when(
      str_to_lower(Status) == "cancer" ~ "Cancer",
      str_to_lower(Status) == "not cancer" ~ "Not cancer",
      TRUE ~ Status
    )
  )

# -------------------------------
# Merge motif data with metadata
# -------------------------------
plot_data <- motif_mat %>%
  left_join(meta, by = "SampleID") %>%
  filter(Status %in% c("Cancer", "Not cancer"))

# Set motif order manually (like in figure)
plot_data$motif <- factor(plot_data$motif, levels = selected_motifs)

# -------------------------------
# Generate the plot
# -------------------------------
p <- ggplot(plot_data, aes(x = motif, y = Frequency, fill = Status)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("Cancer" = "#e41a1c", "Not cancer" = "#377eb8")) +
  labs(
    title = "",
    x = "",
    y = "Motif Frequency (%)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("fem_output/selected_motifs_boxplot.pdf", p, width = 8, height = 5)
cat("✅ Selected motif boxplot saved to fem_output/selected_motifs_boxplot.pdf\n")
