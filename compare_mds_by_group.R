#!/usr/bin/env Rscript

# ------------------------------------------------------
# Script: compare_mds_by_group.R
# Purpose: Compare MDS between Cancer and Not Cancer,
#          color all subtypes, and report median/IQR values
# Author: Zozan Ismail, 2025
# ------------------------------------------------------

# Load libraries
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

# --- File paths ---
metadata_file <- "Cohortfil_Zozan_FEM_withsubclasses250408.xlsx"
mds_file <- "fem_output/mds_summary.tsv"
output_plot <- "fem_output/MDS_all_subtypes_colored.pdf"

# --- Load data ---
cat("📂 Reading metadata and MDS file...\n")
meta <- read_excel(metadata_file)
mds <- read.table(mds_file, header = TRUE, sep = "\t")

# --- Merge metadata with MDS scores ---
merged <- meta %>%
  mutate(SampleID = as.character(StudieID)) %>%
  inner_join(mds %>% mutate(SampleID = as.character(SampleID)), by = "SampleID")

# --- Define cancer group and assign subtypes ---
merged <- merged %>%
  mutate(
    Cancer_Status = ifelse(tolower(Utfall) == "cancer", "Cancer", "Not Cancer"),
    Subtype = Disease_Stage
  )

# --- Compute median and IQR for each group ---
summary_stats <- merged %>%
  group_by(Cancer_Status) %>%
  summarise(
    Median = round(median(MDS), 5),
    IQR_Lower = round(quantile(MDS, 0.25), 5),
    IQR_Upper = round(quantile(MDS, 0.75), 5),
    .groups = "drop"
  )

cat("📊 Median and IQR of MDS by group:\n")
print(summary_stats)

# --- Improved subtype color palette ---
custom_colors <- c(
  "Localized"     = "#1b9e77",  # teal
  "Metastatic"    = "#d95f02",  # orange
  "Hematopoetic"  = "#7570b3",  # purple
  "Autoimmune"    = "#e7298a",  # pink
  "Inflammatory"  = "#66a61e",  # green
  "Infectious"    = "#e6ab02",  # yellow-brown
  "No diagnosis"  = "#a6761d",  # brown
  "Other"         = "#666666"   # dark grey
)

# --- Plot: MDS by group with subtype colors ---
p <- ggplot(merged, aes(x = Cancer_Status, y = MDS)) +
  geom_boxplot(aes(fill = Cancer_Status), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(color = Subtype), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = custom_colors, na.value = "grey50") +
  labs(
    title = "",
    x = "",
    y = "Motif Diversity Score (MDS)",
    color = "Subtype"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  stat_compare_means(aes(group = Cancer_Status), method = "wilcox.test", label = "p.format")

# --- Save to PDF ---
cat("💾 Saving plot to", output_plot, "\n")
ggsave(output_plot, p, width = 8, height = 5.5)

cat("✅ Done! Check", output_plot, "for results.\n")
