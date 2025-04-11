#!/usr/bin/env Rscript

# ------------------------------------------------------
# Script: compare_mds_by_group.R
# Purpose: Compare MDS between Cancer and Not Cancer,
#          and color all subtypes (cancer + non-cancer)
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
    Cancer_Status = ifelse(Utfall == "Cancer", "Cancer", "Not Cancer"),
    Subtype = ifelse(Utfall == "Cancer", Disease_Stage, Disease_Stage)  # unified field
  )

# --- Full custom color palette ---
custom_colors <- c(
  # Cancer subtypes
  "Localized"     = "#66c2a5",
  "Metastatic"    = "#984ea3",
  "Hematopoetic"  = "#ff7f00",

  # Non-cancer subtypes
  "Autoimmune"    = "#e41a1c",
  "Inflammatory"  = "#377eb8",
  "Infectious"    = "#a6d854",
  "No diagnosis"  = "#b3e2cd",
  "Other"         = "#fbb4ae"
)

# --- Plot: MDS by group with subtype colors ---
p <- ggplot(merged, aes(x = Cancer_Status, y = MDS)) +
  geom_boxplot(aes(fill = Cancer_Status), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(color = Subtype), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = custom_colors, na.value = "grey50") +
  labs(
    title = "MDS: Cancer vs Not Cancer (All Subtypes Colored)",
    x = "",
    y = "Motif Diversity Score",
    color = "Subtype"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  stat_compare_means(aes(group = Cancer_Status), method = "wilcox.test", label = "p.format")

# --- Save to PDF ---
cat("💾 Saving plot to", output_plot, "\n")
ggsave(output_plot, p, width = 8, height = 5.5)

cat("✅ Done! Check", output_plot, "for results.\n")
