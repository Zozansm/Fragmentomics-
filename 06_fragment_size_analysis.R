#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script: 06_fragment_size_analysis.R
# Purpose: Analyze cfDNA fragment size distribution between cancer and non-cancer groups
#          Includes summary statistics, statistical testing, and visualization.
# Author: Zozan Ismail, 2025
# -----------------------------------------------------------

# Load libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
})

# ------------------------------
# Configuration
# ------------------------------
bedpe_dir <- "bedpe"         # Folder with *.bedpe files
metadata_file <- "Cohortfil_Zozan_FEM_withsubclasses250408.xlsx"
output_dir <- "fem_output"
summary_file <- file.path(output_dir, "fragment_size_summary.tsv")

# Create output folder if missing
dir.create(output_dir, showWarnings = FALSE)

# ------------------------------
# 1. Read and summarize fragment sizes
# ------------------------------

cat("📂 Reading BEDPE files and summarizing fragment sizes...\n")

bedpe_files <- list.files(bedpe_dir, pattern = ".bedpe$", full.names = TRUE)

fragment_summary_list <- lapply(bedpe_files, function(file) {
  sample_id <- str_extract(basename(file), "^\\d{4}")
  df <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)

  sizes <- df$X11  # Assuming fragment size is in column 11

  tibble(
    SampleID = sample_id,
    Mean_Size = mean(sizes, na.rm = TRUE),
    Median_Size = median(sizes, na.rm = TRUE),
    IQR_Size = IQR(sizes, na.rm = TRUE),
    Percent_Short = mean(sizes < 150, na.rm = TRUE) * 100  # % fragments <150bp
  )
})

fragment_summary <- bind_rows(fragment_summary_list)

# Save summary table
write_tsv(fragment_summary, summary_file)

cat("✅ Fragment size summary written to:", summary_file, "\n")

# ------------------------------
# 2. Merge with metadata
# ------------------------------
meta <- read_excel(metadata_file) %>%
  rename(SampleID = StudieID, Cancer_Status = Utfall) %>%
  mutate(
    SampleID = as.character(SampleID),
    Cancer_Status = case_when(
      tolower(Cancer_Status) == "cancer" ~ "Cancer",
      tolower(Cancer_Status) == "not cancer" ~ "Non-cancer",
      TRUE ~ Cancer_Status
    )
  )

merged_data <- fragment_summary %>%
  inner_join(meta, by = "SampleID")

# ------------------------------
# 3. Statistical Tests
# ------------------------------

cat("\n🔬 Statistical tests between Cancer and Non-cancer groups:\n")

# Wilcoxon test for Median Size
wilcox_median <- wilcox.test(Median_Size ~ Cancer_Status, data = merged_data)

cat("\n- Wilcoxon test (Median fragment size):\n")
print(wilcox_median)

# Wilcoxon test for %<150bp
wilcox_short <- wilcox.test(Percent_Short ~ Cancer_Status, data = merged_data)

cat("\n- Wilcoxon test (Percent short fragments <150bp):\n")
print(wilcox_short)

# ------------------------------
# 4. Visualization
# ------------------------------

# Kernel Density Plot
p_density <- ggplot(merged_data, aes(x = Median_Size, fill = Cancer_Status, color = Cancer_Status)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 166, linetype = "dashed", color = "gray") +
  labs(
    title = "cfDNA Fragment Size Density",
    x = "Median Fragment Size (bp)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "fragment_size_density.pdf"), p_density, width = 8, height = 5)

# Violin Plot: Median Sizes
p_violin <- ggplot(merged_data, aes(x = Cancer_Status, y = Median_Size, fill = Cancer_Status)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Median Fragment Size by Group", x = "Group", y = "Median Size (bp)")

ggsave(file.path(output_dir, "fragment_size_violin_median.pdf"), p_violin, width = 6, height = 5)

# Violin Plot: Percent short fragments
p_short <- ggplot(merged_data, aes(x = Cancer_Status, y = Percent_Short, fill = Cancer_Status)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Percentage of Fragments <150bp by Group", x = "Group", y = "% Short Fragments")

ggsave(file.path(output_dir, "fragment_size_violin_percent_short.pdf"), p_short, width = 6, height = 5)

cat("\n✅ Plots saved to:", output_dir, "\n")

cat("\n🎯 Fragment size analysis completed.\n")
