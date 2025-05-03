#!/usr/bin/env Rscript

# -------------------------------------------------------------
# Script: analyze_fragment_size.R
# Purpose: Analyze cfDNA fragment size distributions across samples
# Output: Density plots, median sizes, % <150bp, group comparisons
# Author: Zozan S., 2025
# -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# -------------------------------
# Config
# -------------------------------
bedpe_dir <- "bedpe"
meta_file <- "Cohortfil_Zozan_FEM_withsubclasses250408.xlsx"
output_dir <- "fem_output"
dir.create(output_dir, showWarnings = FALSE)

# -------------------------------
# Load metadata
# -------------------------------
meta <- read_excel(meta_file) %>%
  rename(SampleID = StudieID, Group = Utfall) %>%
  mutate(
    SampleID = as.character(SampleID),
    Group = if_else(tolower(Group) == "cancer", "Cancer", "Not cancer")
  )

# -------------------------------
# Read BEDPE files and extract sizes
# -------------------------------
read_bedpe_sizes <- function(file) {
  sample_id <- str_extract(basename(file), "^\\d{4}")
  df <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  size_col <- ifelse(ncol(df) >= 11, 11, ncol(df))  # fallback
  tibble(SampleID = sample_id, Size = as.numeric(df[[size_col]]))
}

bedpe_files <- list.files(bedpe_dir, pattern = "\\.bedpe$", full.names = TRUE)
all_sizes <- map_dfr(bedpe_files, read_bedpe_sizes)

# Merge with metadata
all_sizes <- left_join(all_sizes, meta, by = "SampleID") %>%
  filter(!is.na(Size) & !is.na(Group))

# -------------------------------
# Plot fragment size density
# -------------------------------
ggplot(all_sizes, aes(x = Size, color = Group)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  coord_cartesian(xlim = c(50, 300)) +
  scale_color_manual(values = c("Cancer" = "#e41a1c", "Not cancer" = "#377eb8")) +
  labs(
    title = "cfDNA Fragment Size Distribution",
    x = "Fragment Size (bp)",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "fragment_size_density.pdf"), width = 8, height = 5)

# -------------------------------
# Summary: Median and % <150 bp
# -------------------------------
summary_df <- all_sizes %>%
  group_by(SampleID, Group) %>%
  summarise(
    median_size = median(Size),
    percent_lt150 = mean(Size < 150) * 100,
    .groups = "drop"
  )

write_tsv(summary_df, file.path(output_dir, "fragment_size_summary.tsv"))

# -------------------------------
# Plot: Boxplot of median sizes
# -------------------------------
p1 <- ggplot(summary_df, aes(x = Group, y = median_size, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Cancer" = "#e41a1c", "Not cancer" = "#377eb8")) +
  labs(y = "Median Fragment Size (bp)", x = NULL) +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "boxplot_median_size.pdf"), p1, width = 6, height = 4)

# -------------------------------
# Plot: Boxplot of % <150 bp
# -------------------------------
p2 <- ggplot(summary_df, aes(x = Group, y = percent_lt150, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Cancer" = "#e41a1c", "Not cancer" = "#377eb8")) +
  labs(y = "Percent Fragments < 150bp", x = NULL) +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "boxplot_percent_lt150bp.pdf"), p2, width = 6, height = 4)

cat("✅ Fragment size analysis complete. Output saved to fem_output/\n")

