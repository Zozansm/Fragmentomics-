#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# Script: 04_build_motif_matrix.R
# Purpose: Normalize 4-mer motif counts to relative frequencies (%),
#          and construct a merged motif frequency matrix across samples.
# Author: [Your Name]
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ------------------------------
# Configuration
# ------------------------------
input_dir <- "fa2"                                  # Folder with *_motif_count.txt files
output_file <- "fem_output/motif_frequency_matrix.tsv"  # Output file

dir.create(dirname(output_file), showWarnings = FALSE)

# ------------------------------
# Generate canonical list of all 256 4-mers
# ------------------------------
valid_4mers <- apply(
  expand.grid(rep(list(c("A", "C", "G", "T")), 4)),
  1,
  paste0,
  collapse = ""
)

# ------------------------------
# Initialize empty list to store vectors
# ------------------------------
motif_matrix_list <- list()

# ------------------------------
# Process each sample file
# ------------------------------
motif_files <- list.files(input_dir, pattern = "_motif_count.txt$", full.names = TRUE)

cat("🔍 Normalizing and merging", length(motif_files), "motif count files...\n")

for (file in motif_files) {
  sample_id <- str_extract(basename(file), "^\\d{4}")
  
  df <- read.table(file, stringsAsFactors = FALSE, col.names = c("count", "motif")) %>%
    filter(str_detect(motif, "^[ACGT]{4}$")) %>%
    mutate(motif = toupper(motif))
  
  # Pad missing motifs
  full <- setNames(rep(0, length(valid_4mers)), valid_4mers)
  full[df$motif] <- df$count
  
  # Normalize to percentages
  total <- sum(full)
  norm_freq <- (full / total) * 100
  
  motif_matrix_list[[sample_id]] <- norm_freq
}

# ------------------------------
# Build matrix
# ------------------------------
motif_matrix <- as.data.frame(do.call(cbind, motif_matrix_list))
motif_matrix <- tibble::rownames_to_column(motif_matrix, "motif")
motif_matrix <- motif_matrix[order(motif_matrix$motif), ]  # Ensure motif row order is consistent
motif_matrix <- motif_matrix[, c("motif", sort(colnames(motif_matrix)[-1]))]  # Sort columns

# ------------------------------
# Save result
# ------------------------------
write.table(
  motif_matrix,
  output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("✅ Motif matrix written to:", output_file, "\n")

