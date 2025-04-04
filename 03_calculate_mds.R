#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# Script: 03_calculate_mds.R
# Purpose: Calculate Motif Diversity Score (MDS) for each sample
#          based on motif count files (4-mer frequencies)
#
# Output: A TSV file with two columns: SampleID and MDS
# Author: [Your Name], [Date]
# -------------------------------------------------------------------

# Load necessary libraries (quietly to avoid clutter)
suppressPackageStartupMessages({
  library(entropy)   # For Shannon entropy calculation
  library(dplyr)     # For data manipulation
  library(stringr)   # For string processing
})

# ------------------------------
# Configuration
# ------------------------------
motif_dir <- "fa2"                                # Directory with *_motif_count.txt files
output_file <- "fem_output/mds_summary.tsv"       # Output summary file

# Create output folder if it doesn't exist
dir.create(dirname(output_file), showWarnings = FALSE)

# List all motif count files in the directory
motif_files <- list.files(
  motif_dir,
  pattern = "_motif_count.txt$",
  full.names = TRUE
)

# Generate all valid 4-mer motifs (256 total)
valid_4mers <- apply(expand.grid(rep(list(c("A", "C", "G", "T")), 4)), 1, paste0, collapse = "")

# ------------------------------
# Function: calculate_mds
# ------------------------------
calculate_mds <- function(file) {
  # Read motif counts
  data <- read.table(file, stringsAsFactors = FALSE)
  colnames(data) <- c("count", "motif")

  # Filter for valid ACGT-only 4-mers
  data <- data[grepl("^[ACGT]{4}$", toupper(data$motif)), ]
  data$motif <- toupper(data$motif)

  # Pad missing motifs with count = 0
  full_counts <- setNames(rep(0, length(valid_4mers)), valid_4mers)
  full_counts[data$motif] <- data$count

  # Convert counts to probabilities
  probs <- full_counts / sum(full_counts)

  # Compute normalized Shannon entropy (MDS)
  mds <- entropy(probs, unit = "log2") / log2(256)

  # Extract only the first 4 digits as sample ID
  sample_id <- str_extract(basename(file), "^\\d{4}")

  return(data.frame(SampleID = sample_id, MDS = round(mds, 5)))
}

# ------------------------------
# Main Execution
# ------------------------------
cat("🔍 Calculating MDS for all samples...\n")

# Apply the function to all motif count files
mds_results <- bind_rows(lapply(motif_files, calculate_mds))

# Remove any duplicated SampleIDs (just in case)
mds_results <- mds_results[!duplicated(mds_results$SampleID), ]

# Save result
write.table(
  mds_results,
  output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("✅ MDS summary written to:", output_file, "\n")
