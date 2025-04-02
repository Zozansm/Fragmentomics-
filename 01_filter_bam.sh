#!/bin/bash

# ---------------------------------------------
# Script: 01_filter_bam.sh
# Description: Filters BAM files for cfDNA analysis.
#              Extracts properly paired, autosomal reads with MAPQ ≥ 30,
#              generates BEDPE with fragment size, and extracts 5′ end 4-mers.
#
# Usage: ./01_filter_bam.sh /path/to/bam_directory
# ---------------------------------------------

# Check if input directory is provided
if [ -z "$1" ]; then
  echo "❌ Error: Please provide the path to the BAM directory."
  echo "Usage: $0 /path/to/bam_directory"
  exit 1
fi

input_dir=$1
nmer=4  # Length of 5′ end to extract

# Create output directories
mkdir -p filtered_bam        # For filtered BAM files
mkdir -p bedpe               # For BEDPE with fragment sizes
mkdir -p five_prime_ends     # For R1 and R2 4bp BED files

# Loop through only sorted BAM files (these contain full data)
for BAM in "$input_dir"/*.sorted.bam; do
  [ -e "$BAM" ] || continue  # Skip if no .bam files found

  NAME=$(basename "$BAM")           # Extract filename (e.g., 7002.sorted.bam)
  BASE=${NAME//.bam/}               # Remove .bam extension

  echo "🔍 Processing: $BASE"

  # Set output file paths
  FILTERED_BAM="filtered_bam/${BASE}.filtered.bam"
  BEDPE_OUT="bedpe/${BASE}.filtered.bedpe"
  R1_OUT="five_prime_ends/${BASE}_4bp_r1.bed"
  R2_OUT="five_prime_ends/${BASE}_4bp_r2.bed"

  # STEP 1: Filter properly paired, autosomal reads with MAPQ ≥ 30
  samtools view -bh -f 2 -F 3844 -q 30 "$BAM" $(seq 1 22 | sed 's/^/chr/') > "$FILTERED_BAM"

  # STEP 2: Convert to BEDPE and calculate fragment size
  samtools sort -n "$FILTERED_BAM" | \
    bedtools bamtobed -i - -bedpe | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $6-$2}' > "$BEDPE_OUT"

  # STEP 3: Extract 4bp 5′ end sequences for R1 and R2
  awk -v nmer="$nmer" 'OFS="\t" {print $1, $2, $2+nmer, $7, $8, $9, $11, $12}' "$BEDPE_OUT" > "$R1_OUT"
  awk -v nmer="$nmer" 'OFS="\t" {print $4, $6-nmer, $6, $7, $8, $10, $11, $12}' "$BEDPE_OUT" > "$R2_OUT"

  # Log output
  echo "✅ Completed: $BASE"
  echo "   ➤ Filtered BAM:       $FILTERED_BAM"
  echo "   ➤ BEDPE + sizes:      $BEDPE_OUT"
  echo "   ➤ 5′ ends (R1/R2):     $R1_OUT, $R2_OUT"
done
