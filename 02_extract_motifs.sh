#!/bin/bash

# Usage: ./02_extract_motifs.sh <REF_FASTA>
# Example: ./02_extract_motifs.sh /path/to/HumControls.fasta

set -euo pipefail

REF=$1
NMER=4
THREADS=6  # Number of parallel jobs, tune this based on your CPU

mkdir -p fa fa2 logs

# Function to process one sample
process_sample() {
    PREFIX="$1"
    R1_BED="five_prime_ends/${PREFIX}_${NMER}bp_r1.bed"
    R2_BED="five_prime_ends/${PREFIX}_${NMER}bp_r2.bed"
    OUT_MERGED="fa2/${PREFIX}_merged_${NMER}bp_motif.bed"
    OUT_COUNTS="fa2/${PREFIX}_motif_count.txt"

    # Skip if output already exists
    if [[ -f "$OUT_COUNTS" ]]; then
        echo "✅ Skipping $PREFIX (already processed)"
        return
    fi

    echo "📖 Extracting sequences from $PREFIX..."

    for BED_FILE in "$R1_BED" "$R2_BED"; do
        if [[ ! -f "$BED_FILE" ]]; then
            echo "❌ Missing BED file: $BED_FILE" >&2
            return
        fi

        ID=$(basename -s .bed "$BED_FILE")
        echo "📘 Processing: $ID"

        bedtools getfasta -fi "$REF" -bed "$BED_FILE" -s -bedOut -fo | \
            awk 'OFS="\t" {print $1, $2, $3, $6, $7, $8, toupper($9)}' > "fa/${ID}_fa.bed"
    done

    echo "📊 Merging and counting motifs for $PREFIX..."

    cat "fa/${PREFIX}_${NMER}bp_r1_fa.bed" "fa/${PREFIX}_${NMER}bp_r2_fa.bed" > "$OUT_MERGED"

    # Assuming motif is in column 6
    awk '{print toupper($6)}' "$OUT_MERGED" | sort | uniq -c | sort -k1,1n > "$OUT_COUNTS"

    echo "✅ Done with $PREFIX"
}

export -f process_sample
export REF NMER

# Extract all sample prefixes based on R1 BED files
PREFIXES=$(find five_prime_ends -name "*_${NMER}bp_r1.bed" | sed -E "s|five_prime_ends/||" | sed -E "s/_${NMER}bp_r1.bed//")

# Run in parallel
echo "$PREFIXES" | parallel -j "$THREADS" process_sample {} 2>&1 | tee logs/extract_motifs.log

echo "🎉 All samples processed."
