#!/bin/bash
# =====================================================
# Script: run_featureCounts.sh
# Purpose: Run featureCounts on all alignment files (BAM/SAM) in a directory and save one output file per alignment file
# Usage: bash run_featureCounts.sh <input_dir> <annotation_file> <output_dir>
# =====================================================

set -e  # Stop if any command fails

# -------------------------
# Parse input arguments
# -------------------------
INPUT_DIR=$1
ANNOTATION=$2
OUTPUT_DIR=$3

# -------------------------
# Check arguments
# -------------------------
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_dir> <annotation_file> <output_dir>"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' not found!"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo "Error: Annotation file '$ANNOTATION' not found!"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -------------------------
# Find alignment files (BAM or SAM)
# -------------------------
mapfile -t ALIGN_FILES < <(find "$INPUT_DIR" -type f \( -name "*.bam" -o 
-name "*.sam" \) | sort)

if [ ${#ALIGN_FILES[@]} -eq 0 ]; then
    echo "Error: No BAM or SAM files found in '$INPUT_DIR'."
    exit 1
fi

# -------------------------
# Run featureCounts on each file individually
# -------------------------
for FILE in "${ALIGN_FILES[@]}"; do
    BASENAME=$(basename "$FILE")
    # Remove extension for output file
    SAMPLE_NAME="${BASENAME%%.*}"
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_featureCounts.txt"

    featureCounts \
        -T 8 \
        -a "$ANNOTATION" \
        -o "$OUTPUT_FILE" \
        "$FILE"
done

echo "featureCounts completed for all files!"

