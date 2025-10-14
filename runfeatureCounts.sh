#!/bin/bash
# =====================================================
# Script: run_featureCounts.sh
# Purpose: Run featureCounts on all alignment files (BAM/SAM) in a 
directory
# Usage: bash run_featureCounts.sh <input_dir> <annotation_file> 
<output_dir>
# =====================================================

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
# Run featureCounts
# -------------------------
featureCounts \
    -T 8 \                                      # Use 8 threads (adjust as 
needed)
    -a "$ANNOTATION" \                          # Annotation file 
(GTF/GFF)
    -o "$OUTPUT_DIR/featureCounts_results.txt" \ # Output file
    "${ALIGN_FILES[@]}"                         # All alignment files


