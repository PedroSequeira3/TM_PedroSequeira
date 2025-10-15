#!/bin/bash
# ============================================
# Script to run Salmon quantification in paired-end mode
# Usage: ./runSalmon.sh <input_dir> <output_dir>
# Salmon index must exist at ~/SALMON_INDEX or will be built automatically
# ============================================

INPUT_DIR=$1
OUTPUT_DIR=$2
SALMON_INDEX="$HOME/SALMON_INDEX"

# --- Check usage ---
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory '$INPUT_DIR' not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Check if Salmon index exists ---
if [ ! -d "$SALMON_INDEX" ] || [ -z "$(ls -A "$SALMON_INDEX" 2>/dev/null)" 
]; then
    echo "No Salmon index found in '$SALMON_INDEX', creating one."

    SETUP_SCRIPT="$(dirname "$0")/setupSalmon.sh"
    if [ ! -x "$SETUP_SCRIPT" ]; then
        echo "Error: setupSalmon.sh not found or not executable in $(dirname "$0")"
        echo "Please make sure it's in the same directory and has execute permission."
        exit 1
    fi

    echo "Provide the path to the transcriptome FASTA file:"
    read TRANSCRIPTOME_FASTA

    "$SETUP_SCRIPT" "$TRANSCRIPTOME_FASTA"
else
    echo "Using existing Salmon index at '$SALMON_INDEX'"
fi

# --- Run Salmon on paired-end reads ---
for R1 in "$INPUT_DIR"/*_R1*.fastq* "$INPUT_DIR"/*_R1*.fq* 
"$INPUT_DIR"/*_1*.fastq* "$INPUT_DIR"/*_1*.fq*; do
    [ -e "$R1" ] || continue  # skip if none match
    SAMPLE=$(basename "$R1" | sed -E 's/_R1.*|_1.*//')
    R2=$(find "$INPUT_DIR" -type f \( -name "${SAMPLE}_R2*.fastq*" -o 
-name "${SAMPLE}_2*.fastq*" \) | head -n 1)

    if [ -z "$R2" ]; then
        continue
    fi

    OUT_SAMPLE_DIR="$OUTPUT_DIR/${SAMPLE}_salmon"
    mkdir -p "$OUT_SAMPLE_DIR"

    salmon quant -i "$SALMON_INDEX" \
                 -l IU \
                 --threads 15 \
                 --validateMappings \
                 -1 "$R1" -2 "$R2" \
                 -o "$OUT_SAMPLE_DIR" \
                 2> "$OUT_SAMPLE_DIR/${SAMPLE}_log.txt"

done

echo "All samples processed successfully!"

