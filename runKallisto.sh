#!/bin/bash
# ============================================
# Script to run Kallisto quantification in paired-end mode
# Usage: ./runKallisto.sh <input_dir> <output_dir>
# Kallisto index must exist at ~/KALLISTO_INDEX or will be built 
automatically
# ============================================

INPUT_DIR=$1
OUTPUT_DIR=$2
KALLISTO_INDEX="$HOME/KALLISTO_INDEX/kallisto_index"

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

# --- Check if Kallisto index exists ---
if [ ! -f "$KALLISTO_INDEX" ]; then
    echo "No Kallisto index found at '$KALLISTO_INDEX', creating one."

    SETUP_SCRIPT="$(dirname "$0")/setupKallisto.sh"
    if [ ! -x "$SETUP_SCRIPT" ]; then
        echo "Error: setupKallisto.sh not found or not executable in $(dirname "$0")"
        echo "Please make sure it's in the same directory and has execute permission."
        exit 1
    fi

    echo "Provide the path to the transcriptome FASTA file:"
    read TRANSCRIPTOME_FASTA

    "$SETUP_SCRIPT" "$TRANSCRIPTOME_FASTA"
else
    echo "Using existing Kallisto index at '$KALLISTO_INDEX'"
fi

# --- Run Kallisto on paired-end reads ---
for R1 in "$INPUT_DIR"/*_R1*.fastq* "$INPUT_DIR"/*_R1*.fq* 
"$INPUT_DIR"/*_1*.fastq* "$INPUT_DIR"/*_1*.fq*; do
    [ -e "$R1" ] || continue  # skip if none match
    SAMPLE=$(basename "$R1" | sed -E 's/_R1.*|_1.*//')
    R2=$(find "$INPUT_DIR" -type f \( -name "${SAMPLE}_R2*.fastq*" -o 
-name "${SAMPLE}_2*.fastq*" \) | head -n 1)

    if [ -z "$R2" ]; then
        continue
    fi

    OUT_SAMPLE_DIR="$OUTPUT_DIR/${SAMPLE}_kallisto"
    mkdir -p "$OUT_SAMPLE_DIR"

    kallisto quant -i "$KALLISTO_INDEX" \
                   -t 20 \
                   -o "$OUT_SAMPLE_DIR" \
                   "$R1" "$R2" \
                   2> "$OUT_SAMPLE_DIR/${SAMPLE}_log.txt"

done

echo "All samples processed successfully!"

