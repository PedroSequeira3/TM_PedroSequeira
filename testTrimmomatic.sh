#!/bin/bash

# ==========================================================
# Trimmomatic batch script (paired-only)
# Usage:
#   ./testTrimmomatic.sh INPUT_DIR OUTPUT_DIR INTERVAL QUALITY LENGTH
#
# Keeps only paired reads (_1P and _2P).
# ==========================================================

# Check for correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR INTERVAL QUALITY LENGTH"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
INTERVAL="$3"
QUALITY="$4"
LENGTH="$5"

# Create output dir if needed
mkdir -p "$OUTPUT_DIR"

# Build a sorted array of files
FILES=($(ls "$INPUT_DIR"/*.fq.gz | sort))

# Check if number of files is even
if (( ${#FILES[@]} % 2 != 0 )); then
    echo "Error: Odd number of FASTQ files found. Please check pairing."
    exit 1
fi

ADAPTERS="/opt/tools/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

for ((i=0; i<${#FILES[@]}; i+=2)); do
    R1="${FILES[$i]}"
    R2="${FILES[$i+1]}"
    BASE=$(basename "$R1" .fq.gz)

    echo "Processing: $R1 and $R2 -> $BASE"
    
OUT1="$OUTPUT_DIR/${BASE}_Q${QUALITY}_L${LENGTH}_S${INTERVAL}_R1P.fq.gz"
    
OUT2="$OUTPUT_DIR/${BASE}_Q${QUALITY}_L${LENGTH}_S${INTERVAL}_R2P.fq.gz"

    trimmomatic PE -threads 15 -phred33 \
        "$R1" "$R2" \
        "$OUT1" /dev/null \
        "$OUT2" /dev/null \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 \
        SLIDINGWINDOW:"$INTERVAL":"$QUALITY" MINLEN:"$LENGTH" \
        > "$OUTPUT_DIR/${BASE}_trimmomatic.log" 2>&1
done

