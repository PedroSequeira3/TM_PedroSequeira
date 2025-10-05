#!/bin/bash
# ==========================================================
# FastQ Screen Auto-Setup and Run Script
# ==========================================================
# Usage: ./runFastQScreen.sh INPUT_DIR OUTPUT_DIR
# ==========================================================

# --- Check arguments ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# --- Configuration ---
THREADS=8
HOME_GENOME_DIR="$HOME/FastQ_Screen_Genomes"
CONFIG_FILE="$HOME_GENOME_DIR/fastq_screen.conf"

# --- Check or generate FastQ_Screen_Genomes directory ---
if [ ! -d "$HOME_GENOME_DIR" ]; then
    echo “No FastQ_Screen_Genomes directory found in \$HOME."
    
    fastq_screen --get_genomes

    # Check success
    if [ ! -d "$HOME_GENOME_DIR" ]; then
        echo "Failed to generate FastQ_Screen_Genomes directory."
        exit 1
    fi
    echo "FastQ_Screen_Genomes directory created at: $HOME_GENOME_DIR"
else
    echo "Found FastQ_Screen_Genomes directory in \$HOME."
fi

# --- Check config file ---
if [ ! -f "$CONFIG_FILE" ]; then
    echo "No fastq_screen.conf found in FastQ_Screen_Genomes.”
exit 1
else
    echo "Using existing config file: $CONFIG_FILE"
fi

# --- Prepare output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Find FASTQ files ---
FILES=($(ls "$INPUT_DIR"/*.fq.gz 2>/dev/null))

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No FASTQ files found in $INPUT_DIR"
    exit 1
fi

# --- Run FastQ Screen on each file ---
for FILE in "${FILES[@]}"; do
    BASENAME=$(basename "$FILE" .fq.gz)

    fastq_screen \
        --conf "$CONFIG_FILE" \
        --threads "$THREADS" \
        --outdir "$OUTPUT_DIR" \
        "$FILE" > "$OUTPUT_DIR/${BASENAME}_fastqscreen.log" 2>&1

done

