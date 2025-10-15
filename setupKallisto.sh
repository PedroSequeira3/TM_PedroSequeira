#!/bin/bash

# ============================================
# Script to build Kallisto transcriptome index
# Usage: ./setupKALLISTO.sh <transcriptome_fasta>
# ============================================

TRANSCRIPTOME_FASTA=$1
OUTPUT_DIR="$HOME/KALLISTO_INDEX"
INDEX_NAME="kallisto_index.idx"

# --- Check arguments ---
if [ -z "$TRANSCRIPTOME_FASTA" ]; then
    echo "Usage: $0 <transcriptome_fasta>"
    exit 1
fi

if [ ! -f "$TRANSCRIPTOME_FASTA" ]; then
    echo "Error: transcriptome FASTA file not found at 
$TRANSCRIPTOME_FASTA"
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Build the Kallisto index ---
kallisto index -i "$OUTPUT_DIR/$INDEX_NAME" "$TRANSCRIPTOME_FASTA"

echo "Kallisto index successfully generated in: $OUTPUT_DIR"

