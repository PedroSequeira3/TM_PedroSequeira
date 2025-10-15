#!/bin/bash

# ============================================
# Script to build Salmon transcriptome index
# Usage: ./setupSALMON.sh <transcriptome_fasta>
# ============================================

TRANSCRIPTOME_FASTA=$1
OUTPUT_DIR="$HOME/SALMON_INDEX"

# --- Check arguments ---
if [ -z "$TRANSCRIPTOME_FASTA" ]; then
    echo "Usage: $0 <transcriptome_fasta>"
    exit 1
fi

if [ ! -f "$TRANSCRIPTOME_FASTA" ]; then
    echo "Error: transcriptome FASTA file not found at $TRANSCRIPTOME_FASTA"
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Build the Salmon index ---
salmon index -t "$TRANSCRIPTOME_FASTA" -i "$OUTPUT_DIR" -p 15

echo "Salmon index successfully generated in: $OUTPUT_DIR"


