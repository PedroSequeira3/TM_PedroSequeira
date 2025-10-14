#!/bin/bash

# ============================================
# Script to build HISAT2 genome index
# Usage: ./setupHISAT2.sh <genome_fasta>
# ============================================

GENOME_FASTA=$1
OUTPUT_DIR="$HOME/HISAT2_INDEX"
BASENAME="HISAT2_INDEX"

# --- Check arguments ---
if [ -z "$GENOME_FASTA" ]; then
    echo "Usage: $0 <genome_fasta>"
    exit 1
fi

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: genome FASTA file not found at $GENOME_FASTA"
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Build the HISAT2 index ---
hisat2-build -p 15 "$GENOME_FASTA" "$OUTPUT_DIR/$BASENAME"

echo "HISAT2 index successfully generated in: $OUTPUT_DIR"


