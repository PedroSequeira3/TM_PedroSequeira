#!/bin/bash

# ============================================
# Script to build Bowtie2 genome index
# Usage: ./setupBOWTIE2.sh <genome_fasta>
# ============================================

GENOME_FASTA=$1
OUTPUT_DIR="$HOME/BOWTIE2_INDEX"
BASENAME="BOWTIE2_INDEX"

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

# --- Build the Bowtie2 index ---
bowtie2-build --threads 15 "$GENOME_FASTA" "$OUTPUT_DIR/$BASENAME"

echo "Bowtie2 index successfully generated in: $OUTPUT_DIR"


