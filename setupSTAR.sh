#!/bin/bash

# ============================================
# Script to build STAR genome index
# Usage: ./setupSTAR.sh <genome_fasta> <annotation_gtf>
# ============================================

GENOME_FASTA=$1
ANNOTATION_GTF=$2
OUTPUT_DIR="$HOME/STAR_INDEX"

# --- Check arguments ---
if [ -z "$GENOME_FASTA" ] || [ -z "$ANNOTATION_GTF" ]; then
    echo "Usage: $0 <genome_fasta> <annotation_gtf>"
    exit 1
fi

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: genome FASTA file not found at $GENOME_FASTA"
    exit 1
fi

if [ ! -f "$ANNOTATION_GTF" ]; then
    echo "Error: annotation GTF file not found at $ANNOTATION_GTF"
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Build the STAR index ---
STAR --runThreadN 15 \
     --runMode genomeGenerate \
     --genomeDir "$OUTPUT_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$ANNOTATION_GTF" \
     --sjdbOverhang 100 \
     --limitGenomeGenerateRAM 40000000000

echo "STAR index successfully generated in: $OUTPUT_DIR"

