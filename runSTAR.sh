#!/bin/bash
# ============================================
# Script to run STAR in paired-end mode
# Usage: ./setupSTAR.sh <input_dir> <output_dir>
# STAR index must exist at ~/STAR_INDEX or will be built automatically
# ============================================

INPUT_DIR=$1
OUTPUT_DIR=$2
STAR_INDEX="$HOME/STAR_INDEX"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory '$INPUT_DIR' not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Check if STAR index exists ---
if [ ! -d "$STAR_INDEX" ] || [ -z "$(ls -A "$STAR_INDEX" 2>/dev/null)" ]; 
then
    echo "No STAR index found in '$STAR_INDEX, creating one.”

    SETUP_SCRIPT="$(dirname "$0”)/setupSTAR.sh"
    if [ ! -x "$SETUP_SCRIPT" ]; then
        echo "Error: setupSTAR.sh not found or not executable in $(dirname "$0")"
        echo "Please make sure it's in the same directory and has execute permission."
        exit 1
    fi

    # Ask for the genome and annotation files
    echo "Provide the path to the genome FASTA file:"
    read GENOME_FASTA
    echo "Provide the path to the annotation GTF file:"
    read GTF_FILE

    # Run the setup script to build the index
    "$SETUP_SCRIPT" "$GENOME_FASTA" "$GTF_FILE"
else
    echo "Using existing STAR index at '$STAR_INDEX'"
fi

# --- Run STAR on paired-end reads ---
for R1 in "$INPUT_DIR"/*_R1*.fastq* "$INPUT_DIR"/*_R1*.fq* \
          "$INPUT_DIR"/*_P1*.fastq* "$INPUT_DIR"/*_P1*.fq*; do
    [ -e "$R1" ] || continue  # skip if none match
    SAMPLE=$(basename "$R1" | sed -E 's/_R1.*|_1.*//')
    R2=$(find "$INPUT_DIR" -type f \( -name "${SAMPLE}_R2*.fastq*" -o 
-name "${SAMPLE}_2*.fastq*" \) | head -n 1)

    if [ -z "$R2" ]; then
        continue
    fi

    STAR --genomeDir "$STAR_INDEX" \
         --readFilesCommand zcat \
         --readFilesIn "$R1" "$R2" \
         --runThreadN 15 \
         --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --outStd Log

done

echo "All samples processed successfully!"

