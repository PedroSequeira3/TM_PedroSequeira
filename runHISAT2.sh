#!/bin/bash
# ============================================
# Script to run HISAT2 in paired-end mode
# Usage: ./runHISAT2.sh <input_dir> <output_dir>
# HISAT2 index must exist at ~/HISAT2_INDEX or will be built automatically
# ============================================

INPUT_DIR=$1
OUTPUT_DIR=$2
HISAT2_INDEX="$HOME/HISAT2_INDEX"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory '$INPUT_DIR' not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Check if HISAT2 index exists ---
if [ ! -d "$HISAT2_INDEX" ] || [ -z "$(ls -A "$HISAT2_INDEX" 2>/dev/null)" 
]; then
    echo "No HISAT2 index found in '$HISAT2_INDEX', creating one."

    SETUP_SCRIPT="$(dirname "$0")/setupHISAT2.sh"
    if [ ! -x "$SETUP_SCRIPT" ]; then
        echo "Error: setupHISAT2.sh not found or not executable in $(dirname "$0")"
        echo "Please make sure it's in the same directory and has execute permission."
        exit 1
    fi

    echo "Provide the path to the genome FASTA file:"
    read GENOME_FASTA
    echo "Provide the path to the annotation GTF file (optional, press Enter to skip):"
    read GTF_FILE

    "$SETUP_SCRIPT" "$GENOME_FASTA" "$GTF_FILE"
else
    echo "Using existing HISAT2 index at '$HISAT2_INDEX'"
fi

# --- Run HISAT2 on paired-end reads ---
for R1 in "$INPUT_DIR"/*_R1*.fastq* "$INPUT_DIR"/*_R1*.fq* 
"$INPUT_DIR"/*_1*.fastq* "$INPUT_DIR"/*_1*.fq*; do
    [ -e "$R1" ] || continue  # skip if none match
    SAMPLE=$(basename "$R1" | sed -E 's/_R1.*|_1.*//')
    R2=$(find "$INPUT_DIR" -type f \( -name "${SAMPLE}_R2*.fastq*" -o 
-name "${SAMPLE}_2*.fastq*" \) | head -n 1)

    if [ -z "$R2" ]; then
        continue
    fi

    hisat2 -x "$HISAT2_INDEX/bowtie_index" \
           -1 "$R1" -2 "$R2" \
           -p 15 \
           -S "$OUTPUT_DIR/${SAMPLE}.sam" \
           2> "$OUTPUT_DIR/${SAMPLE}_log.txt"

    # Convert SAM to sorted BAM
    samtools view -@ 15 -bS "$OUTPUT_DIR/${SAMPLE}.sam" | samtools sort -@ 15 -o "$OUTPUT_DIR/${SAMPLE}_sorted.bam"
    samtools index "$OUTPUT_DIR/${SAMPLE}_sorted.bam"
    rm "$OUTPUT_DIR/${SAMPLE}.sam"

    echo "Finished $SAMPLE"
done

echo "All samples processed successfully!"

