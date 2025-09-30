#!/bin/bash

# Usage: ./run_fastqc.sh <input_dir> <output_dir>

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC for FASTQ files
for file in "$INPUT_DIR"/*.{fq,fq.gz,fastq,fastq.gz}; do
    [ -e "$file" ] || continue 
    echo "Processing: $file"
    fastqc "$file" -o "$OUTPUT_DIR"
done


