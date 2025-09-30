#!/bin/bash

# Usage: ./run_multiqc.sh <input_dir> <output_dir>

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run MultiQC
multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

