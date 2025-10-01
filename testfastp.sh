#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR INTERVAL QUALITY LENGTH"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
INTERVAL="$3"
QUALITY="$4"
LENGTH="$5"

# Create output dir if needed
mkdir -p "$OUTPUT_DIR"

# Build a sorted array of files
FILES=($(ls "$INPUT_DIR"/*.fq.gz | sort))

# Check if number of files is even
if (( ${#FILES[@]} % 2 != 0 )); then
    echo "Error: Odd number of FASTQ files found. Make sure files are 
paired."
    exit 1
fi

# Loop two by two
for ((i=0; i<${#FILES[@]}; i+=2)); do
    R1="${FILES[$i]}"
    R2="${FILES[$i+1]}"

    # Create a sample base name (take common prefix of the two files)
    BASE=$(basename "$R1" .fq.gz)

    echo "Processing: $R1 and $R2 -> $BASE"

    fastp -i "$R1" -I "$R2" \
    -o 
"$OUTPUT_DIR"/"${BASE}"_Q"$QUALITY"_L"$LENGTH"_S"$INTERVAL"_R1.fq.gz \
    -O 
"$OUTPUT_DIR"/"${BASE}"_Q"$QUALITY"_L"$LENGTH"_S"$INTERVAL"_R2.fq.gz \
    --detect_adapter_for_pe \
    --length_required "$LENGTH" --low_complexity_filter \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
    --cut_right --cut_right_window_size "$INTERVAL" 
--cut_right_mean_quality "$QUALITY" \
    --html 
"$OUTPUT_DIR"/"${BASE}"_Q"$QUALITY"_L"$LENGTH"_S"$INTERVAL"_report.html \
    --thread 15 > "$OUTPUT_DIR"/"${BASE}"_fastp.log 2>&1
done

