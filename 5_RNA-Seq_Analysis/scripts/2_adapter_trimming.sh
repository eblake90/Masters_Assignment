#!/bin/bash

# Default paths
DEFAULT_INPUT_FILE="data/0_fastq_input/ENCFF493KQW.fastq"
DEFAULT_OUTPUT_DIR="data/2_trimmed"
DEFAULT_ADAPTER_FILE="$CONDA_PREFIX/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa"

# Parse command line arguments
INPUT_FILE="${1:-$DEFAULT_INPUT_FILE}"
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT_DIR}"
ADAPTER_FILE="${3:-$DEFAULT_ADAPTER_FILE}"

# Set paths
OUTPUT_FILE="$OUTPUT_DIR/ENCFF493KQW_trimmed.fastq"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting adapter trimming..."

# Run Trimmomatic
# Remove TruSeq3 adapters: 2 mismatches allowed, palindrome disabled (using single-end reads), min score 10
# Remove low-quality bases (Q≤3) from read start
# Remove low-quality bases (Q≤3) from read end
# Cut when 4-base window average quality drops below 15
# Discard reads shorter than 95bp after trimming
trimmomatic SE \
    "$INPUT_FILE" \
    "$OUTPUT_FILE" \
    ILLUMINACLIP:"$ADAPTER_FILE":2:0:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:95

echo ""
echo "Trimming Statistics:"
echo "Original reads:"
wc -l "$INPUT_FILE" | awk '{print $1/4}' | cut -d'.' -f1
echo ""
echo "Trimmed reads:"
wc -l "$OUTPUT_FILE" | awk '{print $1/4}' | cut -d'.' -f1