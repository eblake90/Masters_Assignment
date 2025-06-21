#!/bin/bash

# FastQC Quality Control Script for RNA-seq data
# This script performs quality control and identifies adapter sequences

# Default paths
DEFAULT_INPUT_FILE="data"
DEFAULT_OUTPUT_DIR="data/1_fastqc_results"

# Parse command line arguments
INPUT_FILE="${1:-$DEFAULT_INPUT_FILE}"
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT_DIR}"

# Creating output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Running FastQC on FASTQ file
echo "Processing: $(basename "$INPUT_FILE")"
fastqc "$INPUT_FILE" --outdir "$OUTPUT_DIR"

echo ""
echo "FastQC analysis complete!"
echo "Results saved to: $OUTPUT_DIR"