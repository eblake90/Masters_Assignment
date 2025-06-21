#!/bin/bash

# Default paths
DEFAULT_SAM_FILE="data/3-3_aligned/chr21_22_Aligned.out.sam"
DEFAULT_GTF_FILE="data/4_annotations/chr21_22_annotation.gtf"
DEFAULT_OUTPUT_DIR="data/5_coverage"

# Parse command line arguments
SAM_FILE="${1:-$DEFAULT_SAM_FILE}"
GTF_FILE="${2:-$DEFAULT_GTF_FILE}"
OUTPUT_DIR="${3:-$DEFAULT_OUTPUT_DIR}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Counting reads per gene with featureCounts..."
featureCounts -a "$GTF_FILE" -o "$OUTPUT_DIR/gene_coverage.txt" "$SAM_FILE"