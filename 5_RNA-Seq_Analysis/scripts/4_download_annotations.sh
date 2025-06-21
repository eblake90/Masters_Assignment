#!/bin/bash

# Default paths
DEFAULT_ANNOTATIONS_DIR="data/4_annotations"

# Parse command line arguments
ANNOTATIONS_DIR="${1:-$DEFAULT_ANNOTATIONS_DIR}"

# Create directory
mkdir -p "$ANNOTATIONS_DIR"

echo "Downloading human gene annotations..."

# Download GTF annotation file from Ensembl (Release 114, GRCh38)
if [ ! -f "$ANNOTATIONS_DIR/Homo_sapiens.GRCh38.114.gtf.gz" ]; then
    wget -q https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz -O "$ANNOTATIONS_DIR/Homo_sapiens.GRCh38.114.gtf.gz"
fi

# Extract Chr21 and Chr22 specific annotations
gunzip -c "$ANNOTATIONS_DIR/Homo_sapiens.GRCh38.114.gtf.gz" | grep -E "^(21|22)" > "$ANNOTATIONS_DIR/chr21_22_annotation.gtf"