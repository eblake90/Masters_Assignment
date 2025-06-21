#!/bin/bash

# Default paths
DEFAULT_TRIMMED_READS="data/2_trimmed/ENCFF493KQW_trimmed.fastq"
DEFAULT_REFERENCE_DIR="data/3-1_reference"
DEFAULT_GENOME_DIR="data/3-2_star_index"
DEFAULT_OUTPUT_DIR="data/3-3_aligned"
DEFAULT_THREADS=10

# Parse command line arguments
TRIMMED_READS="${1:-$DEFAULT_TRIMMED_READS}"
REFERENCE_DIR="${2:-$DEFAULT_REFERENCE_DIR}"
GENOME_DIR="${3:-$DEFAULT_GENOME_DIR}"
OUTPUT_DIR="${4:-$DEFAULT_OUTPUT_DIR}"
THREADS="${5:-$DEFAULT_THREADS}"

# Create directories
mkdir -p "$GENOME_DIR" "$OUTPUT_DIR" "$REFERENCE_DIR"

echo "Starting STAR alignment for Chr21 and Chr22..."

# Step 1: Download reference genome (Chr21 and Chr22 only)
# Download Chr21 and Chr22 FASTA files using Ensembl URLs
if [ ! -f "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.21.fa" ]; then
    wget -q https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz -O "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
    gunzip "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
fi

if [ ! -f "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.22.fa" ]; then
    wget -q https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz -O "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
    gunzip "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
fi

# Combine Chr21 and Chr22 into single reference
cat "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.21.fa" "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.chromosome.22.fa" > "$REFERENCE_DIR/chr21_22_combined.fa"

echo "Starting STAR genome indexing"
# Step 2: Generate STAR genome index
STAR --runMode genomeGenerate \
    --runThreadN $DEFAULT_THREADS \
    --genomeDir "$GENOME_DIR" \
    --genomeSAindexNbases 12 \
    --genomeFastaFiles "$REFERENCE_DIR/chr21_22_combined.fa"

# Step 3: Align reads using STAR
echo "Starting STAR alignment"
STAR --runMode alignReads \
    --runThreadN $DEFAULT_THREADS \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$TRIMMED_READS" \
    --outFileNamePrefix "$OUTPUT_DIR/chr21_22_"