#!/bin/bash

# Master script to run the entire RNA-seq analysis pipeline

# Default values
DEFAULT_DATA_DIR="data"
DEFAULT_FASTQ_URL="https://www.encodeproject.org/files/ENCFF493KQW/@@download/ENCFF493KQW.fastq.gz"
DEFAULT_THREADS=10

# Set parameters with defaults (can be overridden by command line arguments)
DATA_DIR="${1:-${DEFAULT_DATA_DIR}}"
FASTQ_URL="${2:-${DEFAULT_FASTQ_URL}}"
THREADS="${3:-${DEFAULT_THREADS}}"

# Directory structure
mkdir -p "${DATA_DIR}/0_fastq_input"
mkdir -p "${DATA_DIR}/1_fastqc_results"
mkdir -p "${DATA_DIR}/2_trimmed"
mkdir -p "${DATA_DIR}/3-1_reference"
mkdir -p "${DATA_DIR}/3-2_star_index"
mkdir -p "${DATA_DIR}/3-3_aligned"
mkdir -p "${DATA_DIR}/4_annotations"
mkdir -p "${DATA_DIR}/5_coverage"
mkdir -p "${DATA_DIR}/6_normalized"

# File paths
FASTQ_INPUT_DIR="${DATA_DIR}/0_fastq_input"
INPUT_FASTQ="${FASTQ_INPUT_DIR}/ENCFF493KQW.fastq"
FASTQC_OUTPUT="${DATA_DIR}/1_fastqc_results"
TRIMMED_OUTPUT="${DATA_DIR}/2_trimmed"
TRIMMED_FASTQ="${TRIMMED_OUTPUT}/ENCFF493KQW_trimmed.fastq"
REFERENCE_DIR="${DATA_DIR}/3-1_reference"
GENOME_INDEX_DIR="${DATA_DIR}/3-2_star_index"
ALIGNMENT_OUTPUT="${DATA_DIR}/3-3_aligned"
ANNOTATIONS_DIR="${DATA_DIR}/4_annotations"
COVERAGE_OUTPUT="${DATA_DIR}/5_coverage"
SAM_FILE="${ALIGNMENT_OUTPUT}/chr21_22_Aligned.out.sam"
GTF_FILE="${ANNOTATIONS_DIR}/chr21_22_annotation.gtf"

# Adapter file path
ADAPTER_FILE="$CONDA_PREFIX/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa"

echo "Starting RNA-seq analysis pipeline for Chr21 and Chr22..."

# Step 0: Download FASTQ file
echo "Running Step 0: Downloading FASTQ file from ENCODE..."
python3 scripts/0_fastq_downloader.py \
  --url "$FASTQ_URL" \
  --output_dir "$FASTQ_INPUT_DIR"

echo "======================================"
echo " "
echo " "

# Step 1: Quality Control with FastQC
echo "Running Step 1: Quality control with FastQC..."
bash scripts/1_QC.sh \
  "$INPUT_FASTQ" \
  "$FASTQC_OUTPUT"

echo "======================================"
echo " "
echo " "

# Step 2: Adapter trimming with Trimmomatic
echo "Running Step 2: Adapter trimming with Trimmomatic..."
bash scripts/2_adapter_trimming.sh \
  "$INPUT_FASTQ" \
  "$TRIMMED_OUTPUT" \
  "$ADAPTER_FILE"

echo "======================================"
echo " "
echo " "

# Step 3: STAR alignment
echo "Running Step 3: STAR alignment..."
bash scripts/3_star_alignment.sh \
  "$TRIMMED_FASTQ" \
  "$REFERENCE_DIR" \
  "$GENOME_INDEX_DIR" \
  "$ALIGNMENT_OUTPUT" \
   $THREADS
echo "======================================"
echo " "
echo " "

# Step 4: Download annotations
echo "Running Step 4: Downloading gene annotations..."
bash scripts/4_download_annotations.sh \
  "$ANNOTATIONS_DIR"

echo "======================================"
echo " "
echo " "

# Step 5: Gene coverage counting
echo "Running Step 5: Counting gene coverage with featureCounts..."
bash scripts/5_gene_coverage.sh \
  "$SAM_FILE" \
  "$GTF_FILE" \
  "$COVERAGE_OUTPUT"

echo "======================================"
echo " "
echo " "

# Step 6: Expression normalization
echo "Running Step 6: Normalizing gene expression..."
python3 scripts/6_normalise_expression.py \
  --coverage_file "${COVERAGE_OUTPUT}/gene_coverage.txt" \
  --gtf_file "$GTF_FILE" \
  --output_dir "${DATA_DIR}/6_normalized"

echo "Pipeline completed successfully!"