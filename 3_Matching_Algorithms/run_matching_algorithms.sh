#!/bin/bash

# Master script to run the entire DNA sequence matching algorithms pipeline

# Parse command line arguments with defaults
DEFAULT_DATA_DIR="data"

# Set parameters with defaults (can be overridden by command line arguments)
DATA_DIR="${1:-${DEFAULT_DATA_DIR}}"
SEED="${2:-42}"  # Random seed for reproducibility
MOTIF_LENGTHS="${3:-15,50,120}"  # Comma-separated motif lengths to generate
FAMILIES_PER_LENGTH="${4:-20}"  # How many families per length
MAX_K="${5:-5}"  # Maximum Hamming-distance mutants to generate

# Create derived paths
INPUT_DATA_DIR="${DATA_DIR}/1_dna_sequence_data"
OUTPUT_RESULTS_DIR="${DATA_DIR}/2_matching_results"
ANALYSIS_DIR="${DATA_DIR}/3_result_analysis"

# Directory structure
mkdir -p "${INPUT_DATA_DIR}"
mkdir -p "${OUTPUT_RESULTS_DIR}"
mkdir -p "${ANALYSIS_DIR}"

# Parameters for approximate matching
MAX_MISMATCHES_RANGE=(0 1 2 3 4 5)  # Runs pigeon-hole with these mismatch values

echo "Starting DNA sequence matching algorithms pipeline..."

# Step 1: Generate DNA sequence with motifs
echo "Running Step 1: Generating DNA sequence with motifs..."
python3 scripts/0_dna_sequence_generator.py \
  --output_dir "${DATA_DIR}/1_dna_sequence_data" \
  --seed ${SEED} \
  --motif_lengths ${MOTIF_LENGTHS} \
  --families_per_length ${FAMILIES_PER_LENGTH} \
  --max_k ${MAX_K}

echo "======================================"
echo " "
echo " "

# Step 2: Run naive string matching
echo "Running Step 2: Running naive string matching..."
python3 scripts/1_naive_matching.py \
  --input_fasta "${DATA_DIR}/1_dna_sequence_data/GRCh38_chr22_with_motifs.fa" \
  --motif_positions "${DATA_DIR}/1_dna_sequence_data/chr22_motif_positions.csv" \
  --output_csv "${DATA_DIR}/2_matching_results/naive_search_results.csv"

echo "======================================"
echo " "
echo " "

# Step 3: Run Boyer-Moore string matching
echo "Running Step 3: Running Boyer-Moore string matching..."
python3 scripts/2_boyer_moore_matching.py \
  --input_fasta "${DATA_DIR}/1_dna_sequence_data/GRCh38_chr22_with_motifs.fa" \
  --motif_positions "${DATA_DIR}/1_dna_sequence_data/chr22_motif_positions.csv" \
  --output_csv "${DATA_DIR}/2_matching_results/boyer_moore_search_results.csv"

echo "======================================"
echo " "
echo " "

# Step 4: Run pigeon-hole approximate matching with different mismatch values
echo "Running Step 4: Running pigeon-hole approximate matching with multiple mismatch values..."

for MISMATCH in "${MAX_MISMATCHES_RANGE[@]}"; do
  echo "  Running with max_mismatches=${MISMATCH}..."
  python3 scripts/3_pigeon_hole_matching.py \
    --input_fasta "${DATA_DIR}/1_dna_sequence_data/GRCh38_chr22_with_motifs.fa" \
    --motif_positions "${DATA_DIR}/1_dna_sequence_data/chr22_motif_positions.csv" \
    --output_dir "${DATA_DIR}/2_matching_results" \
    --max_mismatches ${MISMATCH}

  echo "  Completed run with max_mismatches=${MISMATCH}"
  echo "  --------------------------------------"
done

# Step 5: Compile and analyze results
echo "Running Step 5: Compiling results..."
python3 scripts/4_compile_results.py \
  --naive_results "${DATA_DIR}/2_matching_results/naive_search_results.csv" \
  --boyer_moore_results "${DATA_DIR}/2_matching_results/boyer_moore_search_results.csv" \
  --pigeon_hole_dir "${DATA_DIR}/2_matching_results" \
  --output_file "${DATA_DIR}/3_result_analysis/algorithm_comparison.csv" \
  --max_mismatch_level ${MAX_K}

echo "Running Step 6: Generating summary statistics..."
python3 scripts/5_summary_stats.py \
  --input_file "${DATA_DIR}/3_result_analysis/algorithm_comparison.csv" \
  --output_dir "${DATA_DIR}/3_result_analysis" \
  --max_mismatch_level ${MAX_K}

echo "Pipeline completed successfully!"