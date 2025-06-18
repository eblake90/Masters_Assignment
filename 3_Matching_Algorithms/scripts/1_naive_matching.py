#!/usr/bin/env python3

import argparse
import csv
import time
import os
import pandas as pd


def parse_arguments():
    """Parse command-line arguments for naive string matching."""
    parser = argparse.ArgumentParser(
        description='Naive string matching for DNA sequences'
    )
    parser.add_argument(
        '--input_fasta',
        default=("data/1_dna_sequence_data/GRCh38_chr22_with_motifs.fa"),
        help='FASTA file path'
    )
    parser.add_argument(
        '--motif_positions',
        default=("data/1_dna_sequence_data/chr22_motif_positions.csv"),
        help='Motif positions CSV'
    )
    parser.add_argument(
        '--output_csv',
        default=("data/2_matching_results/naive_search_results.csv"),
        help='Output CSV file'
    )
    return parser.parse_args()


def read_fasta(file_path):
    """
    Read a FASTA file and return sequence.

    Args:
        file_path (str): Path to FASTA file

    Returns:
        str: DNA sequence from FASTA file
    """
    print(f"Reading FASTA: {file_path}")
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines[1:])
    print(f"Loaded {len(sequence)} nucleotides")
    return sequence


def read_motifs(file_path):
    """
    Read motifs from CSV file, filtering for _ham0 motifs only.

    Args:
        file_path (str): Path to CSV file containing motif data

    Returns:
        list: List of dictionaries containing ham0 motif information
    """
    print(f"Reading motifs: {file_path}")
    motifs = []
    ham0_motifs = []

    with open(file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            motifs.append({
                'name': row['motif_name'],
                'start': int(row['start']),
                'sequence': row['sequence']
            })

            # Filter to only include motifs with "_ham0" in their name
            if "_ham0" in row['motif_name']:
                ham0_motifs.append({
                    'name': row['motif_name'],
                    'start': int(row['start']),
                    'sequence': row['sequence']
                })

    print(f"    Loaded {len(motifs)} total motifs")
    print(f"    Filtered to {len(ham0_motifs)} ham0 motifs")
    return ham0_motifs


def naive_search(chromosome, motif):
    """
    Naive string matching.

    Args:
        chromosome (str): DNA sequence to search in
        motif (str): DNA motif to search for

    Returns:
        list: List of starting positions where motif matches
    """
    matches = []
    n, m = len(chromosome), len(motif)

    # Try every possible starting position in chromosome
    for i in range(n - m + 1):
        # Assume we have a match until proven otherwise
        match = True
        # Check each character of motif against chromosome at current position
        for j in range(m):
            # If any character doesn't match, this position is not a match
            if chromosome[i + j] != motif[j]:
                match = False
                break

        # If we checked all characters without finding mismatch,
        # record this position
        if match:
            matches.append(i)

    return matches


def run_search(sequence, motifs):
    """
    Run search for each motif.

    Args:
        sequence (str): DNA sequence to search in
        motifs (list): List of motif dictionaries to search for

    Returns:
        dict: Dictionary containing search results for each motif
    """
    results = {}
    for motif in motifs:
        name = motif['name']
        motif_seq = motif['sequence']

        start_time = time.time()
        matches = naive_search(sequence, motif_seq)
        execution_time = time.time() - start_time

        results[name] = {
            'sequence': motif_seq,
            'matches': matches,
            'execution_time': execution_time
        }

        print(f"Motif '{name}': Found {len(matches)} matches in "
              f"{execution_time:.6f}s")

    return results


def verify_results(results, expected_motifs):
    """
    Verify if results match expected positions.

    Args:
        results (dict): Search results from run_search function
        expected_motifs (list): List of expected motif dictionaries

    Returns:
        dict: Dictionary containing verification results for each motif
    """
    verification = {}
    for motif in expected_motifs:
        name = motif['name']
        expected_pos = motif['start']

        if name in results:
            found_pos = results[name]['matches']
            found = expected_pos in found_pos

            verification[name] = {
                'expected_position': expected_pos,
                'found_positions': found_pos,
                'correctly_found': found
            }

    return verification


def save_results(results, verification, output_file):
    """Save results to CSV file using pandas."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Prepare data for DataFrame
    output_data = []
    for name, search_data in results.items():
        if name in verification:
            verify_data = verification[name]

            output_data.append({
                'motif_name': name,
                'motif_sequence': search_data['sequence'],
                'expected_position': verify_data['expected_position'],
                'found_positions': ','.join(map(str, search_data['matches'])),
                'correctly_found': verify_data['correctly_found'],
                'execution_time': f"{search_data['execution_time']:.6f}"
            })

    # Create DataFrame and save
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")


def main():
    args = parse_arguments()

    # Read sequence and motifs
    sequence = read_fasta(args.input_fasta)
    motifs = read_motifs(args.motif_positions)

    # Run search
    search_results = run_search(sequence, motifs)

    # Verify results
    verification = verify_results(search_results, motifs)

    # Save results
    save_results(search_results, verification, args.output_csv)

    print("Naive search completed!")


if __name__ == "__main__":
    main()