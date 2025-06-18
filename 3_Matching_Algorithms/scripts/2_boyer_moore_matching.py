#!/usr/bin/env python3

import argparse
import csv
import time
import os
import pandas as pd


def parse_arguments():
    """Parse command-line arguments for Boyer-moore matching."""
    parser = argparse.ArgumentParser(
        description='Boyer-moore matching for DNA sequences'
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
        default=("data/2_matching_results/"
                 "boyer_moore_search_results.csv"),
        help='Output CSV file'
    )
    return parser.parse_args()


def read_fasta(file_path):
    """
    Read a FASTA file and return sequence.

    Args:
        file_path (str): Path to the FASTA file

    Returns:
        str: DNA sequence from the FASTA file
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
        file_path (str): Path to the CSV file containing motif data

    Returns:
        list: List of dictionaries containing ham0 motif data
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


def preprocess_bad_char(motif):
    """
    Preprocess motif for the bad character rule.

    Args:
        motif (str): DNA motif sequence

    Returns:
        dict: Dictionary mapping characters to their rightmost positions
    """
    bad_char = {}
    for i in range(len(motif)):
        bad_char[motif[i]] = i
    return bad_char


def preprocess_good_suffix(motif):
    """
    Preprocess motif for the good suffix rule.

    Args:
        motif (str): DNA motif sequence

    Returns:
        list: Array of shift distances for good suffix rule
    """
    m = len(motif)

    # Initialize shift array with maximum possible shift (pattern length)
    # Each position starts with worst-case shift distance
    shift = [m] * m

    # Initialize suffix array to store border lengths
    suffix = [0] * m
    suffix[m - 1] = m  # Last position matches entire remaining pattern

    # Compute suffix match lengths for each position
    # Working backwards through the pattern to find matching suffixes
    for i in range(m - 2, -1, -1):
        j = i
        # Check how far this position matches with the pattern's suffix
        while j >= 0 and motif[j] == motif[m - 1 - i + j]:
            j -= 1
        # Store the length of the matching suffix at position i
        suffix[i] = i - j

    # Compute shift distances using suffix information
    # For each position, calculate minimum safe shift distance
    for i in range(m):
        # Find position where this suffix could align next
        j = m - 1 - suffix[i]
        if j < m:
            # Use minimum of current shift or new calculated shift
            shift[j] = min(shift[j], m - 1 - i)

    return shift


def boyer_moore_search(chromosome, motif):
    """
    Boyer-Moore string matching with bad character and good suffix rules.

    Args:
        chromosome (str): DNA sequence to search in
        motif (str): DNA motif to search for

    Returns:
        list: List of positions where motif matches in chromosome
    """
    matches = []
    n, m = len(chromosome), len(motif)

    # Preprocess motif to build lookup tables for fast shifting
    # bad_char: rightmost position of each character in motif
    # good_suffix: safe shift distances when suffix matches but prefix doesn't
    bad_char = preprocess_bad_char(motif)
    good_suffix = preprocess_good_suffix(motif)

    i = 0

    # Try alignments until motif can't fit in remaining chromosome
    while i <= n - m:
        # Start comparing from end of motif (right-to-left comparison)
        j = m - 1

        # Compare motif with chromosome from right to left
        while j >= 0:
            # If characters don't match, stop comparing at this alignment
            if motif[j] != chromosome[i + j]:
                break
            # Move to next character to the left
            j -= 1

        # Check if we found a complete match
        if j < 0:
            # j went below 0, meaning all characters matched
            matches.append(i)
            i += 1
        else:
            # We hit a mismatch at position j
            # Calculate how far we can safely shift using both heuristics

            # Bad character rule: shift based on mismatched character
            # If mismatched char appears earlier in motif, align with that
            # If it doesn't appear, shift past it entirely
            bc_shift = max(1, j - bad_char.get(chromosome[i + j], -1))

            # Good suffix rule: shift based on matched suffix (if any)
            # Uses precomputed safe distances from good_suffix table
            gs_shift = good_suffix[j]

            # Take the larger shift to ensure we don't miss any matches
            i += max(bc_shift, gs_shift)

    return matches


def run_search(sequence, motifs):
    """
    Run Boyer-Moore search for each motif.

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
        matches = boyer_moore_search(sequence, motif_seq)
        execution_time = time.time() - start_time

        results[name] = {
            'sequence': motif_seq,
            'matches': matches,
            'execution_time': execution_time
        }

        print(f"Motif '{name}': Found {len(matches)} matches "
              f"in {execution_time:.6f}s")

    return results


def verify_results(results, expected_motifs):
    """
    Verify if results match expected positions.

    Args:
        results (dict): Search results from run_search function
        expected_motifs (list): List of expected motif data

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
            false_positives = [pos for pos in found_pos
                               if pos != expected_pos]

            verification[name] = {
                'expected_position': expected_pos,
                'found_positions': found_pos,
                'correctly_found': found,
                'false_positives_count': len(false_positives)
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
    """Main function to execute Boyer-Moore DNA sequence matching."""
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

    print("Boyer-Moore search completed!")


if __name__ == "__main__":
    main()