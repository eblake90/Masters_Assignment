#!/usr/bin/env python3

import argparse
import csv
import time
import pandas as pd
import os
from collections import defaultdict


def parse_arguments():
    """Parse command-line arguments for pigeon-hole string matching."""
    parser = argparse.ArgumentParser(
        description='Pigeon-hole matching for DNA sequences'
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
        '--output_dir',
        default=("data/2_matching_results"),
        help='Output directory for CSV files'
    )
    parser.add_argument(
        '--max_mismatches',
        type=int,
        default=5,
        help='Maximum number of mismatches allowed'
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
    Read motifs from CSV file, filtering for _ham0 motifs for search.
    Keep track of all ham variants for verification.

    Args:
        file_path (str): Path to the CSV file containing motif data

    Returns:
        tuple: (ham0_motifs list, motif_families dict)
    """
    print(f"Reading motifs: {file_path}")
    ham0_motifs = []
    all_motifs = []

    # Dictionary to group motif variants by family
    motif_families = defaultdict(list)

    with open(file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            motif_data = {
                'name': row['motif_name'],
                'start': int(row['start']),
                'sequence': row['sequence']
            }
            all_motifs.append(motif_data)

            # Extract family name (strip "_ham*" part)
            family_prefix = row['motif_name'].split('_ham')[0]
            motif_families[family_prefix].append(motif_data)

            # Filter to only include motifs with "_ham0" in their name
            # for search
            if "_ham0" in row['motif_name']:
                ham0_motifs.append(motif_data)

    print(f"    Loaded {len(all_motifs)} total motifs")
    print(f"    Filtered to {len(ham0_motifs)} ham0 motifs for search")
    print(f"    Found {len(motif_families)} motif families")

    return ham0_motifs, motif_families


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


def pigeon_hole_search(chromosome, motif, max_mismatches):
    """
    Implement approximate string matching using the pigeon-hole principle.

    Args:
        chromosome (str): DNA sequence to search in
        motif (str): DNA motif to search for
        max_mismatches (int): Maximum allowed mismatches

    Returns:
        tuple: (matches list, detailed_matches list)
    """
    matches = []
    detailed_matches = []  # (position, mismatch_count, mismatch_positions)
    n, m = len(chromosome), len(motif)

    # Apply pigeon-hole principle: divide motif into k+1 segments
    # If motif has ≤k mismatches, at least one segment must match exactly
    num_segments = max_mismatches + 1
    segment_length = m // num_segments

    # Create segments from the motif
    # Each segment gets a starting position and the actual sequence
    segments = []
    for i in range(num_segments):
        start = i * segment_length
        # Last segment gets any remaining characters if motif length doesn't
        # divide evenly
        end = start + segment_length if i < num_segments - 1 else m
        segments.append((start, motif[start:end]))

    # Store candidate positions where full motif might match
    # Using set to automatically handle duplicates from overlapping segments
    potential_matches = set()

    # Search for exact matches of each segment in the chromosome
    # Any segment match suggests the full motif might be nearby
    for segment_start, segment in segments:
        # Setup rightmost position of each character in segment
        bad_char = preprocess_bad_char(segment)

        i = 0
        while i <= n - len(segment):
            # Compare segment with chromosome using right-to-left approach
            j = len(segment) - 1
            match = True

            # Check if segment matches at current position
            while j >= 0:
                if segment[j] != chromosome[i + j]:
                    match = False
                    break
                j -= 1

            if match:
                # Found exact segment match - calculate where full motif
                # would start
                full_motif_start = i - segment_start

                # Only consider positions where full motif would fit in
                # chromosome
                if 0 <= full_motif_start <= n - m:
                    potential_matches.add(full_motif_start)
                i += 1
            else:
                # No match - use bad character rule to skip positions
                char = chromosome[i + j]
                i += max(1, j - bad_char.get(char, -1))

    # Verify each potential position by checking the complete motif
    # Allow up to max_mismatches differences
    for pos in potential_matches:
        mismatches = 0
        mismatch_positions = []  # Track where mismatches occur

        for j in range(m):
            # Check bounds and compare characters
            if pos + j < n and motif[j] != chromosome[pos + j]:
                mismatches += 1
                mismatch_positions.append(j)

                # Stop if we exceed allowed mismatches
                if mismatches > max_mismatches:
                    break

        # Accept this position if mismatches are within allowed limit
        if mismatches <= max_mismatches:
            matches.append(pos)
            detailed_matches.append((pos, mismatches, mismatch_positions))

    # Sort results by chromosome position for consistent output
    matches.sort()
    detailed_matches.sort()

    return matches, detailed_matches


def run_search(sequence, motifs, max_mismatches):
    """
    Run pigeon-hole search for each motif.

    Args:
        sequence (str): DNA sequence to search in
        motifs (list): List of motif dictionaries to search for
        max_mismatches (int): Maximum allowed mismatches

    Returns:
        dict: Dictionary containing search results for each motif
    """
    results = {}
    for motif in motifs:
        name = motif['name']
        motif_seq = motif['sequence']

        start_time = time.time()
        matches, detailed_matches = pigeon_hole_search(
            sequence, motif_seq, max_mismatches
        )
        execution_time = time.time() - start_time

        results[name] = {
            'sequence': motif_seq,
            'matches': matches,
            'detailed_matches': detailed_matches,
            'execution_time': execution_time,
            'max_mismatches': max_mismatches
        }

        print(f"Motif '{name}': Found {len(matches)} matches with "
              f"{max_mismatches} mismatches in {execution_time:.6f}s")

    return results


def verify_results(results, motif_families, max_mismatches):
    """
    Verify if results match expected positions with allowed mismatches.
    Check if ham variants are found.

    Args:
        results (dict): Search results from run_search function
        motif_families (dict): Dictionary of motif families
        max_mismatches (int): Maximum allowed mismatches

    Returns:
        dict: Dictionary containing verification results for each motif
    """
    verification = {}

    for ham0_name, search_data in results.items():
        # Extract family prefix (removing "_ham0")
        family_prefix = ham0_name.split('_ham')[0]

        # Get all family variants
        family_variants = motif_families.get(family_prefix, [])

        # Find positions from search
        found_positions = search_data['matches']

        # Initialize verification for this motif
        verification[ham0_name] = {
            'expected_position': None,
            'found_positions': found_positions,
            'correctly_found': False,
            'false_positives_count': 0,
            'ham_variants_found': {}
        }

        # Check each variant in the family
        for variant in family_variants:
            variant_name = variant['name']
            expected_pos = variant['start']

            # Check if this variant position was found in search results
            found = False
            closest_position = None
            closest_distance = float('inf')

            for pos in found_positions:
                distance = abs(pos - expected_pos)
                if distance <= max_mismatches:
                    found = True
                    if distance < closest_distance:
                        closest_position = pos
                        closest_distance = distance

            # For ham0 (base motif), set main verification fields
            if variant_name == ham0_name:
                verification[ham0_name]['expected_position'] = expected_pos
                verification[ham0_name]['correctly_found'] = found
                verification[ham0_name]['closest_position'] = closest_position

                # Count positions not matching the original ham0 as false
                # positives
                false_positive_count = 0

                for position in found_positions:
                    distance_from_expected = abs(position - expected_pos)

                    if distance_from_expected > max_mismatches:
                        false_positive_count += 1

                verification[ham0_name]['false_positives_count'] = (
                    false_positive_count
                )

            # For all variants, store in ham_variants_found
            ham_suffix = (variant_name.split('_ham')[1]
                          if '_ham' in variant_name else 'unknown')
            verification[ham0_name]['ham_variants_found'][
                f'ham{ham_suffix}'
            ] = {
                'name': variant_name,
                'position': expected_pos,
                'found': found,
                'closest_match': closest_position
            }

    return verification


def save_results(results, verification, output_file, max_mismatches):
    """Save results to CSV file using pandas."""
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Prepare data for DataFrame
    output_data = []
    for name, search_data in results.items():
        if name in verification:
            verify_data = verification[name]
            ham_variants = verify_data['ham_variants_found']

            # Prepare row data
            row_data = {
                'motif_name': name,
                'motif_sequence': search_data['sequence'],
                'expected_position': verify_data['expected_position'],
                'found_positions': ','.join(map(str, search_data['matches'])),
                'correctly_found': verify_data['correctly_found'],
                'closest_position': verify_data.get('closest_position', 'N/A'),
                'false_positives_count': verify_data['false_positives_count'],
                'execution_time': f"{search_data['execution_time']:.6f}",
                'max_mismatches': max_mismatches
            }

            # Add ham variant information
            for i in range(1, 6):  # Supports ham1 through ham5 ONLY
                ham_key = f'ham{i}'
                if ham_key in ham_variants:
                    row_data[f'{ham_key}_found'] = ham_variants[ham_key]['found']
                else:
                    row_data[f'{ham_key}_found'] = 'N/A'

            output_data.append(row_data)

    # Create DataFrame and save
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")


def main():
    """Main function to execute pigeon-hole DNA sequence matching."""
    args = parse_arguments()

    # Create output filename based on max_mismatches
    output_file = os.path.join(
        args.output_dir,
        f"pigeon_hole_search_results_max-mismatch-{args.max_mismatches}.csv"
    )

    # Read sequence and motifs
    sequence = read_fasta(args.input_fasta)
    ham0_motifs, motif_families = read_motifs(args.motif_positions)

    # Run search
    search_results = run_search(sequence, ham0_motifs, args.max_mismatches)

    # Verify results
    verification = verify_results(search_results, motif_families,
                                  args.max_mismatches)

    # Save results
    save_results(search_results, verification, output_file,
                 args.max_mismatches)

    print(f"Pigeon-hole search completed with {args.max_mismatches} "
          f"mismatches!")
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()