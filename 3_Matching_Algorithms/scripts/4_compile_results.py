#!/usr/bin/env python3

import argparse
import pandas as pd
import os
from collections import defaultdict


def parse_arguments():
    """Parse command-line arguments for compiling algorithm results."""
    parser = argparse.ArgumentParser(
        description='Compile and compare results from DNA sequence matching '
                    'algorithms'
    )
    parser.add_argument(
        '--naive_results',
        default='data/2_matching_results/naive_search_results.csv',
        help='Path to naive search results CSV'
    )
    parser.add_argument(
        '--boyer_moore_results',
        default='data/2_matching_results/boyer_moore_search_results.csv',
        help='Path to Boyer-Moore search results CSV'
    )
    parser.add_argument(
        '--pigeon_hole_dir',
        default='data/2_matching_results',
        help='Directory containing pigeon-hole results files'
    )
    parser.add_argument(
        '--output_file',
        default='data/3_result_analysis/algorithm_comparison.csv',
        help='Output path for comparison CSV'
    )
    parser.add_argument(
        '--max_mismatch_level',
        type=int,
        default=5,
        help='Maximum mismatch level to include (default: 5)'
    )
    return parser.parse_args()


def read_basic_results(file_path, algorithm_name):
    """
    Read results from naive or Boyer-Moore CSV files.

    Args:
        file_path (str): Path to the CSV file containing algorithm results
        algorithm_name (str): Name of the algorithm for logging purposes

    Returns:
        dict: Dictionary mapping motif names to result dictionaries
    """
    print(f"Reading {algorithm_name} results from: {file_path}")

    df = pd.read_csv(file_path)
    results = {}

    for _, row in df.iterrows():
        motif_name = row['motif_name']
        correctly_found = str(row['correctly_found']).lower() == 'true'
        execution_time = float(row['execution_time'])
        results[motif_name] = {
            'found': correctly_found,
            'time': execution_time
        }

    print(f"    Loaded {len(results)} motif results for {algorithm_name}")
    return results


def read_pigeon_hole_results(pigeon_hole_dir, max_mismatch_level):
    """
    Read pigeon-hole results from multiple CSV files for different mismatch levels.

    Args:
        pigeon_hole_dir (str): Directory containing pigeon-hole result files
        max_mismatch_level (int): Maximum mismatch level to process

    Returns:
        tuple: (results dict, available_mismatch_levels list)
    """
    print(f"Reading pigeon-hole results from: {pigeon_hole_dir}")
    print(f"    Looking for mismatch levels 0 to {max_mismatch_level}")

    results = defaultdict(dict)
    available_mismatch_levels = []

    for mismatch_level in range(max_mismatch_level + 1):
        file_path = os.path.join(
            pigeon_hole_dir,
            f'pigeon_hole_search_results_max-mismatch-{mismatch_level}.csv'
        )
        df = pd.read_csv(file_path)
        available_mismatch_levels.append(mismatch_level)
        motif_count = 0

        for _, row in df.iterrows():
            motif_name = row['motif_name']
            correctly_found = str(row['correctly_found']).lower() == 'true'
            execution_time = float(row['execution_time'])

            # Extract ham variant information
            ham_variants = {}
            for ham_num in range(1, 6):  # ham1 through ham5
                ham_key = f'ham{ham_num}_found'
                if ham_key in row:
                    ham_value = str(row[ham_key])
                    if ham_value.lower() in ['true', 'false']:
                        ham_variants[f'ham{ham_num}'] = (
                                ham_value.lower() == 'true'
                        )

            results[motif_name][mismatch_level] = {
                'found': correctly_found,
                'time': execution_time,
                'ham_variants': ham_variants
            }
            motif_count += 1

        print(f"    Loaded {motif_count} motif results for mismatch level "
              f"{mismatch_level}")

    print(f"    Successfully loaded {len(available_mismatch_levels)} "
          f"mismatch levels: {available_mismatch_levels}")
    return dict(results), available_mismatch_levels


def create_comparison_table(naive_results, boyer_moore_results,
                            pigeon_hole_results, available_mismatch_levels):
    """
    Create a unified comparison table from all algorithm results.

    Args:
        naive_results (dict): Results from naive algorithm
        boyer_moore_results (dict): Results from Boyer-Moore algorithm
        pigeon_hole_results (dict): Results from pigeon-hole algorithm
        available_mismatch_levels (list): List of available mismatch levels

    Returns:
        list: List of dictionaries containing comparison data for each motif
    """
    print("Creating unified comparison table...")

    # Get all unique motif names
    all_motif_names = set()
    all_motif_names.update(naive_results.keys())
    all_motif_names.update(boyer_moore_results.keys())
    all_motif_names.update(pigeon_hole_results.keys())

    sorted_motif_names = sorted(all_motif_names)
    print(f"    Found {len(sorted_motif_names)} unique motifs")

    # Prepare data for DataFrame
    output_data = []

    for motif_name in sorted_motif_names:
        row_data = {'motif_name': motif_name}

        # Add naive results
        row_data['naive_found'] = naive_results[motif_name]['found']
        row_data['naive_time'] = f"{naive_results[motif_name]['time']:.3f}"

        # Add Boyer-Moore results
        row_data['boyer_moore_found'] = (
            boyer_moore_results[motif_name]['found']
        )
        row_data['boyer_moore_time'] = (
            f"{boyer_moore_results[motif_name]['time']:.3f}"
        )

        # Add pigeon-hole results for each mismatch level
        for mismatch_level in sorted(available_mismatch_levels):
            mm_found_field = f'pigeon_hole_mm{mismatch_level}_found'
            mm_time_field = f'pigeon_hole_mm{mismatch_level}_time'

            if (motif_name in pigeon_hole_results and
                    mismatch_level in pigeon_hole_results[motif_name]):

                ph_data = pigeon_hole_results[motif_name][mismatch_level]
                row_data[mm_found_field] = ph_data['found']
                row_data[mm_time_field] = f"{ph_data['time']:.3f}"

                # Add ham variant data
                ham_variants = ph_data['ham_variants']
                for ham_num in range(1, 6):
                    ham_field = (f'pigeon_hole_mm{mismatch_level}_ham'
                                 f'{ham_num}_found')
                    if f'ham{ham_num}' in ham_variants:
                        row_data[ham_field] = ham_variants[f'ham{ham_num}']

        output_data.append(row_data)

    print(f"    Created comparison data for {len(output_data)} motifs")
    return output_data


def save_comparison_table(output_data, output_file):
    """Save the comparison table to CSV file."""
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Create DataFrame and save
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, index=False)


def main():
    """Main function to execute algorithm results compilation."""

    # Parse command-line arguments
    args = parse_arguments()

    # Read naive search results
    naive_results = read_basic_results(args.naive_results, "naive")

    # Read Boyer-Moore search results
    boyer_moore_results = read_basic_results(args.boyer_moore_results,
                                             "Boyer-Moore")

    # Read pigeon-hole search results
    pigeon_hole_results, available_mismatch_levels = (
        read_pigeon_hole_results(args.pigeon_hole_dir,
                                 args.max_mismatch_level)
    )

    # Create comparison table
    output_data = create_comparison_table(
        naive_results, boyer_moore_results, pigeon_hole_results,
        available_mismatch_levels
    )

    # Save comparison table
    save_comparison_table(output_data, args.output_file)
    print("\nAlgorithm results compilation completed successfully!")
    print(f"Comparison table saved to: {args.output_file}")


if __name__ == "__main__":
    main()