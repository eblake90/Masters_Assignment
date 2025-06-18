#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import math
import os


def parse_arguments():
    """
    Parse command-line arguments for summary statistics generation.
    """
    parser = argparse.ArgumentParser(
        description='Generate comprehensive summary statistics for '
                    'algorithm comparison results'
    )
    parser.add_argument(
        '--input_file',
        default="data/3_result_analysis/algorithm_comparison.csv",
        help='Algorithm comparison CSV file'
    )
    parser.add_argument(
        '--output_dir',
        default="data/3_result_analysis",
        help='Output directory for summary files'
    )
    parser.add_argument(
        '--max_mismatch_level',
        type=int,
        default=5,
        help='Maximum mismatch level used in comparison (default: 5)'
    )
    return parser.parse_args()


def extract_motif_length(motif_name):
    """
    Extract motif length from motif name using regex pattern.

    Args:
        motif_name (str): Name of the motif containing length information

    Returns:
        int or None: Extracted motif length, None if not found
    """
    match = re.search(r'_len(\d+)_', motif_name)
    if match:
        return int(match.group(1))
    return None


def calculate_overall_statistics(df, max_mismatch_level):
    """
    Calculate overall statistics for all algorithms.

    Args:
        df (pd.DataFrame): DataFrame containing algorithm comparison results
        max_mismatch_level (int): Maximum mismatch level to process

    Returns:
        list: List of dictionaries containing overall statistics for each algorithm
    """
    print("Calculating overall statistics...")

    # Define algorithm columns to analyze
    algorithms = ['naive', 'boyer_moore']

    # Add pigeon-hole algorithms for available mismatch levels
    for mm in range(max_mismatch_level + 1):
        ph_found_col = f'pigeon_hole_mm{mm}_found'
        if ph_found_col in df.columns:
            algorithms.append(f'pigeon_hole_mm{mm}')

    output_data = []
    naive_avg_time = None

    for alg in algorithms:
        found_col = f'{alg}_found'
        time_col = f'{alg}_time'

        # Direct calculations
        success_rate = df[found_col].mean() * 100
        times = df[time_col].astype(float)
        avg_time = times.mean()
        median_time = times.median()
        time_std_error = (times.std() / math.sqrt(len(times))
                          if len(times) > 1 else 0.0)

        # Store naive time for speedup calculations
        if alg == 'naive':
            naive_avg_time = avg_time

        # Calculate speedup ratio compared to naive algorithm
        if alg == 'naive':
            speedup_ratio = 1.0
        else:
            speedup_ratio = naive_avg_time / avg_time

        speedup_str = (f"{speedup_ratio:.3f}"
                       if speedup_ratio != 'N/A' else 'N/A')

        output_data.append({
            'algorithm': alg,
            'success_rate': f"{success_rate:.1f}%",
            'avg_execution_time': f"{avg_time:.6f}",
            'time_std_error': f"{time_std_error:.6f}",
            'median_execution_time': f"{median_time:.6f}",
            'speedup_ratio': speedup_str,
            'total_motifs': len(df)
        })

    return output_data


def calculate_ham_variant_statistics(df, max_mismatch_level):
    """
    Calculate ham variant statistics for pigeon-hole algorithms.

    Args:
        df (pd.DataFrame): DataFrame containing algorithm comparison results
        max_mismatch_level (int): Maximum mismatch level to process

    Returns:
        list: List of dictionaries containing ham variant statistics
    """
    print("Calculating ham variant statistics...")

    output_data = []

    for mm in range(max_mismatch_level + 1):
        alg_name = f'pigeon_hole_mm{mm}'
        row_data = {'algorithm': alg_name}

        # Check for ham variants (ham1 through ham5)
        for ham in range(1, 6):
            ham_col = f'pigeon_hole_mm{mm}_ham{ham}_found'

            if ham_col in df.columns:
                # Convert to boolean and calculate success rate
                ham_bool = df[ham_col].astype(str).str.lower() == 'true'
                success_rate = (ham_bool.sum() / len(df)) * 100
                row_data[f'ham{ham}_success_rate'] = f"{success_rate:.1f}%"
            else:
                row_data[f'ham{ham}_success_rate'] = "Column Missing"

        output_data.append(row_data)

    return output_data


def calculate_length_based_statistics(df, max_mismatch_level):
    """
    Calculate statistics grouped by motif length.

    Args:
        df (pd.DataFrame): DataFrame containing algorithm comparison results
        max_mismatch_level (int): Maximum mismatch level to process

    Returns:
        list: List of dictionaries containing length-based statistics
    """
    print("Calculating length-based statistics...")

    # Define algorithm columns to analyze
    algorithms = ['naive', 'boyer_moore']

    # Add pigeon-hole algorithms for available mismatch levels
    for mm in range(max_mismatch_level + 1):
        ph_found_col = f'pigeon_hole_mm{mm}_found'
        if ph_found_col in df.columns:
            algorithms.append(f'pigeon_hole_mm{mm}')

    output_data = []

    # Group by motif length
    for length in sorted(df['motif_length'].unique()):
        length_df = df[df['motif_length'] == length].copy()
        naive_avg_time = None

        for alg in algorithms:
            found_col = f'{alg}_found'
            time_col = f'{alg}_time'

            # Convert to string and check for 'True'
            success_rate = df[found_col].mean() * 100
            times = length_df[time_col].astype(float)
            avg_time = times.mean()
            median_time = times.median()
            time_std_error = (times.std() / math.sqrt(len(times))
                              if len(times) > 1 else 0.0)

            # Store naive time for speedup calculations
            if alg == 'naive':
                naive_avg_time = avg_time

            # Calculate speedup ratio compared to naive algorithm
            if alg == 'naive':
                speedup_ratio = 1.0
            else:
                speedup_ratio = naive_avg_time / avg_time

            output_data.append({
                'motif_length': length,
                'algorithm': alg,
                'success_rate': f"{success_rate:.1f}%",
                'avg_execution_time': f"{avg_time:.6f}",
                'time_std_error': f"{time_std_error:.6f}",
                'median_execution_time': f"{median_time:.6f}",
                'speedup_ratio': f"{speedup_ratio:.3f}",
                'total_motifs': len(length_df)
            })

    return output_data


def save_summary_statistics(overall_data, ham_data, length_data, output_dir):
    """
    Save all summary statistics to CSV files using pandas.
    """
    print("Saving summary statistics...")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save overall statistics
    overall_file = os.path.join(output_dir, "summary_overall.csv")
    overall_df = pd.DataFrame(overall_data)
    overall_df.to_csv(overall_file, index=False)
    print(f"    Overall statistics saved to: {overall_file}")

    # Save ham variant statistics
    ham_file = os.path.join(output_dir, "summary_ham_variants.csv")
    ham_df = pd.DataFrame(ham_data)
    ham_df.to_csv(ham_file, index=False)
    print(f"    Ham variant statistics saved to: {ham_file}")

    # Save length-based statistics
    length_file = os.path.join(output_dir, "summary_by_length.csv")
    length_df = pd.DataFrame(length_data)
    length_df.to_csv(length_file, index=False)
    print(f"    Length-based statistics saved to: {length_file}")


def main():
    """Main function to execute summary statistics generation."""
    print("Starting summary statistics generation...")

    # Parse command-line arguments
    args = parse_arguments()

    # Read comparison data
    print(f"Reading comparison data from: {args.input_file}")
    df = pd.read_csv(args.input_file)

    # Extract motif length for each row
    df['motif_length'] = df['motif_name'].apply(extract_motif_length)
    print(f"Loaded {len(df)} motif records")

    # Calculate statistics
    overall_data = calculate_overall_statistics(df, args.max_mismatch_level)
    ham_data = calculate_ham_variant_statistics(df, args.max_mismatch_level)
    length_data = calculate_length_based_statistics(df,
                                                    args.max_mismatch_level)

    # Save results
    save_summary_statistics(overall_data, ham_data, length_data,
                            args.output_dir)

    print(f"\nSummary statistics generation completed!")
    print(f"Output files saved to: {args.output_dir}")


if __name__ == "__main__":
    main()