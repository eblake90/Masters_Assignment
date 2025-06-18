#!/usr/bin/env python3

import argparse
import csv
import os
import re
import pandas as pd


def parse_arguments():
    """
    Parse command-line arguments for extracting protein accessions.
    """
    parser = argparse.ArgumentParser(
        description='Extract protein accessions from methylation data.'
    )
    parser.add_argument(
        '--input_file',
        default='../1_Proteomics_Methylation_Analysis/'
                'data/methylation_results_significant.csv',
        help='Input CSV file with methylation data'
    )
    parser.add_argument(
        '--output_file',
        default='data/1_accessions_input.csv',
        help='Output CSV file for accessions'
    )
    return parser.parse_args()


def extract_accessions_from_dataframe(df):
    """
    Extract unique protein accessions from methylation dataframe.

    Args:
        df (pd.DataFrame): DataFrame containing methylation data with
                          'Protein_Position' column

    Returns:
        set: Set of unique protein accessions found in the data
    """
    accessions = set()

    # Iterate through each protein position entry
    for protein_position in df['Protein_Position']:
        # Split by semicolon in case multiple positions are listed
        for pos in str(protein_position).split(';'):
            # Extract 6-character accession ID from start of position string
            # Pattern matches protein IDs like 'P12345' or 'Q9Y123'
            match = re.match(r'^([A-Z0-9]{6})', pos.strip())
            if match:
                accessions.add(match.group(1))

    return accessions


def save_accessions_to_csv(accessions, output_file):
    """
    Save protein accessions to CSV file.
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Write accessions to CSV file with header
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Accession'])

        # Sort accessions for consistent output
        for accession in sorted(accessions):
            writer.writerow([accession])


def main():
    """
    Main function to execute protein accession extraction pipeline.
    """
    # Parse command-line arguments
    args = parse_arguments()

    # Read methylation data from CSV
    print(f"Reading data from: {args.input_file}")
    df = pd.read_csv(args.input_file)
    print(f"Loaded {len(df)} rows")

    # Extract unique protein accessions
    accessions = extract_accessions_from_dataframe(df)
    print(f"Found {len(accessions)} unique accessions")

    # Save accessions to output CSV
    save_accessions_to_csv(accessions, args.output_file)
    print(f"Accessions saved to: {args.output_file}")


if __name__ == "__main__":
    main()