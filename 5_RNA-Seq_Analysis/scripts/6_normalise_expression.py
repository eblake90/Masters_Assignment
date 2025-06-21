#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd


def parse_arguments():
    """
    Parse command-line arguments for summary statistics generation.
    """
    parser = argparse.ArgumentParser(
        description='Calculate RPM and RPKM normalization from gene '
                    'coverage data'
    )
    parser.add_argument(
        '--coverage_file',
        default=('data/5_coverage/gene_coverage.txt'),
        help='Input file containing gene coverage data from featureCounts'
    )
    parser.add_argument(
        '--gtf_file',
        default=('data/4_annotations/chr21_22_annotation.gtf'),
        help='Input GTF file containing gene annotations'
    )
    parser.add_argument(
        '--output_dir',
        default=('data/6_normalized'),
        help='Output directory for normalized expression files'
    )
    return parser.parse_args()


def extract_gene_lengths(gtf_file):
    """
    Extract gene lengths from GTF file.

    Args:
        gtf_file (str): Path to GTF annotation file

    Returns:
        dict: Dictionary mapping gene IDs to their lengths in base pairs
    """
    print("  Parsing GTF file and extracting gene coordinates...")
    gene_lengths = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines starting with #
            if line.startswith('#'):
                continue

            # Split the line into tab-separated fields
            fields = line.strip().split('\t')

            # Process only gene features
            if len(fields) >= 9 and fields[2] == 'gene':
                # Extract gene_id from the attributes column using regex
                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)

                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    start_pos = int(fields[3])  # Start position
                    end_pos = int(fields[4])    # End position

                    # Calculate gene length
                    gene_length = end_pos - start_pos + 1
                    gene_lengths[gene_id] = gene_length

    print(f"  Found {len(gene_lengths)} genes in GTF file")
    return gene_lengths


def read_coverage_data(coverage_file):
    """
    Read gene coverage data from featureCounts output.

    Args:
        coverage_file (str): Path to featureCounts output file

    Returns:
        pd.DataFrame: DataFrame containing gene coverage data
    """
    # Read coverage data, skipping the first row which contains column headers
    coverage_df = pd.read_csv(
        coverage_file,
        sep='\t',
        skiprows=1,
        names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count']
    )

    # Convert Count column to numeric, handling any errors
    print("  Converting Count column to numeric...")
    coverage_df['Count'] = pd.to_numeric(coverage_df['Count'],
                                         errors='coerce')

    # Check for any NaN values after conversion
    nan_count = coverage_df['Count'].isna().sum()
    if nan_count > 0:
        print(f"  Warning: {nan_count} rows had non-numeric counts "
              f"and were set to NaN")
        # Remove rows with NaN counts
        coverage_df = coverage_df.dropna(subset=['Count'])

    print(f"  Loaded coverage data for {len(coverage_df)} genes")
    return coverage_df


def calculate_normalization(coverage_df, gene_lengths):
    """
    Calculate RPM and RPKM normalization values.

    Args:
        coverage_df (pd.DataFrame): DataFrame with gene coverage data
        gene_lengths (dict): Dictionary mapping gene IDs to lengths

    Returns:
        pd.DataFrame: DataFrame with added RPM and RPKM columns
    """
    print("  Calculating normalization values...")

    # Step 1: Add gene lengths to the coverage dataframe
    coverage_df['GTF_Length'] = coverage_df['Geneid'].map(gene_lengths)

    # Step 2: Calculate total reads mapped to genes for normalization
    # This is the scaling factor for RPM calculation
    total_reads = coverage_df['Count'].sum()
    print(f"  Total reads mapped to genes: {total_reads}")

    # Step 3: Calculate RPM (Reads Per Million)
    coverage_df['RPM'] = (coverage_df['Count'] / total_reads) * 1000000

    # Step 4: Calculate RPKM (Reads Per Kilobase per Million)
    coverage_df['Length_kb'] = coverage_df['GTF_Length'] / 1000
    coverage_df['RPKM'] = ((coverage_df['Count'] / coverage_df['Length_kb'])
                           / total_reads * 1000000)

    print("  RPM and RPKM calculations completed")
    return coverage_df


def create_output_summary(result_df):
    """
    Generate summary statistics and display results.
    """
    print("\n=== NORMALIZATION SUMMARY ===")
    print(f"Total genes analyzed: {len(result_df)}")
    print(f"Genes with reads > 0: {len(result_df[result_df['Raw_Count'] > 0])}")
    print(f"Mean RPM: {result_df['RPM'].mean():.2f}")
    print(f"Mean RPKM: {result_df['RPKM'].mean():.2f}")

    # Display top 5 most expressed genes by RPKM
    print(f"\nTop 5 genes by RPKM expression:")
    top_genes = result_df.nlargest(5, 'RPKM')[
        ['Gene_ID', 'Raw_Count', 'RPM', 'RPKM']
    ]
    print(top_genes.to_string(index=False))


def main():
    """
    Main function to orchestrate the normalization process.
    """
    args = parse_arguments()

    print("Starting gene expression normalization for Chr21 and Chr22...")
    print(f"Processing coverage file: {args.coverage_file}")
    print(f"Using GTF annotations: {args.gtf_file}")

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")

    # Step 1: Extract gene lengths from GTF file
    print("\nStep 1: Extracting gene lengths from GTF annotations")
    gene_lengths = extract_gene_lengths(args.gtf_file)

    # Step 2: Read gene coverage data
    print("\nStep 2: Reading gene coverage data")
    coverage_df = read_coverage_data(args.coverage_file)

    # Step 3: Calculate RPM and RPKM normalization
    print("\nStep 3: Calculating RPM and RPKM normalization")
    normalized_df = calculate_normalization(coverage_df, gene_lengths)

    # Step 4: Prepare final output dataframe
    print("\nStep 4: Preparing output files")
    result_df = normalized_df[
        ['Geneid', 'Chr', 'Count', 'GTF_Length', 'RPM', 'RPKM']
    ].copy()
    result_df.columns = [
        'Gene_ID', 'Chromosome', 'Raw_Count', 'Gene_Length_bp', 'RPM', 'RPKM'
    ]

    # Step 5: Save results to file
    output_file = os.path.join(args.output_dir, 'normalized_expression.txt')
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Results saved to: {output_file}")

    # Step 6: Generate and display summary
    print("\nStep 5: Generating summary statistics")
    create_output_summary(result_df)

    print(f"\nTask 5g completed successfully!")
    print(f"Normalized expression data available in: {args.output_dir}")


if __name__ == '__main__':
    main()