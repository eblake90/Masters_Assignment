#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import re


def parse_arguments():
    """
    Parse command-line arguments for methylation analysis pipeline.
    """
    parser = argparse.ArgumentParser(
        description='Methylation analysis pipeline.'
    )
    parser.add_argument(
        '--input_file',
        default='data/Healthy_vs_Cancer_methyl_peptide_training_set.csv',
        help='Input CSV file'
    )
    parser.add_argument(
        '--output_dir',
        default='data',
        help='Output directory'
    )
    parser.add_argument(
        '--p_value_threshold',
        type=float,
        default=0.05,
        help='P-value threshold for statistical significance (default: 0.05)'
    )
    parser.add_argument(
        '--fold_change_threshold',
        type=float,
        default=1.5,
        help='Fold change threshold for differential methylation '
             '(default: 1.5)'
    )
    return parser.parse_args()


def extract_position_info(position_string):
    """
    Extract protein name and start/end positions from position string.

    Args:
        position_string (str): Position string from CSV data

    Returns:
        tuple or None: (protein, start, end) if valid, None otherwise
    """
    if pd.isna(position_string):
        return None

    pos_match = re.search(r'(\w+) \[(\d+)-(\d+)\]', position_string)
    if not pos_match:
        return None

    protein = pos_match.group(1)
    start = int(pos_match.group(2))
    end = int(pos_match.group(3))

    return protein, start, end


def extract_methylation_positions(modifications, sequence):
    """
    Extract methylation positions from modifications string.

    Args:
        modifications (str): Modifications string from CSV data
        sequence (str): Peptide sequence

    Returns:
        list: List of (amino_acid, position, methylation_type) tuples
    """
    if pd.isna(modifications) or pd.isna(sequence):
        return []

    positions = []
    pattern = r'(\d+)x(Methyl|Dimethyl|Trimethyl) \[([A-Z])(\d+)'
    matches = re.findall(pattern, modifications)

    for count, methyl_type, aa, pos in matches:
        positions.append((aa, int(pos), methyl_type))

    return positions


def is_kp_rp_motif(aa, pos, sequence):
    """
    Check if amino acid is followed by proline (KP/RP motif).

    Args:
        aa (str): Amino acid character
        pos (int): Position in sequence
        sequence (str): Peptide sequence

    Returns:
        bool: True if amino acid is followed by proline
    """
    if aa not in ['K', 'R']:
        return False

    # Remove special characters from sequence
    clean_seq = re.sub(r'[\[\].]', '', sequence)

    if pos + 1 < len(clean_seq) and clean_seq[pos + 1] == 'P':
        return True

    return False


def find_cancer_methylated_peptides(data):
    """
    Find methylated peptides present in cancer samples.

    Args:
        data (pd.DataFrame): Input peptide data

    Returns:
        pd.DataFrame: Filtered methylated peptides from cancer samples
    """
    # Cancer sample columns
    cancer_cols = [
        "Found in Sample: [S4] F4: Cancer, 3",
        "Found in Sample: [S5] F5: Cancer, 2",
        "Found in Sample: [S6] F6: Cancer, 1"
    ]

    # Filter for peptides found in cancer samples
    cancer_mask = ~data[cancer_cols].isin(["", "Not Found"]).all(axis=1)
    cancer_peptides = data[cancer_mask].copy()

    # Filter for methylated peptides
    methylation_pattern = r'\d+x(?:Methyl|Dimethyl|Trimethyl)'
    has_methylation = cancer_peptides['Modifications'].str.contains(
        methylation_pattern, na=False, regex=True
    )

    # Filter for miscleavage > 0 (indicates blocked cleavage)
    miscleavage_filter = cancer_peptides['# Missed Cleavages'] > 0

    # Combine filters
    methylated = cancer_peptides[has_methylation & miscleavage_filter].copy()

    # Remove terminal methylations (not relevant for cleavage blocking)
    valid_indices = []
    for idx, row in methylated.iterrows():
        if not is_terminal_methylation(row['Modifications'],
                                       row['Annotated Sequence']):
            valid_indices.append(idx)

    return methylated.loc[valid_indices].copy()


def is_terminal_methylation(modifications, sequence):
    """
    Check if methylation is at the terminal position.

    Args:
        modifications (str): Modifications string from CSV data
        sequence (str): Peptide sequence

    Returns:
        bool: True if methylation is at terminal position
    """
    if pd.isna(modifications) or pd.isna(sequence):
        return False

    # Extract methylation positions
    pattern = r'\d+x(?:Methyl|Dimethyl|Trimethyl) \[(\w)(\d+)'
    matches = re.findall(pattern, modifications)

    if not matches:
        return False

    # Get clean sequence length
    clean_seq = re.sub(r'[\[\].]', '', sequence)
    seq_length = len(clean_seq)

    # Check if any methylation is at last position
    for aa, pos in matches:
        if int(pos) == seq_length:
            return True

    return False


def build_unmethylated_dict(all_data):
    """
    Build dictionary of unmethylated peptides keyed by position.

    Args:
        all_data (pd.DataFrame): Complete peptide dataset

    Returns:
        dict: Dictionary keyed by (protein, start, end) tuples
    """
    unmethylated_dict = {}
    methylation_pattern = r'\d+x(?:Methyl|Dimethyl|Trimethyl)'

    for idx, row in all_data.iterrows():
        # Skip peptides with methylation modifications
        if (pd.notna(row['Modifications']) and
            re.search(methylation_pattern, row['Modifications'])):
            continue

        # Extract position information
        position_info = extract_position_info(
            row['Positions in Master Proteins']
        )
        if position_info:
            protein, start, end = position_info
            key = (protein, start, end)
            unmethylated_dict[key] = row

    return unmethylated_dict


def fragment_exists(protein, start, end, unmethylated_dict):
    """
    Check if a specific fragment exists in the unmethylated dictionary.

    Args:
        protein (str): Protein name
        start (int): Fragment start position
        end (int): Fragment end position
        unmethylated_dict (dict): Dictionary of unmethylated peptides

    Returns:
        bool: True if fragment exists in dictionary
    """
    fragment_key = (protein, start, end)
    return fragment_key in unmethylated_dict


def is_peptide_valid(peptide_row, unmethylated_dict):
    """
    Check if methylated peptide is valid by finding expected cleavage fragments.

    Args:
        peptide_row (pd.Series): Row data for methylated peptide
        unmethylated_dict (dict): Dictionary of unmethylated peptides

    Returns:
        bool: True if expected cleavage fragments are found
    """
    # Extract peptide position
    position_info = extract_position_info(
        peptide_row['Positions in Master Proteins']
    )
    if not position_info:
        return False

    protein, pep_start, pep_end = position_info

    # Extract methylation positions
    methyl_positions = extract_methylation_positions(
        peptide_row['Modifications'],
        peptide_row['Annotated Sequence']
    )

    if not methyl_positions:
        return False

    # Check each methylation site for validation
    for aa, methyl_pos, methyl_type in methyl_positions:
        # Skip KP/RP motifs (enzyme can't cleave anyway)
        if is_kp_rp_motif(aa, methyl_pos,
                          peptide_row['Annotated Sequence']):
            continue

        # Calculate expected fragment positions if cleavage occurred
        left_fragment_end = pep_start + methyl_pos - 1
        right_fragment_start = pep_start + methyl_pos

        # Check if either expected fragment exists
        if fragment_exists(protein, pep_start, left_fragment_end,
                          unmethylated_dict):
            return True

        if fragment_exists(protein, right_fragment_start, pep_end,
                          unmethylated_dict):
            return True

    return False


def validate_peptides(methylated_peptides, all_data):
    """
    Validate methylated peptides by finding corresponding unmethylated fragments.

    Args:
        methylated_peptides (pd.DataFrame): Methylated peptides to validate
        all_data (pd.DataFrame): Complete dataset for finding controls

    Returns:
        pd.DataFrame: Validated methylated peptides
    """
    # Build dictionary of all unmethylated peptides
    print("Building unmethylated peptide dictionary...")
    unmethylated_dict = build_unmethylated_dict(all_data)
    print(f"Found {len(unmethylated_dict)} unmethylated peptides")

    # Validate each methylated peptide
    validated_indices = []

    for idx, peptide_row in methylated_peptides.iterrows():
        if is_peptide_valid(peptide_row, unmethylated_dict):
            validated_indices.append(idx)

    print(f"Validated {len(validated_indices)} out of "
          f"{len(methylated_peptides)} methylated peptides")
    return methylated_peptides.loc[validated_indices].copy()


def filter_significant_peptides(peptides, p_value_threshold,
                               fold_change_threshold):
    """
    Filter peptides based on statistical significance criteria.

    Args:
        peptides (pd.DataFrame): Input peptides to filter
        p_value_threshold (float): P-value threshold for significance
        fold_change_threshold (float): Fold change threshold

    Returns:
        pd.DataFrame: Significantly changed peptides
    """
    # Check if required columns exist
    ratio_col = 'Abundance Ratio: (Cancer) / (Healthy)'
    pvalue_col = 'Abundance Ratio Adj. P-Value: (Cancer) / (Healthy)'
    quan_col = 'Quan Info'

    # Filter out peptides with no quantification values
    filtered_peptides = peptides.copy()
    if quan_col in filtered_peptides.columns:
        no_quan_mask = filtered_peptides[quan_col] == 'No Quan Values'
        if no_quan_mask.any():
            print(f"Filtering out {no_quan_mask.sum()} peptides with "
                  f"'No Quan Values'")
            filtered_peptides = filtered_peptides[~no_quan_mask].copy()

    # Apply significance criteria
    significant_indices = []

    for idx, row in filtered_peptides.iterrows():
        if pd.notna(row[ratio_col]) and pd.notna(row[pvalue_col]):
            ratio = row[ratio_col]
            pvalue = row[pvalue_col]

            # Check if meets significance criteria
            if (pvalue <= p_value_threshold and
                (ratio >= fold_change_threshold or
                 ratio <= 1 / fold_change_threshold)):
                significant_indices.append(idx)

    significant_peptides = filtered_peptides.loc[significant_indices].copy()
    print(f"Found {len(significant_peptides)} significant peptides out of "
          f"{len(filtered_peptides)} total")

    return significant_peptides


def create_output(peptides, all_data, output_dir, p_value_threshold,
                 fold_change_threshold):
    """
    Create output files with control and cancer sequences.

    Args:
        peptides (pd.DataFrame): Validated peptides to output
        all_data (pd.DataFrame): Complete dataset for finding controls
        output_dir (str): Directory path for output files
        p_value_threshold (float): P-value threshold for significance
        fold_change_threshold (float): Fold change threshold for significance

    Returns:
        None
    """
    # Build unmethylated peptides dictionary for finding control sequences
    unmethylated_dict = {}
    methylation_pattern = r'\d+x(?:Methyl|Dimethyl|Trimethyl)'

    for idx, row in all_data.iterrows():
        if (pd.notna(row['Modifications']) and
            re.search(methylation_pattern, row['Modifications'])):
            continue

        pos_match = re.search(r'(\w+) \[(\d+)-(\d+)\]',
                             row['Positions in Master Proteins'])
        if pos_match:
            protein = pos_match.group(1)
            start = int(pos_match.group(2))
            end = int(pos_match.group(3))
            key = (protein, start, end)
            unmethylated_dict[key] = row

    # Prepare output data
    output_data = []

    for idx, row in peptides.iterrows():
        # Extract sample information
        healthy_samples = []
        cancer_samples = []

        healthy_cols = [
            "Found in Sample: [S1] F1: Healthy, 3",
            "Found in Sample: [S2] F2: Healthy, 2",
            "Found in Sample: [S3] F3: Healthy, 1"
        ]
        cancer_cols = [
            "Found in Sample: [S4] F4: Cancer, 3",
            "Found in Sample: [S5] F5: Cancer, 2",
            "Found in Sample: [S6] F6: Cancer, 1"
        ]

        for col in healthy_cols:
            if row.get(col, "Not Found") not in ["", "Not Found"]:
                sample_id = re.search(r'\[(\w+)\]', col).group(1)
                healthy_samples.append(sample_id)

        for col in cancer_cols:
            if row.get(col, "Not Found") not in ["", "Not Found"]:
                sample_id = re.search(r'\[(\w+)\]', col).group(1)
                cancer_samples.append(sample_id)

        # Get cancer sequence (methylated) with highlighting
        cancer_sequence = re.sub(r'[\[\].]', '', row['Annotated Sequence'])

        # Highlight methylated amino acids with modification type
        if pd.notna(row['Modifications']):
            # Extract methylation positions and types
            methyl_pattern = (r'(\d+)x(Methyl|Dimethyl|Trimethyl) '
                             r'\[([A-Z])(\d+)')
            methyl_matches = re.findall(methyl_pattern, row['Modifications'])

            # Map methylation types to letters
            methyl_map = {'Methyl': 'M', 'Dimethyl': 'D', 'Trimethyl': 'T'}

            seq_list = list(cancer_sequence)

            # Replace each methylated amino acid
            for count, methyl_type, aa, pos in methyl_matches:
                list_index = int(pos)  # Convert to 0-indexed
                if 0 <= list_index < len(seq_list):
                    methyl_letter = methyl_map.get(methyl_type)
                    seq_list[list_index] = (f"{{{methyl_letter}-"
                                          f"{seq_list[list_index]}}}")

            cancer_sequence = ''.join(seq_list)

        # Add position brackets
        cancer_pos = row['Positions in Master Proteins']
        cancer_pos_match = re.search(r'\[(\d+)-(\d+)\]', cancer_pos)
        if cancer_pos_match:
            cancer_start = cancer_pos_match.group(1)
            cancer_end = cancer_pos_match.group(2)
            cancer_sequence_with_pos = (f"[{cancer_start}]-{cancer_sequence}"
                                      f"-[{cancer_end}]")
        else:
            cancer_sequence_with_pos = cancer_sequence

        # Find overlapping control sequences (unmethylated)
        pos_match = re.search(r'(\w+) \[(\d+)-(\d+)\]',
                             row['Positions in Master Proteins'])
        control_sequences = []
        overlapping_positions = []

        if pos_match:
            protein = pos_match.group(1)
            methyl_start = int(pos_match.group(2))
            methyl_end = int(pos_match.group(3))

            # Find all unmethylated peptides that overlap with methylated
            # peptide position
            for (prot, start, end), unmeth_row in unmethylated_dict.items():
                if prot == protein:
                    # Check if positions overlap
                    if not (end < methyl_start or start > methyl_end):
                        frag_seq = re.sub(r'[\[\].]', '',
                                        unmeth_row['Annotated Sequence'])
                        control_sequences.append(frag_seq)
                        overlapping_positions.extend([start, end])

            # Stitch sequences together with position range
            if control_sequences and overlapping_positions:
                min_pos = min(overlapping_positions)
                max_pos = max(overlapping_positions)
                stitched_sequence = " | ".join(control_sequences)
                control_sequence = (f"[{min_pos}]-{stitched_sequence}"
                                  f"-[{max_pos}]")
            else:
                control_sequence = "No overlapping fragments found"
        else:
            control_sequence = "No overlapping fragments found"

        # Handle abundance values
        abundance_healthy = row.get('Abundances (Grouped): Healthy', '')
        abundance_cancer = row.get('Abundances (Grouped): Cancer', '')

        # Check if values are missing, empty, or NaN
        if pd.isna(abundance_healthy) or abundance_healthy == '':
            abundance_healthy = "None"

        if pd.isna(abundance_cancer) or abundance_cancer == '':
            abundance_cancer = "None"

        # Create output row
        output_row = {
            'Control_Sequence': control_sequence,
            'Cancer_Sequence': cancer_sequence_with_pos,
            'Modifications': row['Modifications'],
            'Protein_Position': row['Positions in Master Proteins'],
            'Miscleavages': row['# Missed Cleavages'],
            'Healthy_Samples': ("; ".join(healthy_samples)
                              if healthy_samples else "None"),
            'Cancer_Samples': ("; ".join(cancer_samples)
                             if cancer_samples else "None"),
            'Abundance_Healthy': abundance_healthy,
            'Abundance_Cancer': abundance_cancer,
            'PSMs': row['# PSMs']
        }

        output_data.append(output_row)

    # Save all results
    output_df = pd.DataFrame(output_data)
    output_file = os.path.join(output_dir, "methylation_results.csv")
    output_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")

    # Filter and save significant results
    significant_peptides = filter_significant_peptides(
        peptides, p_value_threshold, fold_change_threshold
    )

    if len(significant_peptides) > 0:
        # Create output data for significant peptides only
        significant_output_data = []
        significant_indices = significant_peptides.index

        # Check if required columns exist for categorization
        ratio_col = 'Abundance Ratio: (Cancer) / (Healthy)'
        pvalue_col = 'Abundance Ratio Adj. P-Value: (Cancer) / (Healthy)'

        for i, output_row in enumerate(output_data):
            peptide_idx = peptides.index[i]
            if peptide_idx in significant_indices:
                # Add methylation category for significant peptides
                peptide_row = peptides.loc[peptide_idx]
                if (pd.notna(peptide_row[ratio_col]) and
                    pd.notna(peptide_row[pvalue_col])):
                    ratio = peptide_row[ratio_col]
                    pvalue = peptide_row[pvalue_col]

                    if pvalue <= p_value_threshold:
                        if ratio >= fold_change_threshold:
                            category = "Cancer_Enriched"
                        elif ratio <= 1 / fold_change_threshold:
                            category = "Healthy_Enriched"
                        else:
                            category = "Non_differential"
                    else:
                        category = "Non_differential"
                else:
                    category = "Unknown"

                # Add category to output row
                output_row_with_category = output_row.copy()
                output_row_with_category['Methylation_Category'] = category
                significant_output_data.append(output_row_with_category)

        significant_output_df = pd.DataFrame(significant_output_data)
        significant_output_file = os.path.join(
            output_dir, "methylation_results_significant.csv"
        )
        significant_output_df.to_csv(significant_output_file, index=False)
        print(f"Significant results saved to: {significant_output_file}")
    else:
        print("No significant peptides found - "
              "methylation_results_significant.csv not created")

    # Create summary
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"Methylation Analysis Summary\n")
        f.write(f"Total validated peptides: {len(peptides)}\n")
        f.write(f"Significant peptides (p<{p_value_threshold}, "
                f"FC>{fold_change_threshold}): "
                f"{len(significant_peptides) if len(significant_peptides) > 0 else 0}\n")

        # Add category breakdown if significant peptides exist
        if len(significant_peptides) > 0 and 'significant_output_df' in locals():
            f.write(f"\nSignificant Peptide Categories:\n")
            category_counts = significant_output_df['Methylation_Category'].value_counts()
            for category, count in category_counts.items():
                f.write(f"  {category}: {count}\n")

    print(f"Summary saved to: {summary_file}")


def main():
    """
    Main function to execute methylation analysis pipeline.
    """
    args = parse_arguments()

    input_file = args.input_file
    output_dir = args.output_dir
    p_value_threshold = args.p_value_threshold
    fold_change_threshold = args.fold_change_threshold

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    print(f"Loading data from: {input_file}")
    data = pd.read_csv(input_file)
    print(f"Loaded {len(data)} peptides")

    # Find methylated peptides in cancer samples
    print("Finding methylated peptides in cancer samples...")
    methylated_peptides = find_cancer_methylated_peptides(data)
    print(f"Found {len(methylated_peptides)} methylated peptides")

    # Validate peptides by finding expected cleavage fragments
    print("Validating peptides...")
    validated_peptides = validate_peptides(methylated_peptides, data)
    print(f"Validated {len(validated_peptides)} peptides")

    # Create output files
    print("Creating output files...")
    create_output(validated_peptides, data, output_dir, p_value_threshold,
                 fold_change_threshold)

    print("Analysis complete!")


if __name__ == "__main__":
    main()