#!/usr/bin/env python3

import argparse
import math
import os
import re


def parse_arguments():
    """
    Parse command-line arguments for mini-BLAST alignment.
    """
    parser = argparse.ArgumentParser(
        description='perform mini-BLAST alignment between protein sequences'
    )
    parser.add_argument(
        '--query_file',
        default=('data/input/query.fa'),
        help='input fasta file containing query protein sequence'
    )
    parser.add_argument(
        '--target_file',
        default=('data/input/target.fa'),
        help='input fasta file containing target protein sequences'
    )
    parser.add_argument(
        '--output_file',
        default=('data/output/mini_blast_results.txt'),
        help='output file path for alignment results'
    )
    parser.add_argument(
        '--threshold_score',
        default=11,
        type=int,
        help='threshold score for extending hits'
    )
    parser.add_argument(
        '--neighborhood_threshold',
        default=13,
        type=int,
        help='threshold for neighborhood word scoring'
    )
    return parser.parse_args()


def create_blosum62():
    """
    Create BLOSUM62 substitution matrix based on Henikoff & Henikoff 1992.

    Args:
        None

    Returns:
        dict: Dictionary mapping amino acid pairs to substitution scores
    """
    blosum62 = {}
    amino_acids = ['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q',
                   'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

    # BLOSUM62 scores matrix - symmetric matrix stored as upper triangle
    scores = [
        [9],  # C - Cysteine
        [-1, 4],  # S - Serine
        [-1, 1, 5],  # T - Threonine
        [-3, -1, -1, 7],  # P - Proline
        [0, 1, 0, -1, 4],  # A - Alanine
        [-3, 0, -2, -2, 0, 6],  # G - Glycine
        [-3, 1, 0, -2, -2, 0, 6],  # N - Asparagine
        [-3, 0, -1, -1, -2, -1, 1, 6],  # D - Aspartic Acid
        [-4, 0, -1, -1, -1, -2, 0, 2, 5],  # E - Glutamic Acid
        [-3, 0, -1, -1, -1, -2, 0, 0, 2, 5],  # Q - Glutamine
        [-3, -1, -2, -2, -2, -2, 1, -1, 0, 0, 8],  # H - Histidine
        [-3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5],  # R - Arginine
        [-3, 0, -1, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5],  # K - Lysine
        [-1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5],  # M
        [-1, -2, -1, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4],  # I
        [-1, -2, -1, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4],  # L
        [-1, -2, 0, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4],  # V
        [-2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6],  # F
        [-2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7],  # Y
        [-2, -3, -2, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2, 11]  # W
    ]

    # Build symmetric matrix from upper triangle
    for i, aa1 in enumerate(amino_acids):
        for j, aa2 in enumerate(amino_acids[:i + 1]):
            score = scores[i][j]
            blosum62[(aa1, aa2)] = score
            blosum62[(aa2, aa1)] = score  # Ensure matrix is symmetric

    return blosum62


def parse_fasta(filename):
    """
    Read FASTA file and return dictionary of sequences.

    Args:
        filename (str): Path to FASTA file to parse

    Returns:
        dict: Dictionary mapping headers to sequences
    """
    sequences = {}
    current_header = ''
    current_sequence = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence if it exists
                if current_header:
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

    # Save the last sequence
    if current_header:
        sequences[current_header] = ''.join(current_sequence)

    return sequences


def filter_low_complexity(sequence):
    """
    Remove low complexity regions (runs of 4+ identical amino acids).

    Args:
        sequence (str): Protein sequence to filter

    Returns:
        str: Filtered sequence with low complexity regions removed
    """
    # Regex pattern matches runs of 4 or more identical amino acids
    # (.) captures any single character
    # \1{3,} matches the same character 3 or more additional times
    low_complexity_pattern = r'(.)\1{3,}'

    # Remove matched patterns by replacing with empty string
    filtered_sequence = re.sub(low_complexity_pattern, '', sequence)

    return filtered_sequence


def generate_words(sequence, word_size):
    """
    Generate overlapping words (k-mers) from sequence.

    Args:
        sequence (str): Input protein sequence
        word_size (int): Length of words to generate

    Returns:
        tuple: (words list, positions dict mapping words to positions)
    """
    words = []
    positions = {}  # Track positions of each word in sequence

    # Generate words using sliding window approach
    for i in range(len(sequence) - word_size + 1):
        word = sequence[i:i + word_size]
        words.append(word)

        # Track positions for later hit extension
        if word not in positions:
            positions[word] = []
        positions[word].append(i)

    return words, positions


def score_word_pair(word1, word2, blosum62):
    """
    Calculate alignment score between two words using BLOSUM62.

    Args:
        word1 (str): First word to compare
        word2 (str): Second word to compare
        blosum62 (dict): BLOSUM62 substitution matrix

    Returns:
        int: Total alignment score for word pair
    """
    score = 0
    for i in range(len(word1)):
        # Get substitution score, default to -4 if pair not found
        score += blosum62.get((word1[i], word2[i]), -4)
    return score


def generate_neighborhood_words(words, blosum62, threshold):
    """
    Generate high-scoring neighborhood words for each query word.

    Args:
        words (list): List of query words
        blosum62 (dict): BLOSUM62 substitution matrix
        threshold (int): Minimum score threshold for neighborhood words

    Returns:
        dict: Dictionary mapping each word to list of similar words
    """
    # Get all amino acids from BLOSUM62 matrix
    amino_acids = set()
    for pair in blosum62.keys():
        if isinstance(pair, tuple) and len(pair) == 2:
            amino_acids.add(pair[0])
            amino_acids.add(pair[1])
    amino_acids = list(amino_acids)

    expanded_words = {}

    # Process each unique word (avoid duplicates)
    for word in set(words):
        expanded_words[word] = [word]  # Always include original word

        # Generate all possible 3-word combinations
        for a1 in amino_acids:
            for a2 in amino_acids:
                for a3 in amino_acids:
                    candidate = a1 + a2 + a3
                    if candidate == word:
                        continue  # Skip original word

                    # Calculate similarity score using BLOSUM62
                    score = score_word_pair(word, candidate, blosum62)

                    # Keep candidate if score meets threshold
                    if score >= threshold:
                        expanded_words[word].append(candidate)

    return expanded_words


def find_word_matches(target_seq, expanded_words, word_size, positions):
    """
    Find exact matches between target sequence and expanded word list.

    Args:
        target_seq (str): Target sequence to search
        expanded_words (dict): Dictionary of query words and similar words
        word_size (int): Length of words
        positions (dict): Dictionary mapping words to their positions

    Returns:
        list: List of (query_position, target_position) tuples
    """
    matches = []

    # Scan through target sequence
    for i in range(len(target_seq) - word_size + 1):
        target_word = target_seq[i:i + word_size]

        # Check if target word matches any expanded query word
        for query_word, similar_words in expanded_words.items():
            if target_word in similar_words:
                # Add all original positions of matching query word
                for q_pos in positions[query_word]:
                    matches.append((q_pos, i))

    return matches


def extend_hit(query, target, q_pos, t_pos, blosum62, word_size,
               threshold_score=11):
    """
    Extend a seed hit in both directions to form High-scoring Segment Pair.

    Args:
        query (str): Query sequence
        target (str): Target sequence
        q_pos (int): Starting position in query
        t_pos (int): Starting position in target
        blosum62 (dict): BLOSUM62 substitution matrix
        word_size (int): Length of seed word
        threshold_score (int): Threshold for extension termination

    Returns:
        dict: HSP information including positions, score, and alignments
    """
    # Initialize with seed match boundaries
    left_q = q_pos
    right_q = q_pos + word_size - 1
    left_t = t_pos
    right_t = t_pos + word_size - 1

    # Calculate initial score for seed word
    current_score = 0
    for i in range(word_size):
        aa1 = query[q_pos + i]
        aa2 = target[t_pos + i]
        current_score += blosum62.get((aa1, aa2), -4)

    best_score = current_score
    best_endpoints = (left_q, right_q, left_t, right_t)

    # Extend to the right without gaps
    max_score_right = current_score
    while right_q < len(query) - 1 and right_t < len(target) - 1:
        right_q += 1
        right_t += 1

        # Score this amino acid pair
        aa1 = query[right_q]
        aa2 = target[right_t]
        pair_score = blosum62.get((aa1, aa2), -4)
        current_score += pair_score

        # Track maximum score during extension
        if current_score > max_score_right:
            max_score_right = current_score

        # Update the best alignment if score improved
        if current_score > best_score:
            best_score = current_score
            best_endpoints = (left_q, right_q, left_t, right_t)

        # Stop if score drops too much from maximum (X-drop criterion)
        if max_score_right - current_score > threshold_score:
            break

    # Reset score for left extension
    current_score = best_score
    best_endpoints = (left_q, right_q, left_t, right_t)

    # Extend to the left without gaps
    max_score_left = current_score
    while left_q > 0 and left_t > 0:
        left_q -= 1
        left_t -= 1

        # Score this amino acid pair
        aa1 = query[left_q]
        aa2 = target[left_t]
        pair_score = blosum62.get((aa1, aa2), -4)
        current_score += pair_score

        # Track maximum score during extension
        if current_score > max_score_left:
            max_score_left = current_score

        # Update the best alignment if score improved
        if current_score > best_score:
            best_score = current_score
            best_endpoints = (left_q, right_q, left_t, right_t)

        # Stop if score drops too much from maximum
        if max_score_left - current_score > threshold_score:
            break

    # Return the best HSP found
    return {
        'query_start': best_endpoints[0],
        'query_end': best_endpoints[1],
        'target_start': best_endpoints[2],
        'target_end': best_endpoints[3],
        'score': best_score,
        'q_align': query[best_endpoints[0]:best_endpoints[1] + 1],
        't_align': target[best_endpoints[2]:best_endpoints[3] + 1],
        'length': best_endpoints[1] - best_endpoints[0] + 1
    }


def calculate_e_value(score, query_length, db_length, k=0.134,
                      lambda_val=0.3176):
    """
    Calculate E-value for statistical significance assessment.

    Args:
        score (int): Raw alignment score
        query_length (int): Length of query sequence
        db_length (int): Length of database sequence
        k (float): Karlin-Altschul parameter K
        lambda_val (float): Karlin-Altschul parameter lambda

    Returns:
        float: E-value for the alignment score
    """
    # Use parameters from BLAST paper for ungapped protein alignments
    return k * query_length * db_length * math.exp(-lambda_val * score)


def run_mini_blast(query, targets, word_size=3, threshold_score=11,
                   neighborhood_threshold=13):
    """
    Execute mini-BLAST algorithm to compare query against target sequences.

    Args:
        query (str): Query protein sequence
        targets (dict): Dictionary of target sequences
        word_size (int): Size of initial seed words
        threshold_score (int): Minimum score for significant HSPs
        neighborhood_threshold (int): Threshold for neighborhood word scoring

    Returns:
        list: List of (target_header, best_hsp) tuples ranked by E-value
    """
    print("Running mini-BLAST algorithm...")

    # Step 1: Create BLOSUM62 scoring matrix
    blosum62 = create_blosum62()

    # Step 2: Filter low complexity regions from query
    print("Step 1: Filtering low complexity regions")
    filtered_query = filter_low_complexity(query)

    # Step 3: Generate overlapping words from filtered query
    print("Step 2: Generating words from query")
    query_words, positions = generate_words(filtered_query, word_size)

    # Step 4: Generate neighborhood words with similar scoring
    print(f"Step 3: Generating neighborhood words "
          f"(threshold: {neighborhood_threshold})")
    expanded_words = generate_neighborhood_words(query_words, blosum62,
                                                 neighborhood_threshold)

    # Process each target sequence
    results = {}

    for target_header, target_seq in targets.items():
        print(f"Processing target: {target_header}")
        results[target_header] = []

        # Step 5: Find word matches between query and target
        matches = find_word_matches(target_seq, expanded_words, word_size,
                                    positions)
        print(f"  Found {len(matches)} initial word matches")

        # Step 6: Extend each hit to form HSPs
        for q_pos, t_pos in matches:
            hsp = extend_hit(filtered_query, target_seq, q_pos, t_pos,
                             blosum62, word_size, threshold_score)

            # Keep only significant HSPs
            if hsp['score'] > threshold_score:
                # Calculate E-value for statistical significance
                hsp['e_value'] = calculate_e_value(hsp['score'],
                                                   len(filtered_query),
                                                   len(target_seq))
                results[target_header].append(hsp)

        print(f"  Found {len(results[target_header])} significant HSPs")

        # Keep only the best HSP per target (lowest E-value)
        if results[target_header]:
            results[target_header] = sorted(results[target_header],
                                            key=lambda x: x['e_value'])[0]
            print(f"  Best alignment E-value: "
                  f"{results[target_header]['e_value']:.2e}")
        else:
            # Create empty result for targets with no significant hits
            results[target_header] = {'score': 0, 'e_value': float('inf')}

    # Step 7: Rank all targets by E-value (ascending = better)
    ranked_targets = sorted(results.items(), key=lambda x: x[1]['e_value'])

    return ranked_targets


def create_alignment_output(query, ranked_results):
    """
    Format mini-BLAST results for output display.

    Args:
        query (str): Original query sequence
        ranked_results (list): List of ranked alignment results

    Returns:
        str: Formatted alignment output string
    """
    output = ["Mini-BLAST Results\n"]
    output.append(f"Query sequence length: {len(query)}\n")
    output.append("Ranked results:\n")

    for rank, (target_header, result) in enumerate(ranked_results, 1):
        output.append(f"Rank {rank}: {target_header}")

        if result['score'] > 0:
            output.append(f"  Score: {result['score']}")
            output.append(f"  E-value: {result['e_value']:.2e}")

            # Generate alignment display if alignment data exists
            if ('query_start' in result and 'q_align' in result and
                    't_align' in result):
                q_align = result['q_align']
                t_align = result['t_align']
                blosum = create_blosum62()

                # Generate match line showing identical/similar positions
                match_line = ""
                for i in range(len(q_align)):
                    if i < len(t_align):
                        q_aa, t_aa = q_align[i], t_align[i]
                        score = blosum.get((q_aa, t_aa), -4)

                        if q_aa == t_aa:
                            match_line += "|"  # Identical residues
                        elif score > 0:
                            match_line += "+"  # Similar residues
                        else:
                            match_line += " "  # Different residues

                # Format alignment with proper spacing
                q_start_str = str(result['query_start'])
                t_start_str = str(result['target_start'])

                # Calculate padding for visual alignment
                q_prefix_len = len("  Query: ") + len(q_start_str) + 1
                t_prefix_len = len("  Target: ") + len(t_start_str) + 1
                max_prefix_len = max(q_prefix_len, t_prefix_len)

                q_padding = " " * (max_prefix_len - q_prefix_len)
                t_padding = " " * (max_prefix_len - t_prefix_len)
                m_padding = " " * (max_prefix_len - len("  Match: "))

                # Output aligned sequences
                output.append(f"  Query: {q_start_str} {q_padding}{q_align} "
                              f"{result['query_end']}")
                output.append(f"  Match: {m_padding}{match_line}")
                output.append(f"  Target: {t_start_str} {t_padding}{t_align} "
                              f"{result['target_end']}")

                # Calculate and display identity/similarity statistics
                identity = sum(1 for i in range(len(match_line))
                               if match_line[i] == '|')
                similarity = identity + sum(1 for i in range(len(match_line))
                                            if match_line[i] == '+')
                length = len(match_line)

                identity_pct = (identity / length) * 100
                similarity_pct = (similarity / length) * 100

                output.append(f"  Identity: {identity}/{length} "
                              f"({identity_pct:.1f}%)")
                output.append(f"  Similarity: {similarity}/{length} "
                              f"({similarity_pct:.1f}%)")

                # Store metrics for potential later use
                result['identity_pct'] = identity_pct
                result['similarity_pct'] = similarity_pct
        else:
            output.append("  No significant alignments found")

        output.append("")

    return "\n".join(output)


def main():
    """
    Main function to execute mini-BLAST protein sequence alignment.
    """
    args = parse_arguments()

    # Read input sequences
    query_sequences = parse_fasta(args.query_file)
    target_sequences = parse_fasta(args.target_file)

    # Get query sequence for processing
    query_header = list(query_sequences.keys())[0]
    query_seq = query_sequences[query_header]
    print(f"Processing query: {query_header}")
    print(f"Found {len(target_sequences)} sequences in target database")

    # Execute mini-BLAST algorithm
    ranked_results = run_mini_blast(
        query_seq,
        target_sequences,
        3,  # word_size
        args.threshold_score,
        args.neighborhood_threshold
    )

    # Generate formatted output
    output = create_alignment_output(query_seq, ranked_results)

    # Save results to file
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(args.output_file, 'w') as f:
        f.write(output)
    print(f"Results saved to {args.output_file}")


if __name__ == '__main__':
    main()