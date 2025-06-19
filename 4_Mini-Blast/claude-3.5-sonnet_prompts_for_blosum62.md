# Note
This prompt discussion was created during the CIS5002 assignment in creating a blosum62 alignment function.
Within this task, as confirmed by Prof. Jean Paul Ebejer, I am able to reuse parts of my previous code to
aid in the MMB5009 assignment.
The script I am recycling was created with the aid of Claude 3.5 using the following prompt discussion.


# Prompt 1
when matching the ORF protein sequences, is this the best method to find the best possible match to the reference?

``` python
# finds the best position to align orf_sequence sequence with reference
def find_best_alignment(reference, orf_sequence):
    best_score = 0
    best_position = 0

    # trying each possible position on the reference sequence using the orf sequence
    for i in range(len(reference) - len(orf_sequence) + 1):
        # getting the slice of reference sequence to compare
        ref_slice = reference[i:i + len(orf_sequence)]

        # zipping the sequences together for comparison
        paired_sequences = zip(ref_slice, orf_sequence)

        # counting the matching characters relative to position
        matches = sum(1 for x, y in paired_sequences if x == y)
        if matches > best_score:
            best_score = matches
            best_position = i

    match_percentage = (best_score / len(orf_sequence)) * 100
    return best_position, match_percentage
```

## Caluda reply
Let me analyze the current alignment approach and suggest improvements.
<answer_operator>
Looking at sequence alignment methods for protein sequences, particularly for ORFs against reference sequences.
The current method uses a basic sliding window approach with exact character matching, which has several limitations for protein sequence alignment:

No consideration of amino acid similarity/conservation
No gap penalties
No consideration of protein domains or motifs
Only finds local alignments without consideration of global alignment potential

Here's a suggested improved version:
``` python
from Bio.SubsMat import MatrixInfo as matlist

def find_best_alignment(reference, orf_sequence, gap_open=-10, gap_extend=-0.5):
    """
    Performs local alignment using BLOSUM62 scoring matrix

    Args:
        reference: Reference protein sequence
        orf_sequence: ORF translated sequence
        gap_open: Gap opening penalty (default: -10)
        gap_extend: Gap extension penalty (default: -0.5)

    Returns:
        best_position: Best starting position in reference
        match_percentage: Percentage of matching/similar residues
        alignment: Tuple of aligned sequences with gaps
    """
    # Use BLOSUM62 substitution matrix
    blosum62 = matlist.blosum62

    # Initialize scoring matrix
    rows = len(orf_sequence) + 1
    cols = len(reference) + 1
    score_matrix = [[0] * cols for _ in range(rows)]
    traceback = [[None] * cols for _ in range(rows)]

    # Track maximum score
    max_score = 0
    max_pos = (0, 0)

    # Fill scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # Get amino acid pair score from BLOSUM62
            pair = (orf_sequence[i-1], reference[j-1])
            if pair in blosum62:
                match_score = blosum62[pair]
            else:
                pair = (pair[1], pair[0])
                match_score = blosum62.get(pair, -4)  # Default to -4 for unknown pairs

            # Calculate scores for different moves
            diagonal = score_matrix[i-1][j-1] + match_score
            up = score_matrix[i-1][j] + (gap_extend if i > 1 and traceback[i-1][j] == 'up' else gap_open)
            left = score_matrix[i][j-1] + (gap_extend if j > 1 and traceback[i][j-1] == 'left' else gap_open)

            # Find best score
            score_matrix[i][j] = max(0, diagonal, up, left)

            # Track maximum
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

            # Record traceback
            if score_matrix[i][j] == 0:
                traceback[i][j] = None
            elif score_matrix[i][j] == diagonal:
                traceback[i][j] = 'diagonal'
            elif score_matrix[i][j] == up:
                traceback[i][j] = 'up'
            else:
                traceback[i][j] = 'left'

    # Traceback to get alignment
    align1 = []
    align2 = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback[i][j] == 'diagonal':
            align1.append(orf_sequence[i-1])
            align2.append(reference[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'up':
            align1.append(orf_sequence[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(reference[j-1])
            j -= 1

    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    # Calculate similarity percentage including similar amino acids
    similar_count = 0
    for a, b in zip(align1, align2):
        if a != '-' and b != '-':
            pair = (a, b)
            if pair in blosum62 and blosum62[pair] > 0:
                similar_count += 1

    match_percentage = (similar_count / len(align1)) * 100

    return j, match_percentage, (align1, align2)
```

Key improvements:

Uses BLOSUM62 substitution matrix to account for amino acid similarities
Implements gap opening and extension penalties
Returns aligned sequences with gaps
Considers both exact matches and similar amino acids in scoring
Uses dynamic programming for optimal local alignment

You'll need to install Biopython (pip install biopython) for the BLOSUM62 matrix.
Y
</answer_operator>

## Edward's Reply
do not use bio or any major packages, please explain it each step extensively. use a simulated example on what is happening (i.e. make a reference sequence and ORF sequence)

## Clauda's Reply
Let me break down protein sequence alignment from scratch with a clear example.

<answer_operator>
Let's rebuild the alignment function using basic Python. I'll demonstrate with sample sequences to show each step.

```python
def create_blosum62():
    """Creates a simplified BLOSUM62 scoring matrix as a dictionary"""
    # This is a reduced BLOSUM62 matrix focusing on key amino acid relationships
    blosum62 = {
        # Identical matches
        ('A', 'A'): 4,  ('R', 'R'): 5,  ('N', 'N'): 6,  ('D', 'D'): 6,
        ('C', 'C'): 9,  ('Q', 'Q'): 5,  ('E', 'E'): 5,  ('G', 'G'): 6,
        ('H', 'H'): 8,  ('I', 'I'): 4,  ('L', 'L'): 4,  ('K', 'K'): 5,
        ('M', 'M'): 5,  ('F', 'F'): 6,  ('P', 'P'): 7,  ('S', 'S'): 4,
        ('T', 'T'): 5,  ('W', 'W'): 11, ('Y', 'Y'): 7,  ('V', 'V'): 4,

        # Common positive scoring pairs (similar amino acids)
        ('S', 'T'): 1,  ('T', 'S'): 1,  # Serine/Threonine
        ('D', 'E'): 2,  ('E', 'D'): 2,  # Aspartic/Glutamic acid
        ('K', 'R'): 2,  ('R', 'K'): 2,  # Lysine/Arginine
        ('I', 'L'): 2,  ('L', 'I'): 2,  # Isoleucine/Leucine
        ('F', 'Y'): 3,  ('Y', 'F'): 3   # Phenylalanine/Tyrosine
    }
    return blosum62

def find_best_protein_alignment(reference, orf_sequence, gap_open=-10, gap_extend=-0.5):
    """
    Performs local alignment of protein sequences using a scoring matrix

    Args:
        reference: Reference protein sequence
        orf_sequence: ORF translated sequence
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
    """
    # Get scoring matrix
    blosum = create_blosum62()

    # Initialize scoring matrix
    rows = len(orf_sequence) + 1
    cols = len(reference) + 1
    score_matrix = [[0] * cols for _ in range(rows)]
    traceback = [[None] * cols for _ in range(rows)]

    # Track best score
    max_score = 0
    max_pos = (0, 0)

    # Fill scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            # Get amino acid pair
            aa1 = orf_sequence[i-1]
            aa2 = reference[j-1]

            # Get match score
            pair = (aa1, aa2)
            match_score = blosum.get(pair, -4)  # Default to -4 for dissimilar amino acids

            # Calculate possible moves
            diagonal = score_matrix[i-1][j-1] + match_score
            up = score_matrix[i-1][j] + (gap_extend if i > 1 and traceback[i-1][j] == 'up' else gap_open)
            left = score_matrix[i][j-1] + (gap_extend if j > 1 and traceback[i][j-1] == 'left' else gap_open)

            # Find best score
            score_matrix[i][j] = max(0, diagonal, up, left)

            # Update maximum if needed
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

            # Record move for traceback
            if score_matrix[i][j] == 0:
                traceback[i][j] = None
            elif score_matrix[i][j] == diagonal:
                traceback[i][j] = 'diagonal'
            elif score_matrix[i][j] == up:
                traceback[i][j] = 'up'
            else:
                traceback[i][j] = 'left'

    # Traceback to get alignment
    align1, align2 = [], []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback[i][j] == 'diagonal':
            align1.append(orf_sequence[i-1])
            align2.append(reference[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'up':
            align1.append(orf_sequence[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(reference[j-1])
            j -= 1

    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    # Calculate similarity statistics
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    similars = sum(1 for a, b in zip(align1, align2)
                  if (a, b) in blosum and blosum[(a, b)] > 0 and a != '-' and b != '-')

    alignment_length = len(align1)
    match_percentage = (matches / alignment_length) * 100
    similarity_percentage = (similars / alignment_length) * 100

    return j, match_percentage, similarity_percentage, (align1, align2)

# Let's demonstrate with an example:
def demonstrate_alignment():
    # Sample sequences
    reference = "MTVRLGATASGLWLFHQGAWRRLPGPGRLRLGALRSRFLTRVRNVQPRAGPGGPRGPQR"
    orf_seq = "MTVRLGATASGLWL"  # First part of reference (perfect match)

    print("Reference:", reference)
    print("ORF Sequence:", orf_seq)
    print("\nRunning alignment...\n")

    start_pos, match_pct, sim_pct, (align1, align2) = find_best_protein_alignment(reference, orf_seq)

    print("Alignment Result:")
    print("Start position in reference:", start_pos)
    print("Match percentage:", f"{match_pct:.1f}%")
    print("Similarity percentage:", f"{sim_pct:.1f}%")
    print("\nAlignment visualization:")
    print("ORF:      ", align1)
    print("Reference:", align2)
    print("\nMatch indicators:")
    print("           " + ''.join('|' if a == b else ':' if (a, b) in create_blosum62() else ' '
                                 for a, b in zip(align1, align2)))

# Run demonstration
demonstrate_alignment()
```

Let's break down what's happening:

1. **Scoring Matrix Creation**:
   - We create a simplified BLOSUM62 matrix that gives scores for amino acid pairs
   - Identical matches get high positive scores
   - Similar amino acids (like D-E or K-R) get lower positive scores
   - Dissimilar pairs get negative scores (default -4)

2. **Dynamic Programming Matrix**:
   - Create a matrix with ORF sequence on rows and reference on columns
   - Each cell represents the best score possible ending at that position
   - Three possible moves at each position:
     * Diagonal (match/mismatch)
     * Up (gap in reference)
     * Left (gap in ORF)

3. **Gap Penalties**:
   - Opening a gap is expensive (default -10)
   - Extending a gap is cheaper (default -0.5)
   - This encourages fewer, longer gaps rather than many short gaps

4. **Traceback**:
   - Start from highest scoring position
   - Follow the path back to reconstruct alignment
   - Record moves to build aligned sequences with gaps

Example Output:
```
Reference: MTVRLGATASGLWLFHQGAWRRLPGPGRLRLGALRSRFLTRVRNVQPRAGPGGPRGPQR
ORF Sequence: MTVRLGATASGLWL

Alignment Result:
Start position in reference: 0
Match percentage: 100.0%
Similarity percentage: 100.0%

Alignment visualization:
ORF:       MTVRLGATASGLWL
Reference: MTVRLGATASGLWL
           ||||||||||||||
```

The example shows a perfect match at the start of the reference. In real cases, you might see:
- Partial matches with gaps (shown as '-')
- Similar amino acids (shown as ':')
- Non-matching positions (shown as spaces)

This approach is more biologically relevant than the original because it:
1. Considers amino acid properties (through BLOSUM62)
2. Handles insertions/deletions properly
3. Differentiates between similar and dissimilar amino acids
4. Provides detailed alignment information

Y

</answer_operator>


## Edward's Reply
explain more on def create_blosum62()


## Clauda's Reply
Let me break down the BLOSUM (BLOcks SUbstitution Matrix) 62 scoring system and how it works in sequence alignments.

<answer_operator>
The BLOSUM62 matrix represents the likelihood of amino acid substitutions based on observed substitutions in conserved protein blocks that share 62% sequence identity. Let's expand the function with more comprehensive scoring:

```python
def create_blosum62():
    """
    Creates BLOSUM62 scoring matrix for protein alignment

    Returns:
        Dictionary with amino acid pair tuples as keys and substitution scores as values

    Scoring Logic:
    - Higher positive scores = More likely/favorable substitutions
    - Lower/negative scores = Less likely/unfavorable substitutions
    - Based on evolutionary conservation patterns
    """
    blosum62 = {
        # Identical matches (strongly conserved)
        ('A', 'A'): 4,   # Alanine
        ('R', 'R'): 5,   # Arginine
        ('N', 'N'): 6,   # Asparagine
        ('D', 'D'): 6,   # Aspartic acid
        ('C', 'C'): 9,   # Cysteine (highest self-score due to disulfide bonds)
        ('Q', 'Q'): 5,   # Glutamine
        ('E', 'E'): 5,   # Glutamic acid
        ('G', 'G'): 6,   # Glycine
        ('H', 'H'): 8,   # Histidine
        ('I', 'I'): 4,   # Isoleucine
        ('L', 'L'): 4,   # Leucine
        ('K', 'K'): 5,   # Lysine
        ('M', 'M'): 5,   # Methionine
        ('F', 'F'): 6,   # Phenylalanine
        ('P', 'P'): 7,   # Proline
        ('S', 'S'): 4,   # Serine
        ('T', 'T'): 5,   # Threonine
        ('W', 'W'): 11,  # Tryptophan (highest due to rarity and size)
        ('Y', 'Y'): 7,   # Tyrosine
        ('V', 'V'): 4,   # Valine

        # Similar physicochemical properties (positive scores)
        # Hydrophobic amino acids
        ('I', 'L'): 2, ('L', 'I'): 2,  # Isoleucine/Leucine
        ('I', 'V'): 3, ('V', 'I'): 3,  # Isoleucine/Valine
        ('L', 'V'): 1, ('V', 'L'): 1,  # Leucine/Valine
        ('L', 'M'): 2, ('M', 'L'): 2,  # Leucine/Methionine

        # Aromatic amino acids
        ('F', 'Y'): 3, ('Y', 'F'): 3,  # Phenylalanine/Tyrosine
        ('W', 'F'): 1, ('F', 'W'): 1,  # Tryptophan/Phenylalanine
        ('W', 'Y'): 2, ('Y', 'W'): 2,  # Tryptophan/Tyrosine

        # Polar/charged amino acids
        ('K', 'R'): 2, ('R', 'K'): 2,  # Lysine/Arginine (positive)
        ('D', 'E'): 2, ('E', 'D'): 2,  # Aspartic acid/Glutamic acid (negative)
        ('N', 'D'): 1, ('D', 'N'): 1,  # Asparagine/Aspartic acid
        ('Q', 'E'): 2, ('E', 'Q'): 2,  # Glutamine/Glutamic acid
        ('S', 'T'): 1, ('T', 'S'): 1,  # Serine/Threonine

        # Small amino acids
        ('A', 'G'): -1, ('G', 'A'): -1,  # Alanine/Glycine
        ('A', 'S'): 1,  ('S', 'A'): 1,   # Alanine/Serine

        # Dissimilar amino acids (negative scores)
        ('W', 'G'): -3, ('G', 'W'): -3,  # Large aromatic vs smallest
        ('K', 'C'): -3, ('C', 'K'): -3,  # Charged vs special
        ('R', 'P'): -2, ('P', 'R'): -2,  # Charged vs structure breaker
    }
    return blosum62

# Example usage demonstration
def demonstrate_blosum62():
    matrix = create_blosum62()

    # Example protein sequence segment
    seq1 = "LVFYW"  # Hydrophobic and aromatic
    seq2 = "IVFFW"  # Similar but not identical

    print("Comparing sequences:")
    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nScoring each position:")

    total_score = 0
    for i, (aa1, aa2) in enumerate(zip(seq1, seq2)):
        pair = (aa1, aa2)
        score = matrix.get(pair, -4)  # Default -4 for unlisted pairs
        total_score += score
        print(f"Position {i+1}: {aa1}-{aa2} = {score}")
        print(f"Reason: ", end="")

        if aa1 == aa2:
            print(f"Identical match ({aa1})")
        elif score > 0:
            print(f"Similar amino acids (positive score)")
        else:
            print(f"Dissimilar amino acids (negative score)")

    print(f"\nTotal alignment score: {total_score}")

# Run demonstration
demonstrate_blosum62()
```

This will output something like:
```
Comparing sequences:
Sequence 1: LVFYW
Sequence 2: IVFFW

Scoring each position:
Position 1: L-I = 2
Reason: Similar amino acids (positive score)
Position 2: V-V = 4
Reason: Identical match (V)
Position 3: F-F = 6
Reason: Identical match (F)
Position 4: Y-F = 3
Reason: Similar amino acids (positive score)
Position 5: W-W = 11
Reason: Identical match (W)

Total alignment score: 26
```

Key Points about BLOSUM62:

1. **Identical Matches**:
   - Highest scores given to identical amino acids
   - Scores vary based on amino acid conservation rates
   - Rare amino acids (like W) have higher self-scores

2. **Similar Amino Acids**:
   - Positive scores for common substitutions
   - Based on physicochemical properties:
     * Hydrophobic (I, L, V, M)
     * Aromatic (F, Y, W)
     * Charged (K, R, D, E)
     * Small (A, G, S)

3. **Dissimilar Amino Acids**:
   - Negative scores for unlikely substitutions
   - Larger penalties for very different properties
   - Default -4 for unlisted pairs

4. **Scoring Logic**:
   - Higher scores = More likely substitutions
   - Based on observed frequencies in related proteins
   - Helps identify biologically meaningful alignments

5. **Applications**:
   - Protein sequence alignment
   - Homology detection
   - Evolutionary relationship assessment
   - Protein structure prediction

Y

</answer_operator>

## Edward's Reply
i found the following paper that assigns the amino acid scores to one another, review the uploaded henikoff-henikoff-1992-amino-acid-substitution-matrices-from-protein-blocks.pdf in your project knowledge

i want you to rewrite def create_blosum62() so as to mirror this papers scoring matrix

## Claude's Reply
I'll analyze the Henikoff & Henikoff 1992 paper and implement their exact BLOSUM62 scoring matrix.

<answer_operator>
Looking at the paper, the BLOSUM62 matrix shown in Figure 2 (page 10917) provides the actual scores used for amino acid substitutions. Let me rewrite the function to exactly match these values:

```python
def create_blosum62():
    """
    Creates BLOSUM62 substitution matrix based on Henikoff & Henikoff 1992 paper.
    All scores are rounded to nearest integer and scaled by factor of 2.

    Returns:
        Dictionary with amino acid pair tuples as keys and substitution scores as values
    """
    # Initialize empty matrix
    blosum62 = {}

    # Define amino acids in order they appear in the paper's matrix
    amino_acids = ['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q',
                   'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

    # Scores from Figure 2 (lower triangle + diagonal)
    scores = [
        # C   S   T   P   A   G   N   D   E   Q   H   R   K   M   I   L   V   F   Y   W
        [ 9], # C
        [-1,  4], # S
        [-1,  1,  5], # T
        [-3, -1, -1,  7], # P
        [ 0,  1,  0, -1,  4], # A
        [-3,  0, -2, -2,  0,  6], # G
        [-3,  1,  0, -2, -2,  0,  6], # N
        [-3,  0, -1, -1, -2, -1,  1,  6], # D
        [-4,  0, -1, -1, -1, -2,  0,  2,  5], # E
        [-3,  0, -1, -1, -1, -2,  0,  0,  2,  5], # Q
        [-3, -1, -2, -2, -2, -2,  1, -1,  0,  0,  8], # H
        [-3, -1, -1, -2, -1, -2,  0, -2,  0,  1,  0,  5], # R
        [-3,  0, -1, -1, -1, -2,  0, -1,  1,  1, -1,  2,  5], # K
        [-1, -1, -1, -2, -1, -3, -2, -3, -2,  0, -2, -1, -1,  5], # M
        [-1, -2, -1, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3,  1,  4], # I
        [-1, -2, -1, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2,  2,  2,  4], # L
        [-1, -2,  0, -2,  0, -3, -3, -3, -2, -2, -3, -3, -2,  1,  3,  1,  4], # V
        [-2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3,  0,  0,  0, -1,  6], # F
        [-2, -2, -2, -3, -2, -3, -2, -3, -2, -1,  2, -2, -2, -1, -1, -1, -1,  3,  7], # Y
        [-2, -3, -2, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3,  1,  2, 11]  # W
    ]

    # Fill the dictionary with scores
    for i, aa1 in enumerate(amino_acids):
        for j, aa2 in enumerate(amino_acids[:i+1]):
            # Get score from lower triangle
            score = scores[i][j]
            # Add both directions to dictionary
            blosum62[(aa1, aa2)] = score
            blosum62[(aa2, aa1)] = score

    return blosum62

# Demonstration function
def demonstrate_blosum62_paper():
    matrix = create_blosum62()

    # Print header
    amino_acids = ['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q',
                   'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

    print("BLOSUM62 Matrix from Henikoff & Henikoff 1992:")
    print("    " + "  ".join(f"{aa:2}" for aa in amino_acids))

    # Print matrix
    for i, aa1 in enumerate(amino_acids):
        row = [aa1]
        for aa2 in amino_acids:
            score = matrix.get((aa1, aa2), 0)
            row.append(f"{score:2}")
        print(" ".join(row))

    # Example comparison
    seq1 = "LVFYW"
    seq2 = "IVFFW"

    print("\nExample alignment scoring:")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")

    total_score = 0
    for i, (aa1, aa2) in enumerate(zip(seq1, seq2)):
        score = matrix[(aa1, aa2)]
        total_score += score
        print(f"Position {i+1}: {aa1}-{aa2} = {score:2}")

    print(f"Total score: {total_score}")

# Run demonstration
demonstrate_blosum62_paper()
```

Key differences from my previous version:

1. **Exact Scores**: All scores now match Figure 2 from the paper exactly

2. **Score Range**:
   - Highest: 11 (W-W)
   - Lowest: -4 (several pairs, e.g., G-I)
   - Most common positive: 4-6 for identical matches

3. **Notable Relationships**:
   - Cysteine (C) has high self-score (9) due to disulfide bonding
   - Tryptophan (W) has highest self-score (11) due to rarity
   - Conservative substitutions (e.g., I-L, D-E) score positively
   - Dissimilar amino acids score negatively

4. **Matrix Structure**:
   - Symmetric (score(A,B) = score(B,A))
   - Based on observed frequencies in protein blocks
   - Clustered at 62% identity threshold

This version is more accurate to the historical record and the actual BLOSUM62 matrix that became a standard in the field.

Y

</answer_operator>

---
## Note

Due to this, I used the method to make 5_BLOSUM62-alignment_per_amino-acid.py

---


