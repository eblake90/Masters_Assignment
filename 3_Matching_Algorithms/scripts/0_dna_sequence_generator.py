#!/usr/bin/env python3

import argparse
import gzip
import os
import random
import shutil
import urllib.request


def parse_arguments():
    """
    Parse command-line arguments for DNA sequence generator.
    """
    parser = argparse.ArgumentParser(
        description='Generate a modified copy of human chromosome 22 with '
                    'synthetic motifs.'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--output_dir',
        default='data/1_dna_sequence_data',
        help='Directory where output files will be saved'
    )
    parser.add_argument(
        '--motif_lengths',
        default='15,50',
        help='Comma-separated motif lengths to generate (default: 15,50)'
    )
    parser.add_argument(
        '--families_per_length',
        type=int,
        default=2,
        help='How many families (unique cores) per length (default: 2)'
    )
    parser.add_argument(
        '--max_k',
        type=int,
        default=5,
        help='Maximum Hamming-distance mutants to generate (default: 5)'
    )
    return parser.parse_args()


def download_chromosome(output_dir):
    """
    Download and return path to GRCh38 chromosome 22 FASTA if not exists.

    Args:
        output_dir (str): Directory where chromosome file should be saved

    Returns:
        str: Path to the downloaded and extracted FASTA file
    """
    url = ("https://ftp.ensembl.org/pub/release-114/fasta/"
           "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz")

    os.makedirs(output_dir, exist_ok=True)

    gz_path = os.path.join(output_dir, "GRCh38_chr22.fa.gz")
    fa_path = os.path.join(output_dir, "GRCh38_chr22.fa")

    # Check if the file is already downloaded
    if os.path.exists(fa_path):
        print(f"    GRCh38_chr22.fa already exists at: {fa_path}")
        print(f"    Skipping download...")
        return fa_path

    print(f"    Downloading chromosome 22...")
    urllib.request.urlretrieve(url, gz_path)

    print(f"    Extracting...")
    with gzip.open(gz_path, "rb") as gz, open(fa_path, "wb") as fa:
        shutil.copyfileobj(gz, fa)

    if os.path.exists(gz_path):
        os.remove(gz_path)

    print(f"    Downloaded and extracted to: {fa_path}")
    return fa_path


def read_fasta(file_path):
    """
    Read a FASTA file and return the header and sequence.

    Args:
        file_path (str): Path to the FASTA file to read

    Returns:
        tuple: (header string, sequence string)
    """
    print(f"\nReading FASTA file from: {file_path}")

    with open(file_path, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().lstrip('>')
    sequence = ''.join(line.strip() for line in lines[1:]).upper()

    print(f"  Loaded chromosome sequence: {len(sequence)} nucleotides")
    return header, sequence


def is_unique(word, reference):
    """
    Check if a word does not occur in the reference sequence.

    Args:
        word (str): DNA sequence to check for uniqueness
        reference (str): Reference sequence to search in

    Returns:
        bool: True if word is unique (not found in reference), False otherwise
    """
    return reference.find(word) == -1


def generate_random_word(length, rng):
    """
    Generate a random DNA word of specified length.

    Args:
        length (int): Length of the DNA sequence to generate
        rng (random.Random): Random number generator instance

    Returns:
        str: Random DNA sequence of specified length
    """
    return ''.join(rng.choice("ACGT") for _ in range(length))


def generate_hamming_variants(core, max_k, rng):
    """
    Generate variants of a core sequence with Hamming distances 1 to max_k.

    Args:
        core (str): Core DNA sequence to generate variants from
        max_k (int): Maximum Hamming distance for variants
        rng (random.Random): Random number generator instance

    Returns:
        list: List of DNA sequences including core and all variants
    """
    variants = [core]
    bases = "ACGT"

    for k in range(1, max_k + 1):
        word = list(core)
        positions = rng.sample(range(len(core)), k)

        for pos in positions:
            original = word[pos]
            options = [b for b in bases if b != original]
            word[pos] = rng.choice(options)

        variants.append(''.join(word))

    return variants


def generate_unique_motifs(reference, lengths, families_per_length, max_k,
                           seed):
    """
    Generate a set of motifs that don't naturally occur in the reference.

    Args:
        reference (str): Reference DNA sequence to avoid duplicates with
        lengths (list): List of motif lengths to generate
        families_per_length (int): Number of motif families per length
        max_k (int): Maximum Hamming distance for variants
        seed (int): Random seed for reproducibility

    Returns:
        list: List of tuples containing (motif_name, motif_sequence)
    """
    print("\nGenerating unique synthetic motifs:")

    rng = random.Random(seed)
    motifs = []
    family_id = 0

    for length in lengths:
        print(f"    Generating motifs of length {length}...")

        for _ in range(families_per_length):
            # Finding a unique core sequence
            while True:
                core = generate_random_word(length, rng)
                if is_unique(core, reference):
                    variants = generate_hamming_variants(core, max_k, rng)
                    if all(is_unique(v, reference) for v in variants[1:]):
                        break

            # Adding all variants to motifs
            for k, seq in enumerate(variants):
                name = f"M{family_id:02d}_len{length}_ham{k}"
                motifs.append((name, seq))

            family_id += 1

    total_motifs = len(motifs)
    print(f"  Generated {total_motifs} unique synthetic motifs")
    return motifs


def insert_motifs(sequence, motifs, seed):
    """
    Insert motifs at non-overlapping positions in the sequence.

    Args:
        sequence (str): DNA sequence to insert motifs into
        motifs (list): List of tuples containing (motif_name, motif_sequence)
        seed (int): Random seed for reproducible positioning

    Returns:
        tuple: (modified_sequence string, inserted_motifs list)
    """
    print("\nInserting motifs into sequence:")

    rng = random.Random(seed)
    seq_list = list(sequence)
    inserted_motifs = []

    for motif_name, motif_seq in motifs:
        motif_len = len(motif_seq)

        # Finding a non-overlapping position
        while True:
            position = rng.randint(0, len(sequence) - motif_len)

            # Checking for overlap with previously inserted motifs
            overlap = False
            for _, start, end, _ in inserted_motifs:
                if position < end and start < position + motif_len:
                    overlap = True
                    break

            if not overlap:
                break

        # Inserting the motif
        seq_list[position:position + motif_len] = motif_seq
        inserted_motifs.append((motif_name, position, position + motif_len,
                                motif_seq))

    print(f"    Inserted {len(inserted_motifs)} motifs at non-overlapping "
          f"positions")
    return ''.join(seq_list), inserted_motifs


def save_output_files(header, sequence, motif_records, output_dir):
    """
    Save the modified sequence and motif positions to files.
    """
    print("\nSaving output files:")

    # Creating output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Saving the modified FASTA file
    fasta_path = os.path.join(output_dir, "GRCh38_chr22_with_motifs.fa")
    with open(fasta_path, "w") as f:
        f.write(f">{header} | with synthetic motifs\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i + 60] + "\n")
    print(f"  Saved modified FASTA to: {fasta_path}")

    # Saving the motif positions CSV
    csv_path = os.path.join(output_dir, "chr22_motif_positions.csv")
    with open(csv_path, "w") as f:
        f.write("motif_name,start,end,sequence\n")
        for name, start, end, motif in motif_records:
            f.write(f"{name},{start},{end},{motif}\n")
    print(f"  Saved motif positions to: {csv_path}")


def main():
    """
    Main function to execute the DNA sequence generation pipeline.
    """
    # Parsing command-line arguments
    args = parse_arguments()

    # Converting comma-separated motif lengths to a list of integers
    motif_lengths = [int(x.strip()) for x in args.motif_lengths.split(",")
                     if x.strip()]

    # Downloading and read chromosome 22
    chr_path = download_chromosome(args.output_dir)
    header, sequence = read_fasta(chr_path)

    # Generating unique motifs
    motifs = generate_unique_motifs(
        sequence,
        motif_lengths,
        args.families_per_length,
        args.max_k,
        args.seed
    )

    # Inserting motifs into the sequence
    modified_sequence, motif_records = insert_motifs(sequence, motifs,
                                                     args.seed)

    # Saving the results
    save_output_files(header, modified_sequence, motif_records,
                      args.output_dir)

    print("\nDNA sequence generation completed successfully!")


if __name__ == "__main__":
    main()