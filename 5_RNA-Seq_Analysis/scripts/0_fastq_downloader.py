#!/usr/bin/env python3

import argparse
import gzip
import os
import shutil
import urllib.request


def parse_arguments():
    """
    Parse command-line arguments for FASTQ file downloader.
    """
    parser = argparse.ArgumentParser(
        description='Download ENCODE FASTQ file for RNA-seq analysis'
    )
    parser.add_argument(
        '--output_dir',
        default=('data/0_fastq_input'),
        help='Output directory where FASTQ file will be saved'
    )
    parser.add_argument(
        '--url',
        default=('https://www.encodeproject.org/files/ENCFF493KQW/'
                 '@@download/ENCFF493KQW.fastq.gz'),
        help='URL to download FASTQ file from'
    )
    return parser.parse_args()


def download_fastq_file(url, output_dir):
    """
    Download and extract ENCODE FASTQ file if not exists.

    Args:
        url (str): URL to download the FASTQ file from
        output_dir (str): Directory where FASTQ file should be saved

    Returns:
        str: Path to the downloaded and extracted FASTQ file
    """
    os.makedirs(output_dir, exist_ok=True)

    gz_path = os.path.join(output_dir, "ENCFF493KQW.fastq.gz")
    fastq_path = os.path.join(output_dir, "ENCFF493KQW.fastq")

    # Check if the file is already downloaded
    if os.path.exists(fastq_path):
        print(f"    ENCFF493KQW.fastq already exists at: {fastq_path}")
        print(f"    Skipping download...")
        return fastq_path

    print(f"    Downloading FASTQ file...")
    urllib.request.urlretrieve(url, gz_path)

    print(f"    Extracting...")
    with gzip.open(gz_path, "rb") as gz, open(fastq_path, "wb") as fastq:
        shutil.copyfileobj(gz, fastq)

    if os.path.exists(gz_path):
        os.remove(gz_path)

    print(f"    Downloaded and extracted to: {fastq_path}")
    return fastq_path


def main():
    """
    Main function to execute FASTQ file download.
    """
    args = parse_arguments()

    # Download FASTQ file
    download_fastq_file(args.url, args.output_dir)

    print("FASTQ download completed successfully!")


if __name__ == "__main__":
    main()