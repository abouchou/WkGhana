#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import Entrez, SeqIO


def download_genbank_records(accession_numbers, output_folder):
    Entrez.email = "christianpaulako@gmail.com"  # Update if needed

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)  # Create output folder if not exists

    for accession in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            filename = output_path / f"{accession}.gb"
            with open(filename, "w") as output_file:
                SeqIO.write(record, output_file, "genbank")
            print(f"[INFO] Downloaded {accession} → {filename}")
        except Exception as e:
            print(f"[ERROR] Could not download {accession}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="""
Download GenBank records using accession numbers from a file.
Each accession should be on a separate line.
The script will fetch .gb files and save them in the specified output directory.
""",
        epilog="""
Example:
  python download_genbank_records.py -acc accessions.txt -out ./genbank_files

Required:
  - A text file with accession numbers (one per line)
  - Internet connection (uses NCBI Entrez API)
"""
    )

    parser.add_argument("--accession", "-acc", required=True,
                        help="File containing GenBank accession numbers (one per line)")
    parser.add_argument("--output", "-out", required=True,
                        help="Output directory to save downloaded .gb files")

    args = parser.parse_args()

    try:
        with open(args.accession) as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"[ERROR] Accession file not found: {args.accession}")
        return

    try:
        download_genbank_records(accession_numbers, args.output)
        print(f"[DONE] All downloads saved in: {args.output}")
    except Exception as e:
        print(f"[FATAL] Script failed: {e}")


if __name__ == "__main__":
    main()
