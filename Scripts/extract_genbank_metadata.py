#!/usr/bin/env python3

import os
import argparse
from pathlib import Path
from Bio import SeqIO


def extract_feature_qualifier(record, feature_type, qualifier):
    for feature in record.features:
        if feature.type == feature_type and qualifier in feature.qualifiers:
            return feature.qualifiers[qualifier][0]
    return ""


def extract_metadata(record):
    organism = record.annotations.get("organism", "")
    molecule_type = record.annotations.get("molecule_type", "")
    strain = extract_feature_qualifier(record, "source", "strain")
    host = extract_feature_qualifier(record, "source", "host")

    taxon_raw = extract_feature_qualifier(record, "source", "db_xref")
    taxon_id = taxon_raw.split(":")[-1] if "taxon" in taxon_raw else ""

    segment = extract_feature_qualifier(record, "source", "segment")
    country = extract_feature_qualifier(record, "source", "country")
    collection_date = extract_feature_qualifier(record, "source", "collection_date")
    note = extract_feature_qualifier(record, "CDS", "note")
    codon_start = extract_feature_qualifier(record, "CDS", "codon_start")
    product = extract_feature_qualifier(record, "CDS", "product")
    protein_id = extract_feature_qualifier(record, "CDS", "protein_id")

    return [
        organism, molecule_type, strain, host, taxon_id, segment, country,
        collection_date, note, codon_start, product, protein_id
    ]


def main(input_dir, output_dir):
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / "parsed_metadata.tsv"

    with open(output_file, 'w') as out_file:
        out_file.write("File\tOrganism\tMolecule_Type\tStrain\tHost\tTaxon_ID\tSegment\tCountry\tCollection_Date\tNote\tCodon_Start\tProduct\tProtein_ID\n")

        for gb_file in sorted(input_path.glob("*.gb")):
            try:
                record = SeqIO.read(gb_file, "genbank")
                metadata = extract_metadata(record)
                out_file.write(f"{gb_file.name}\t" + "\t".join(metadata) + "\n")
            except Exception as e:
                print(f"[ERROR] Failed to parse {gb_file.name}: {e}")

    print(f"[DONE] Metadata written to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
Extract metadata fields from multiple GenBank (.gb) files and save them in a tab-delimited summary table.
This tool is useful for building curated metadata tables from viral or bacterial genome annotations.
""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:

  1. Extract metadata from all GenBank files in a folder:

     python extract_genbank_metadata.py -i ./gb_files -o ./output

  2. Your output file will be: ./output/parsed_metadata.tsv

Expected input:
  - GenBank flat files with `.gb` extension
  - Located in the specified input directory

Extracted fields:
  - Organism, Molecule type, Strain, Host, Taxon ID, Segment, Country, Collection date,
    Note, Codon start, Product, Protein ID
"""
    )

    parser.add_argument(
        "-i", "--input", required=True, metavar="INPUT_DIR",
        help="Directory containing .gb GenBank files"
    )
    parser.add_argument(
        "-o", "--output", required=True, metavar="OUTPUT_DIR",
        help="Directory to save the extracted metadata file"
    )

    args = parser.parse_args()
    main(args.input, args.output)
