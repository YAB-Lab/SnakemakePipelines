#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare concatenated CDS FASTA with modified headers for runPopGen and/or SAPA.")
    parser.add_argument("-i", "--input", required=True, help="Tab-delimited file with species ID and CDS fasta path per line.")
    parser.add_argument("-o", "--output", required=True, help="Output concatenated FASTA file.")
    return parser.parse_args()

def extract_gene_id(header):
    """Extract gene=... from the original FASTA header."""
    parts = header.split(';')
    for part in parts:
        if part.strip().startswith("gene="):
            return part.strip().split("=")[-1]
    return "UNKNOWN"

def modify_header(original_header, species_id, line_id):
    """Simplify header and append species and line ID."""
    header_main = original_header.split()[0]
    gene_id = extract_gene_id(original_header)
    return f"{header_main} gene={gene_id} species={species_id} line={line_id}"

def main():
    args = parse_args()
    output_records = []

    with open(args.input, 'r') as sample_file:
        for line in sample_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            species_id, fasta_path = line.split("\t")
            line_id = os.path.basename(fasta_path).replace(".CDS.fasta", "")

            for record in SeqIO.parse(fasta_path, "fasta"):
                record.id = modify_header(record.description, species_id, line_id)
                record.description = "" 
                output_records.append(record)

    with open(args.output, "w") as out_fasta:
        SeqIO.write(output_records, out_fasta, "fasta")

if __name__ == "__main__":
    main()
