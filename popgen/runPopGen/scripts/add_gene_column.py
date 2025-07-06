#!/usr/bin/env python3
 
import sys
from Bio import SeqIO
import csv

fasta_file = snakemake.input.fasta
mkout_file = snakemake.input.mkout
output_file = snakemake.output[0]

# Step 1: Build transcript → gene mapping, ensuring uniqueness
transcript_to_gene = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    header_parts = record.description.strip().split()
    transcript_id = header_parts[0]  # e.g. rna-XM_070211082.1
    kv_pairs = dict(part.split("=", 1) for part in header_parts[1:] if "=" in part)
    gene = kv_pairs.get("gene", "NA")

    # Only map the transcript once — keep first mapping
    if transcript_id not in transcript_to_gene:
        transcript_to_gene[transcript_id] = gene


# Step 2: Read MKout and write with new "GENE" column
with open(mkout_file, newline='') as infile, open(output_file, "w", newline='') as outfile:
    reader = csv.DictReader(infile, delimiter="\t")
    fieldnames = ["GENE"] + reader.fieldnames
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        transcript = row["TRANSCRIPT"]
        gene_id = transcript_to_gene.get(transcript, "NA")
        row_with_gene = {"GENE": gene_id}
        row_with_gene.update(row)
        writer.writerow(row_with_gene)
