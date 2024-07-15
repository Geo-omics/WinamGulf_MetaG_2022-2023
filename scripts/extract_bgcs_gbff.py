import os
import csv
from Bio import SeqIO

def extract_gene_sequences(csv_file, gbff_folder, output_fasta):
    # Create a dictionary to store gene sequences
    gene_sequences = {}

    # Read the CSV file
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            file_name = row['file'] + '.gbff'
            record_id = row['record_id']
            start = int(row['start'])
            end = int(row['end'])
            new_header = row['new_header']

            # Read the GenBank file
            gbff_file = os.path.join(gbff_folder, file_name)
            for record in SeqIO.parse(gbff_file, 'genbank'):
                if record.id == record_id:
                    # Check if start is 0, if so, set it to 1
                    if start == 0:
                        start = 1
                    # Extract the gene sequence based on start and end coordinates
                    gene_sequence = record.seq[start-1:end]
                    # Store the gene sequence with the new header
                    gene_sequences[new_header] = gene_sequence

    # Write the extracted gene sequences to a single FASTA file
    with open(output_fasta, 'w') as fastafile:
        for header, sequence in gene_sequences.items():
            fastafile.write(f'>{header}\n{sequence}\n')

# Input CSV file
csv_file = '/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/references/parsing_drep95_man_annot.csv'
# Folder containing GenBank files
gbff_folder = '/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/data/quality_bins/antismash7/multismash/example/genbanks_qualitydrep99'
# Output FASTA file
output_fasta = '/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/references/drep95_manannot_bgcs_feb262024.fasta'

# Extract gene sequences and write to FASTA file
extract_gene_sequences(csv_file, gbff_folder, output_fasta)


