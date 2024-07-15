import os
import pandas as pd
from Bio import SeqIO

def extract_contigs_from_assemblies(csv_file, assembly_dir):
    # Read CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file)

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        sample_id = row['sample_id']
        qseqid = row['qseqid']
        new_header = row['new_header']
        
        # Construct path to fasta file
        fasta_file = os.path.join(assembly_dir, f"{sample_id}.fa")
        if not os.path.exists(fasta_file):
            print(f"Fasta file {fasta_file} not found.")
            continue

        # Load corresponding fasta file
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Find contig with matching qseqid
            if record.id == qseqid:
                # Write contig with new header to output fasta file
                output_file = f"{new_header}.fasta"
                with open(output_file, "a") as f:
                    f.write(f">{new_header}\n")
                    f.write(f"{record.seq}\n")

# Usage
csv_file = "/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/mcy-cyr/cyr-blastcontigs/references/cyr_contigs_assembly_blast.csv"
assembly_dir = "/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/n_fixation/assemblies"
extract_contigs_from_assemblies(csv_file, assembly_dir)


