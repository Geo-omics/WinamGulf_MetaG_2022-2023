from Bio import Entrez, SeqIO
import csv


def fetch_genbank_sequence(genbank_id, start, end):
    try:
        Entrez.email = "lnhart@umich.edu"  # Set your email here
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        sequence = record.seq[start-1:end]  # Adjust for 0-based indexing
        return str(sequence)
    except Exception as e:
        print(f"Error fetching sequence for Genbank ID {genbank_id}: {e}")
        return None

def create_fasta(csv_file, output_fasta):
    with open(csv_file, 'r') as csvfile, open(output_fasta, 'w') as fastafile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            genbank_id = row['\ufeffgenbank_id']
            start = int(row['start'])
            end = int(row['end'])
            db_name = row['db_name']

            sequence = fetch_genbank_sequence(genbank_id, start, end)
            if sequence is not None:
                fastafile.write(f">{db_name}\n{sequence}\n")

if __name__ == "__main__":
    csv_file_path = "db_input_genbank.csv"  # Replace with your CSV file path
    output_fasta_path = "output.fasta"  # Replace with your desired output FASTA file path
    create_fasta(csv_file_path, output_fasta_path)


