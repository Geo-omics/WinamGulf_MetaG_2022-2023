def modify_fasta_headers(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                parts = line.split()
                header_part = parts[0]
                new_header = header_part + "_ChlN\n"
                f_out.write(new_header)
            else:
                f_out.write(line)

input_file = 'ChlN_01052024.fna'
output_file = 'ChlN_01052024_LH.fna'
modify_fasta_headers(input_file, output_file)