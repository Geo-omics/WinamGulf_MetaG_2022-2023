def combine_fasta_files(input_files, output_file):
    with open(output_file, 'w') as f_out:
        for input_file in input_files:
            with open(input_file, 'r') as f_in:
                for line in f_in:
                    f_out.write(line)

input_files = ['nifH_01052024_LH.fna', 'nifD_01052024_LH.fna', 'nifK_01052024_LH.fna', 'ChlN_01052024_LH.fna', 'ChlB_01052024_LH.fna']
output_file = 'nifHDK_ChlNB.fasta'
combine_fasta_files(input_files, output_file)