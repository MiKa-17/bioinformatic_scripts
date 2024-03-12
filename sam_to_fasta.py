import sys
import subprocess
import pandas as pd
'''
"sam_to_fasta.py" is a Python script designed to extract genomic sequences from a SAM file using coordinates provided in a BED (Browser Extensible Data) file. It processes each entry in the BED file, utilizing the samtools view command to extract relevant SAM records for specified genomic regions. The script then generates a single FASTA file, consolidating all extracted sequences and including informative headers for further analysis.
'''

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python process_sam.py input.sam input.bed output.fasta")
    sys.exit(1)

input_sam_file = sys.argv[1]
input_bed_file = sys.argv[2]
output_fasta_file = sys.argv[3]
bed_df = pd.read_csv(input_bed_file, sep='\t', header=None, names=['refName', 'start', 'stop'])
print(bed_df)
# Open the output FASTA file
with open(output_fasta_file, "w") as fasta_file:
    # Read each line from the BED file
    with open(input_bed_file, "r") as bed_file:
        for line in bed_file:
            line = line.strip()
            if not line:
                continue

            print(line.split(" "))
            # Split the line into reference name, start, and stop
            reference_name, start, stop = line.split("\t")
            start, stop = int(start), int(stop)

            # Use samtools view to extract relevant SAM records
            samtools_command = f"samtools view {input_sam_file} {reference_name}:{start}-{stop}"
            try:
                samtools_output = subprocess.check_output(samtools_command, shell=True, text=True)
            except subprocess.CalledProcessError:
                print("Error running samtools view for", reference_name, start, stop)
                continue
            header_counter = 1  # Initialize a counter for unique headers
            for sam_line in samtools_output.splitlines():
                fields = sam_line.strip().split("\t")
                read_name = fields[0]
                start_position = int(fields[3])
                sequence = fields[9]

                # Calculate the matching length
                matching_length = stop - start_position

                # Extract the sequence based on the matching length
                if matching_length > 0:
                    sequence = sequence[:matching_length]

                # Write the extracted sequence to the FASTA file
                fasta_file.write(f">{read_name}_{reference_name}_{start_position}-{stop}_{header_counter}\n")
                header_counter += 1
                fasta_file.write(sequence + "\n")

print(f"Processed data and saved to {output_fasta_file}")
