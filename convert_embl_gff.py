import argparse
import os
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="Convert .embl file to .gff")
parser.add_argument("input_file_path", help="Path to the embl file")
parser.add_argument("reference_genome_name", help="name of the corresponding genomic sequence")
args = parser.parse_args()

# Define input file name and output file name with the same path but different extension
input_file_path = args.input_file_path
output_file_path = os.path.splitext(input_file_path)[0] + ".gff3"

# Specify the reference genome name
reference_genome_name = args.reference_genome_name


# Define a regular expression pattern to match positions like "complement(join(7694..8417,1..92))"
position_pattern = re.compile(r'join\{\[(\d+)\:(\d+)\]\(\-\), \[(\d+)\:(\d+)\]\(-\)\}')


# Function to split and process a complex position
def split_and_process_position(position_str, attributes):

    position_string = str(position_str)
    if re.search(position_pattern, str(position_string)):

        # Extract the join components and split them
        join_positions = re.findall(r'\[.*?\]', position_string)
        for join_position in join_positions:
            start, end = map(int, join_position[join_position.find('[')+1:join_position.find(']')].split(':'))
            output_file.write(f"{reference_genome_name}\t.\t{feature_type}\t{start + 1}\t{end}\t.\t.\t.\t{attributes}\n")
    else:
        # Handle regular positions
        start = feature.location.start + 1  # Adjust for 1-based indexing
        end = feature.location.end
        output_file.write(f"{reference_genome_name}\t.\t{feature_type}\t{start}\t{end}\t.\t.\t.\t{attributes}\n")

# Open the .embl file for reading and the .gff3 file for writing
with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
    # Write the reference genome name as a comment in the GFF3 file
    output_file.write(f"##genome-build {reference_genome_name}\n")
    
    # Iterate through .embl records
    for record in SeqIO.parse(input_file, "embl"):
        for feature in record.features:
            feature_type = feature.type
            
            # Extract the position and attributes
            position_str = feature.location
    # Extract the "gene" attribute, or if it doesn't exist, use the "label" attribute
            gene_attribute = feature.qualifiers.get("gene", feature.qualifiers.get("label", [""]))[0]
            
            # Replace spaces with underscores in the gene name
            gene_attribute = gene_attribute.replace(' ', '_')
            
            # Combine all other attributes
            other_attributes = [f"{key}={value}" for key, value in feature.qualifiers.items() if key != "gene" and key != "label"]
            attributes = f"gene={gene_attribute};{';'.join(other_attributes)}"
            
            
            # Process the position based on its format
            split_and_process_position(position_str, attributes)

print(f"Conversion complete. Output written to {output_file_path}")