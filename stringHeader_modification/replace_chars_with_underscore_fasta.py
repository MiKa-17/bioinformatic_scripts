import os
import re
import argparse

# Function to replace spaces, special characters, and consecutive underscores with a single underscore
def replace_chars_with_underscore(header):
    header = re.sub(r' ', '_', header)
    header = re.sub(r'[^A-Za-z0-9\>\n]', '_', header)
    header = re.sub(r'[_]+', '_', header)
    return header


# Create an ArgumentParser to handle command-line arguments
parser = argparse.ArgumentParser(description="Replace spaces and special characters with underscores in FASTA files.")
parser.add_argument("input_dir", help="Input directory containing FASTA files")

# Parse the command-line arguments
args = parser.parse_args()

# Directory containing the FASTA files
directory = args.input_dir

# Iterate through the files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.fa'):
        filepath = os.path.join(directory, filename)
        modified_lines = []
        with open(filepath, 'r') as file:
            lines = file.readlines()
            
            for line in lines:
                if line.startswith(">"):  # Preserve the ">" character at the start of the header
                    modified_header = replace_chars_with_underscore(line)
                    modified_lines.append(modified_header)
                else:
                    modified_lines.append(line)

        
        with open(filepath, 'w') as file:
            file.writelines(modified_lines)

        print(f"Processed: {filename}")
