import pysam
from Bio import SeqIO
import sys 
import re

# Replace 'input.sam' with the path to your SAM file

if len(sys.argv) != 3:
    print("Usage: python parse_sam_to_r_whole_genome_plot.py input.sam input.fasta")
    sys.exit(1)

input_sam_file = sys.argv[1]
input_fasta_file = sys.argv[2] #genome reference
count_threshold = 0 

# Initialize a dictionary to store counts of seqIds
seqId_counts = {}

# Open the SAM file for reading
with pysam.AlignmentFile(input_sam_file, 'r') as samfile:
    for read in samfile:
        # Extract the third column from the SAM file (RNAME)
        seqId = read.reference_name
        
        # Check if seqId starts with "N"
        if seqId != None:
            if seqId.startswith("NC"):
                # Count the occurrences of each seqId
                if seqId in seqId_counts:
                    seqId_counts[seqId] += 1
                else:
                    seqId_counts[seqId] = 1

# Extract seqIds with counts > 100
selected_seqIds = [seqId for seqId, count in seqId_counts.items() if count > count_threshold or seqId.startswith("NC")]
selected_seqIds.sort()  #Sort the selected_seqIds alphabetically
print(selected_seqIds)
# Format the selected seqIds as a comma-separated string
formatted_seqIds = ', '.join(['"' + seqId + '"' for seqId in selected_seqIds])
#formatted_seqIds_for_grep_names = '|'.join(selected_seqIds)

# Initialize a dictionary to store sequence lengths
seq_lengths = {}
seq_headers = {}

# Open the FASTA file for reading
with pysam.FastaFile(input_fasta_file) as fastaf:
    for seqId in selected_seqIds:
        # Get the length of the genomic sequence
        seq_length = len(fastaf.fetch(seqId))
        seq_lengths[seqId] = seq_length

# Open the FASTA file for reading using Biopython's SeqIO
with open(input_fasta_file, 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        recordID = record.id
        full_header = record.description
        
        # Check if the current seqId is in the selected_seqIds
        if recordID in selected_seqIds:
            seq_headers[recordID] = full_header

# Print the full headers
print("Full Headers:")
for seqId, full_header in seq_headers.items():
    print(f"{seqId}: {full_header}")

# Initialize a list to store shortened headers
shortened_headers = []

# Shorten the headers based on rules
for seqId in selected_seqIds:
    full_header = seq_headers[seqId]
    # Initialize the shortened header
    shortened_header = seqId
    
    # Shorten based on "chromosome" and a number
    if "chromosome" in full_header:

        match = re.search(r'chromosome (\d+|X|Y)', full_header)
        if match:
            chromosome_number = match.group(1)
            shortened_header = f"chr_{chromosome_number}"
    elif "scaffold" in full_header:
        match = re.search(r'scaffold_(\d+)', full_header)
        if match:
            number = match.group(1)
            shortened_header = f"scaf_{number}"
    

    #if "chr1" in full_header:
    #    match = re.search(r'chr1_(\d+)', full_header)
    #    if match:
    #        shortened_header = match.group()

    #if "unplaced_scaffold" in full_header: 
    #    match = re.search(r'unplaced_scaffold_(\d+)', full_header)
    #    if match:
    #        number = match.group(1)
    #        shortened_header = f"unScaf_{number}"

    shortened_headers.extend([f'"{shortened_header}"'])

# Combine the shortened headers into the desired format
formatted_shortened_headers = ', '.join(shortened_headers)

# Format the sequence lengths as a comma-separated string
formatted_seq_lengths = ', '.join(map(str, [seq_lengths[seqId] for seqId in selected_seqIds]))

isCircular = {}
for seqId in selected_seqIds:
    isCircular[seqId] = "FALSE"
formatted_isCircular = ', '.join(map(str, [isCircular[seqId] for seqId in selected_seqIds]))

# Print the formatted string
print(f"seqIds <- c(\n   {formatted_seqIds}\n)")

print(f"seqLen <- c(\n  {formatted_seq_lengths}\n)")

print(f"isCircular <- c(\n  {formatted_isCircular}\n)")
# Print the formatted headers
print(f"c(\n{formatted_shortened_headers}\n)")

#print(formatted_seqIds_for_grep_names)