import argparse
import os
import pysam
from collections import Counter

'''
This Python script is designed to process SAM/BAM files (sequence alignment/map files) using the PySAM library, specifically focusing on extracting information related to the Supplementary Alignments (SA) tag. The SA tag contains information about secondary alignments of a read in the SAM/BAM file.
Printing how many reads are mapping to the plasmid and also to the genome (according to the supplementary alignment information). 
Saves the uniq read IDs into an output text file. 

'''


def process_sam_file(input_file):
    sa_values = []
    unique_ids = set()
    unique_ids_all = set()

    # Open the input BAM file using pysam
    with pysam.AlignmentFile(input_file, "rb") as sam_file:
        for read in sam_file:
            # Extract read ID
            read_id_all = read.query_name
            unique_ids_all.add(read_id_all)
            
            try:
                # Try to get the SA tag
                #tag = "SA:Z:(Z|[1-9]|1[0-9]|2[0-9])" # Tag for cairinia moschate ensemble genome
                tag = "SA:Z:" #Tag for human genome
                sa_field = read.get_tag(tag)
                if "N" in sa_field:
                    unique_ids.add(read_id_all)

                    # Process SA tag
                    substrings = sa_field.split(";")
                    for substring in substrings:
                        elements = substring.split(",")
                        if elements and all(element.strip() for element in elements):
                            sa_values.append(elements[0])

            except KeyError:
                # SA tag not found, continue to the next read
                pass

    return sa_values, unique_ids, unique_ids_all

def count_and_print_sa_values(sa_values, unique_ids, unique_ids_all):
    # Convert sa_values list to a set for uniqueness, then back to a list, and sort it
    unique_sa_values = sorted(list(set(sa_values)))

    print(f"Number of reads that are mapping on the plasmid and also on the genome: {len(unique_ids)} of {len(unique_ids_all)}\n")

    # Count occurrences of unique SA values
    sa_counts = Counter(sa_values)

    for sa_value in unique_sa_values:
        print(f"{sa_value} {sa_counts[sa_value]}")

def save_unique_ids(unique_ids, output_ids_file):
    # Save unique IDs to a file if an output file is specified
    if output_ids_file:
        if not os.path.exists(output_ids_file):
            with open(output_ids_file, "w") as id_file:
                id_file.write("\n".join(unique_ids))
        else:
            print(f"Error: Output file '{output_ids_file}' already exists.")

def main():
    parser = argparse.ArgumentParser(description="Process SAM/BAM file and extract SA values")
    parser.add_argument("input_file", help="Input SAM/BAM file")
    parser.add_argument("output_ids_file", nargs='?', default=None, help="Output file for unique IDs")

    args = parser.parse_args()

    sa_values, unique_ids, unique_ids_all = process_sam_file(args.input_file)
    count_and_print_sa_values(sa_values, unique_ids, unique_ids_all)
    save_unique_ids(unique_ids, args.output_ids_file)

if __name__ == "__main__":
    main()



"""

import sys
from collections import Counter
import re
import os

if len(sys.argv) != 3:
    print("Usage: python extract_sa_values.py input.sam output_ids.txt")
    sys.exit(1)


input_file = sys.argv[1]
output_ids_file = sys.argv[2]
sa_values = []
supplimentAligmentID = "SA:Z:"
unique_ids = set() 
unique_ids_all = set()


with open(input_file, "r") as sam_file:
    for line in sam_file:
            header = re.search(r'^@', line)
            if not header:
                columns_all = line.strip().split('\t')
                read_id_all = columns_all[0]
                unique_ids_all.add(read_id_all)
                print(read_id_all)
            sa_match = re.search(r'SA:Z:.*N*\t', line)
            if sa_match:
                columns = line.strip().split('\t')
                read_id = columns[0]
                unique_ids.add(read_id)
                sa_field = sa_match.group(0)
                substrings = sa_field.split("SA:Z:")[1].split(";")
                for substring in substrings:
                    elements = substring.split(",")
                    if elements and all(element.strip() for element in elements):
                        sa_values.append(elements[0])


unique_sa_values = list(set(sa_values))
sa_counts = Counter(sa_values)
unique_sa_values.sort()

print(f"Number of reads that are mapping on the plasmid and also on the genome: {len(unique_ids)} of {len(unique_ids_all)}\n")
for sa_value in unique_sa_values:
    print(f"{sa_value} {sa_counts[sa_value]}")

# Save unique IDs to a file

if not os.path.exists(output_ids_file):
    with open(output_ids_file, "w") as id_file:
        id_file.write("\n".join(unique_ids))
else:
    print(f"Error: Output file '{output_ids_file}' already exists.")
"""
