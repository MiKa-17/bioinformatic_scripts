import pysam
import argparse
import os

def process_sam_file(input_sam_file, output_sam_file, chromosome, start_position, end_position):
    count_reads = 0
    read_ids = set()

    with pysam.AlignmentFile(input_sam_file, "r") as input_sam:
        with pysam.AlignmentFile(output_sam_file, "w", header=input_sam.header) as output_sam:
            for read in input_sam:
                if read.reference_name == chromosome and \
                   start_position <= read.reference_start <= start_position + 10 and \
                   read.reference_end == end_position:
                    count_reads += 1
                    output_sam.write(read)
                    read_ids.add(read.query_name)
    
    return count_reads, read_ids

def write_read_ids(ids_output_file, read_ids):
    with open(ids_output_file, "w") as output_file:
        for read_id in read_ids:
            output_file.write(read_id + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract reads from a SAM file based on specified criteria.")
    parser.add_argument("input_sam_file", help="Input SAM file")
    parser.add_argument("output_sam_file", help="Output SAM file")
    parser.add_argument("chromosome", help="Chromosome name")
    parser.add_argument("start_position", type=int, help="Start position")
    parser.add_argument("end_position", type=int, help="End position")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_sam_file):
        count_reads, read_ids = process_sam_file(args.input_sam_file, args.output_sam_file, args.chromosome, args.start_position, args.end_position)
        
        ids_output_file = args.output_sam_file.replace(".sam", "_ids.txt")
        
        if not os.path.exists(ids_output_file):
            write_read_ids(ids_output_file, read_ids)
        
        print(f"Number of reads mapping only to the region: {count_reads}")
    else:
        print(f"Output file '{args.output_sam_file}' already exists. Skipping processing.")

if __name__ == "__main__":
    main()


"""
import pysam

# Define your input SAM file
input_sam_file = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/3_mapping_to_plasmid_and_genome/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_and_genome.sam"
output_sam_file = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/4_mapping_to_plasmid_genome_insert_region/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_and_genome_insertRegionTTAA.sam"


# Define the specific genomic region you want to extract reads from
chromosome = "PBGGPEx2-0m_TD4_BS14_HC_9418bp"
start_position = 2633  # Replace with your specific start position
end_position = 2996    # Replace with your specific end position

# Initialize a counter for the reads that meet your criteria
count_reads = 0
input_sam = pysam.AlignmentFile(input_sam_file, "r")
output_sam = pysam.AlignmentFile(output_sam_file, "w", header=input_sam.header)
# Initialize a set to store unique read IDs
read_ids = set()

for read in input_sam:
        # Check if the read maps to the specified chromosome
        if read.reference_name == chromosome:
            # Check if the read's mapping position is within the specified region
            if start_position <= read.reference_start <= start_position + 10 and read.reference_end == end_position:
                # Increment the counter for reads that meet the criteria
                count_reads += 1
                output_sam.write(read)
                read_ids.add(read.query_name)

# Print the count of reads that map only to the specific region
print(f"Number of reads mapping only to the region: {count_reads}")
# Close both files
input_sam.close()
output_sam.close()

ids_output_file = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/4_mapping_to_plasmid_genome_insert_region/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_and_genome_insertRegionTTAA_ids.txt"
# Write the unique read IDs to a text file
with open(ids_output_file, "w") as output_file:
    for read_id in read_ids:
        output_file.write(read_id + "\n")

"""