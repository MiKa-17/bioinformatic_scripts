from collections import Counter
import pysam
import os

# Replace these with your file paths
sam_file = "/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/3_mapping/6_reads_on_flyeCons/reads_cr_trimmed_l600_onChr13_on_flyeCons.sam"
reference_file = "/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/7_assembly/10-consensus/consensus.fasta"
output_file = "/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/consSeq_flyCons_reads_optimized.fasta"

# Load the reference genome
with open(reference_file, "r") as reference:
    reference_sequence = "".join(line.strip() for line in reference if not line.startswith(">"))

# Initialize a list to store the consensus sequence
consensus_sequence = list(reference_sequence)

# Open the SAM file and process aligned reads
with pysam.AlignmentFile(sam_file, "r") as sam:
    for read in sam:
        if not read.is_unmapped:
            for i, base in enumerate(read.query_sequence):
                reference_position = read.reference_start + i
                # Update the consensus sequence with the majority base at each position
                if reference_position >= 0 and reference_position < len(consensus_sequence):
                    consensus_sequence[reference_position] = max("ACGTN", key=lambda x: base.count(x))

# Convert the consensus sequence back to a string
consensus_sequence = "".join(consensus_sequence)

# Write the consensus sequence to a FASTA file
with open(output_file, "w") as output:
    output.write(">Consensus\n")
    output.write(consensus_sequence + "\n")