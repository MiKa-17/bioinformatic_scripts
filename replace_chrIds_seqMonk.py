mapping_file = '3_mapping/3_reads_on_plasmid_and_genome/mapping_ids.txt'
input_sam_file = '3_mapping/3_reads_on_plasmid_and_genome/LvthmInsertReads_trimmed_l750_on_plasmid_and_genome_seqMonk.sam'
output_sam_file = '3_mapping/3_reads_on_plasmid_and_genome/LvthmInsertReads_trimmed_l750_on_plasmid_and_genome_seqMonk_newIds.sam'

# Read mapping file
chromosome_mapping = {}
with open(mapping_file, 'r') as f:
    for line in f:
        original, updated = line.strip().split('\t')
        chromosome_mapping[original] = updated

# Process SAM file
with open(input_sam_file, 'r') as f_in, open(output_sam_file, 'w') as f_out:
    for line in f_in:
        if line.startswith('@'):
            # Write header lines unchanged
            f_out.write(line)
        else:
            # Process alignment lines
            parts = line.split('\t')
            chromosome = parts[2]
            if chromosome in chromosome_mapping:
                # Replace chromosome identifier
                parts[2] = chromosome_mapping[chromosome]
            # Write modified line to output SAM file
            f_out.write('\t'.join(parts))

print("SAM file processed. New SAM file saved as:", output_sam_file)
