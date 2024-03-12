import argparse
from Bio import SeqIO

def calculate_strand(start, end):
    """Calculate the strand based on the difference between end and start positions."""
    return '+' if int(end) - int(start) >= 0 else '-'


def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Convert BLAST outfmt 6 to GFF3 format.")
    parser.add_argument("blast_file", help="Path to the BLAST outfmt 6 file")
    parser.add_argument("output_gff3", help="Path to the output GFF3 file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Process BLAST results and create GFF3 entries
    gff3_entries = []
    query_name_count = {}  # To keep track of query_name occurrences

    with open(args.blast_file, "r") as blast_file:
        for line in blast_file:
            fields = line.strip().split("\t")
            # default blast header:
            #"query_name", "subject_name", "percent_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"
            query_name, subject_name, _, _, _, _, _, _, start, end, _, _ = fields
            # Keep track of query_name occurrences
            query_name_count[query_name] = query_name_count.get(query_name, 0) + 1            
            query_name_unique = f"{query_name}_{query_name_count[query_name]}"

            # Calculate the strand based on start and end positions
            strand = calculate_strand(start, end)

            # Create GFF3 entry
            gff3_entry = f"{subject_name}\t.\tblastMatch\t{start}\t{end}\t.\t{strand}\t.\tID={query_name_unique}"
            gff3_entries.append(gff3_entry)

    # Write GFF3 entries to a file
    with open(args.output_gff3, "w") as gff3_file:
        gff3_file.write("##gff-version 3\n")
        gff3_file.write("\n".join(gff3_entries))

if __name__ == "__main__":
    main()
