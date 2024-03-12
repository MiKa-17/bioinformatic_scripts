import sys
import os
import pysam
def sam_to_gff(sam_file):
    # Extract base name (excluding extension) and directory from input file path
    base_name = os.path.splitext(os.path.basename(sam_file))[0]
    output_file = os.path.join(os.path.dirname(sam_file), f"{base_name}_annotation.gff3")

    with pysam.AlignmentFile(sam_file, "r") as samfile, open(output_file, "w") as gff3file:
        for alignment in samfile:
            # Extract relevant information
            seqname = alignment.reference_name
            source = "SAMtoGFF3"
            feature_type = "alignment"
            start = alignment.reference_start + 1  # GFF3 is 1-based
            end = alignment.reference_end
            score = alignment.mapping_quality
            strand = "+" if not alignment.is_reverse else "-"

            # Format the GFF3 entry
            gff3_entry = "\t".join([seqname, source, feature_type, str(start), str(end),
                                    str(score), strand, ".", f'ID={alignment.query_name};gene={alignment.query_name}'])

            # Write to the GFF3 file
            gff3file.write(gff3_entry + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sam_to_gff.py input.sam")
        sys.exit(1)

    sam_file = sys.argv[1]
    sam_to_gff(sam_file)
