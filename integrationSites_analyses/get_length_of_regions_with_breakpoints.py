import pandas as pd
import os
import csv
import re
import io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process BLAST results and annotation data")
    parser.add_argument("positions_file_path", help="Path to file with genomic positions")
    parser.add_argument("annotation_file_path", help="Path to the annotation file")
    args = parser.parse_args()

    annotation_header = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    annotation_lines = re.findall(r'^(?!#).*', open(args.annotation_file_path, 'r').read(), re.MULTILINE)
    annotation_df = pd.read_csv(io.StringIO('\n'.join(annotation_lines)), sep='\t', header=None, names=annotation_header)

    column_names = ["genomicID", "coverage", "position"]
    df = pd.read_csv(args.positions_file_path, sep='\t', header=None, names=column_names)

    data = []
    unique_genomic_regions = df['genomicID'].unique()

    for region in unique_genomic_regions:
        matching_rows = annotation_df[annotation_df["seqname"] == region].sort_values(by=["start"])
        matching_rows_region = matching_rows[matching_rows["feature"] == "region"]

        length_of_region = int(matching_rows_region["end"].iloc[0])
        data.append([region, length_of_region])


    positions_filename = os.path.basename(args.positions_file_path)
    csv_file_path = "4_breakpoints/1_random/genomic_regions_length.csv"

    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows(data)
