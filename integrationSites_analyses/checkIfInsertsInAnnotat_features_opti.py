import pandas as pd
import os
import csv
import re
import io
import argparse

def process_row(row, matching_rows_gene_plus, matching_rows_gene_minus, length_of_region, data, chromosome):

    start_str, end_str = row['Positions'].split("-")
    start1 = int(start_str)
    end1 = int(end_str)

    breakpoint_start_plus, breakpoint_start_minus = start1, length_of_region - start1
    breakpoint_end_plus, breakpoint_end_minus = end1, length_of_region - end1

    classied_position_plus = classify_position((breakpoint_start_plus, breakpoint_end_plus), matching_rows_gene_plus)
    classied_position_minus = classify_position((breakpoint_start_minus,breakpoint_end_minus), matching_rows_gene_minus)


    data.append([chromosome, "plus"] + classied_position_plus)
    data.append([chromosome, "minus"] + classied_position_minus)

def filter_dataframe_by_position(dataframe, position):
    #return dataframe[(dataframe['start'] <= position) & (dataframe['end'] >= position)]
    start, end = position
    return dataframe[(dataframe['start'] >= start) & (dataframe['end'] <= end) | (dataframe['start'] <= start) & (dataframe['end'] >= end)]

def classify_position(position_to_check, gff_df):
    # Filter rows with 'gene' or 'mRNA' features

    gene_mrna_df = gff_df[gff_df['feature'].isin(["gene", 'mRNA', 'transcript'])]
    gene_df = gff_df[gff_df['feature'].isin(['gene'])]

    # Check if the position is within a 'mRNA' or 'gene'
    mrna_transcript_df = gff_df[gff_df['feature'].isin(['mRNA', 'transcript'])]

    matching_mrnas = filter_dataframe_by_position(mrna_transcript_df, position_to_check)
    matching_genes = filter_dataframe_by_position(gene_df, position_to_check)

    # results list: "Position", "AnnotationFound (True/False)", "Feature", "Feature_ID", "Fature_DB", "Classification", "DistanceToTSS"
    result_list = [position_to_check, False, None, None, None, None, None]

    matching_feature = None
    if not matching_mrnas.empty:
        matching_feature = matching_mrnas.iloc[0]
    elif not matching_genes.empty:
        matching_feature = matching_genes.iloc[0]

    if matching_feature is not None:
        matching_feature_id = matching_feature['attribute'].split(';')[0].split('=')[1]
        matching_feature_db = matching_feature['source']
        result_list[1:5] = [True, matching_feature['feature'], matching_feature_id, matching_feature_db]



    if result_list[1]:

        #distance to TSS
        matching = None
        if not matching_genes.empty:
            matching = matching_genes.iloc[0]
        else: 
            matching = matching_feature
        tss = matching['start']
        position_start, position_end = position_to_check
        distance = tss-position_start
        result_list[6] = distance

        # Check if the feature has annotated exons
        exons_df = gff_df[(gff_df['feature'] == 'exon') & (gff_df['attribute'].str.contains(f'Parent={matching_feature_id}'))]

        if not exons_df.empty:
            # Check if the position is within an exon
            is_inside_exon_df = filter_dataframe_by_position(exons_df, position_to_check)

            if not is_inside_exon_df.empty:
                result_list[5] = 'InsideExon'
                cds_df = gff_df[(gff_df['feature'] == 'CDS') & (gff_df['attribute'].str.contains(f'Parent={matching_feature_id}'))]
                inside_cds_df = filter_dataframe_by_position(cds_df, position_to_check)
                if not inside_cds_df.empty:
                    result_list[5] = 'InsideExonCDS'
            else:
                between_exons = any((end1 < position_start) and (start2 > position_end) for end1, start2 in zip(exons_df['end'].shift(0), exons_df['start'].shift(-1)))
                if between_exons:
                    result_list[5] = 'Intron'
                else:
                    # Check if it is between the start of a gene and an exon (5' UTR)
                    is_5utr = any((start1 < position_start) and (start2 > position_start) for start1, start2 in zip(gene_mrna_df['start'].shift(0), exons_df['start'].shift(-1)))
                    
                    if is_5utr:
                        result_list[5] = '5UTR'
                    else:
                        # Check if it is between an exon and the end of a gene (3' UTR)
                        is_3utr = any((end1 > position_end) and (end2 < position_end) for end1, end2 in zip(gene_mrna_df['end'].shift(0), exons_df['end'].shift(-1)))

                        if is_3utr:
                            result_list[5] = '3UTR'
                        else:
                            result_list[5] = 'Unclassified'
        else:
            result_list[5] = 'inGeneButNotInMrna'
            feature_parent_df = gff_df[(gff_df['attribute'].str.contains(f'Parent={matching_feature_id}'))]

            if not feature_parent_df.empty:
                matching_feature_parent = filter_dataframe_by_position(feature_parent_df, position_to_check)

                if not matching_feature_parent.empty:
                    matching_feature_parent = matching_feature_parent.iloc[0]
                    matching_feature_parent_id = matching_feature_parent['attribute'].split(';')[0].split('=')[1]
                    matching_feature_parent_db = matching_feature_parent['source']
                    result_list[1:5] = [True, matching_feature_parent['feature'], matching_feature_parent_id, matching_feature_parent_db]
                    result_list[5] = 'inGeneButNotInMrna'
                else:
                    result_list[5] = '???'
            else:
                result_list[5] = '???'


    else: 
        other_regions_df = gff_df[~gff_df['feature'].isin(['gene', 'mRNA', 'exon', 'transcript'])]

        if not other_regions_df.empty:
            matching_feature = filter_dataframe_by_position(other_regions_df, position_to_check)

            if not matching_feature.empty:
                matching_feature = matching_feature.iloc[0]
                matching_feature_id = matching_feature['attribute'].split(';')[0].split('=')[1]
                matching_feature_db = matching_feature['source']
                result_list[1:5] = [True, matching_feature['feature'], matching_feature_id, matching_feature_db]
                result_list[5] = 'InOtherRegion'
            else:
                result_list[5] = 'NoAnnotatedRegion'
        else:
            result_list[5] = 'NoAnnotatedRegion'


    return result_list





if __name__ == '__main__':

    # Create an argument parser
    parser = argparse.ArgumentParser(description="Process BLAST results and annotation data")
    parser.add_argument("positions_file_path", help="Path to file with genomic positions")
    parser.add_argument("annotation_file_path", help="Path to the annotation file")
    parser.add_argument("output_path", help="Path to save the CSV and plots")

    # Parse command-line arguments
    args = parser.parse_args()

    annotation_header= ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    annotation_lines = re.findall(r'^(?!#).*', open(args.annotation_file_path, 'r').read(), re.MULTILINE)
    annotation_df = pd.read_csv(io.StringIO('\n'.join(annotation_lines)), sep='\t', header=None, names=annotation_header)

    # Define the column names for your BLAST output
    column_names = ["genomicID", "coverage", "Positions"]

    # Read the BLAST output into a DataFrame using the defined column names
    df = pd.read_csv(args.positions_file_path, sep=' ', header=None, names=column_names)

    # Create a list to store data for the CSV file
    data = []

    # Process the DataFrame in a single process
    for index, breakpoint in df.iterrows():
        print("----------------", index)
        chromosome = breakpoint['genomicID']
        #strand_blast = blast_hit['Subject_Strand']
        # Filter annotation DataFrame by matching chromosome
        matching_rows = annotation_df[annotation_df["seqname"] == chromosome].sort_values(by=["start"])

        matching_rows_region = matching_rows[matching_rows["feature"] == "region"]
        length_of_region = int(matching_rows_region["end"].iloc[0])

        
        matching_rows_gene = matching_rows[matching_rows["feature"] != "region"].sort_values(by=["start"])
        matching_rows_gene_plus = matching_rows_gene[matching_rows_gene["strand"] == "+"].sort_values(by=["start"])
        matching_rows_gene_minus = matching_rows_gene[matching_rows_gene["strand"] == "-"].sort_values(by=["start"])
        if not matching_rows_gene_plus.empty and not matching_rows_gene_minus.empty: 
            process_row(breakpoint, matching_rows_gene_plus, matching_rows_gene_minus, length_of_region, data, chromosome)


    positions_filename = os.path.basename(args.positions_file_path)

    
    # Save the data to a CSV file
        # Construct the CSV file path based on the output path and blast file name
    csv_file_path = os.path.join(args.output_path, os.path.splitext(positions_filename)[0] + '_features.csv')
    data_header = ["Chromosome", "Strand", "Position", "AnnotationFound (True/False)", "Feature", "Feature_ID", "Fature_DB", "Classification", "DistanceToTSS"]
    # Write the data to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(data_header)
        csv_writer.writerows(data)
