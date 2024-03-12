import pandas as pd
import os
import csv
import re
import io
import argparse

  
def process_row(breakpointNr, breakpoint, gene_plus, gene_minus, length_of_region, data, chromosome):
    position = breakpoint['position']
    breakpoint_plus, breakpoint_minus = position, length_of_region - position

    if not gene_plus.empty:
        classified_position_plus = classify_position(breakpoint_plus, gene_plus, "plus")
    else: 
        classified_position_plus = [breakpoint_plus, False, None, None, None, "NoAnnotatedRegion", None]
    if not gene_minus.empty:
        classified_position_minus = classify_position(breakpoint_minus, gene_minus, "minus")
    else: 
        classified_position_minus = [breakpoint_minus, False, None, None, None, "NoAnnotatedRegion", None]

    match_on_plus_strand = classified_position_plus[1]
    match_on_minus_strand = classified_position_minus[1]

    # Dictionary to map the classification
    classifications = {
        (True, True): "plus_and_minus",
        (True, False): "only_plus",
        (False, True): "only_minus",
        (False, False): "no_match"
    }

    # Append data with appropriate classification
    classification = classifications[(match_on_plus_strand, match_on_minus_strand)]
    data.append([breakpointNr, chromosome, classification, "plus"] + classified_position_plus)
    data.append([breakpointNr, chromosome, classification, "minus"] + classified_position_minus)



def classify_position(position_to_check, gff_df, strand):
    # Filter mRNA and genes
    matching_mrnas = filter_dataframe_by_position(gff_df[gff_df['feature'].isin(['mRNA', 'transcript'])], position_to_check)
    matching_genes = filter_dataframe_by_position(gff_df[gff_df['feature'] == 'gene'], position_to_check)

    # Initialize result list
    result_list = [position_to_check, False, None, None, None, None, None]

    # Check if matching feature exists
    matching_feature = matching_mrnas.iloc[0] if not matching_mrnas.empty else matching_genes.iloc[0] if not matching_genes.empty else None

    if matching_feature is not None:
        attributes = parse_attribute(matching_feature['attribute'])
        result_list[1:5] = [True, matching_feature['feature'], attributes.get('ID', None), matching_feature['source']]
        tss = matching_feature['start']
        distance = tss - position_to_check
        result_list[6] = distance

        # Check if position is inside an exon
        exons_df = gff_df[(gff_df['feature'] == 'exon') & (gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
        if not exons_df.empty:
            is_inside_exon_df = filter_dataframe_by_position(exons_df, position_to_check)
            if not is_inside_exon_df.empty:
                result_list[5] = 'Exon'
                cds_df = gff_df[(gff_df['feature'] == 'CDS') & (gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
                if not filter_dataframe_by_position(cds_df, position_to_check).empty:
                    result_list[5] = 'ExonCDS'
            else:
                # Determine position relative to gene
                start_gene = matching_genes['start'].min()
                end_gene = matching_genes['end'].max()
                result_list[5] = search_intron_structure(exons_df, position_to_check, start_gene, end_gene, True)
        else:
            # Handle when there are no exons
            feature_parent_df = gff_df[(gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
            if feature_parent_df.empty: 
                pattern_id = f'ID={attributes.get("ID", "")};'
                feature_df = gff_df[(gff_df['attribute'].str.contains(pattern_id))]
                result_list = classify_regions_no_mrna(feature_df, matching_feature, position_to_check, result_list, 'inGeneButNotInMrna', 'Error', True)
            else: 
                result_list = classify_regions_no_mrna(feature_parent_df, matching_feature, position_to_check, result_list, 'inGeneButNotInMrna', 'Error', True)
    else: 
        # Handle when there are no matching features
        other_regions_df = gff_df[~gff_df['feature'].isin(['gene', 'mRNA', 'exon', 'transcript'])]
        classified_region_list = classify_regions_no_mrna(other_regions_df, matching_feature, position_to_check, result_list, 'InOtherRegion', 'NoAnnotatedRegion', False)
        if classified_region_list[2] == "match": 
            classified_region_list[1] = False
            classified_region_list[5] = "NoAnnotatedRegion"

        result_list = classified_region_list
    
    return result_list


def filter_dataframe_by_position(dataframe, position):
    return dataframe[(dataframe['start'] <= position) & (dataframe['end'] >= position)]

def parse_attribute(attribute):
    attribute_parts = attribute.split(';')
    return {part.split('=')[0]: part.split('=')[1] for part in attribute_parts}

def search_intron_structure(df, position, parent_df_start, parent_df_end, is_mrna):
    if any((end1 < position) and (start2 > position) for end1, start2 in zip(df['end'].shift(0), df['start'].shift(-1))):
        if is_mrna:
            return 'Intron'
        else:
            return 'InGeneBetweenTwoRegions'
    elif (parent_df_start < position) and (df['start'].min() > position):
        if is_mrna:
            return '5UTR'
        else:
            return 'StartOfGene'
    elif (parent_df_end > position) and (df['end'].max() < position):
        if is_mrna:
            return '3UTR'
        else: 
            return 'EndOfGene'
    else:
        return 'Error'

def classify_regions_no_mrna(df, parent_df, position, classified_list, matchLabel, noMatchLabel, isGene): 
    if not df.empty:
        matching_df = filter_dataframe_by_position(df, position)
        if not matching_df.empty:
            matching_fature_row = matching_df.iloc[0]
            matching_fature_id = matching_fature_row['attribute'].split(';')[0].split('=')[1]
            matching_fature_db = matching_fature_row['source']
            classified_list[1:5] = [True, matching_fature_row['feature'], matching_fature_id, matching_fature_db]
            classified_list[5] = matchLabel
        else:
            if isGene: 
                matching_fature_parent_row = parent_df.iloc[0]
                matching_fature_parent_id = matching_fature_parent_row['attribute'].split(';')[0].split('=')[1]
                matching_fature_parent_db = matching_fature_parent_row['source']
                classified_list[1:5] = [True, matching_fature_parent_row['feature'], matching_fature_parent_id, matching_fature_parent_db]
                parent_df_start = parent_df['start']
                parent_df_end = parent_df['end']
                classified_list[5] = search_intron_structure(df, position, parent_df_start, parent_df_end, False)
            else: 
                classified_list[5] = noMatchLabel
    else:
        classified_list[5] = noMatchLabel

    return classified_list
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process BLAST results and annotation data")
    parser.add_argument("positions_file_path", help="Path to file with genomic positions")
    parser.add_argument("annotation_file_path", help="Path to the annotation file")
    parser.add_argument("output_path", help="Path to save the CSV and plots")
    args = parser.parse_args()

    annotation_header = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    annotation_lines = re.findall(r'^(?!#).*', open(args.annotation_file_path, 'r').read(), re.MULTILINE)
    annotation_df = pd.read_csv(io.StringIO('\n'.join(annotation_lines)), sep='\t', header=None, names=annotation_header)

    column_names = ["genomicID", "coverage", "position"]
    df = pd.read_csv(args.positions_file_path, sep='\t', header=None, names=column_names)

    data = []
    for index, breakpoint in df.iterrows():
        print(breakpoint)
        chromosome = breakpoint['genomicID']

        matching_rows = annotation_df[annotation_df["seqname"] == chromosome].sort_values(by=["start"])
        matching_rows_region = matching_rows[matching_rows["feature"] == "region"]

        length_of_region = int(matching_rows_region["end"].iloc[0])
    
        matching_rows_gene = matching_rows[matching_rows["feature"] != "region"].sort_values(by=["start"])
        matching_rows_gene_plus = matching_rows_gene[matching_rows_gene["strand"] == "+"].sort_values(by=["start"])
        matching_rows_gene_minus = matching_rows_gene[matching_rows_gene["strand"] == "-"].sort_values(by=["start"])
      
        process_row(index, breakpoint, matching_rows_gene_plus, matching_rows_gene_minus, length_of_region, data, chromosome)

    positions_filename = os.path.basename(args.positions_file_path)
    csv_file_path = os.path.join(args.output_path, os.path.splitext(positions_filename)[0] + '_features.csv')
    data_header = ["BreakpointNr", "Chromosome", "StrandMatch_info", "Strand", "Position", "AnnotationFound(True_False)", "Feature", "Feature_ID", "Fature_DB", "Classification", "DistanceToTSS"]

    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(data_header)
        csv_writer.writerows(data)
