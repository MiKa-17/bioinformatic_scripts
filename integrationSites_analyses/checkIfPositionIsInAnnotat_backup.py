import pandas as pd
import os
import csv
import re
import io
import argparse

def process_row(breakpointNr, breakpoint, gene_plus, gene_minus, length_of_region, data, chromosome):
    position = breakpoint['position']
    #print(position)
    breakpoint_plus, breakpoint_minus = position, length_of_region - position


    classified_position_plus = classify_position(breakpoint_plus, gene_plus, "plus")
    classified_position_minus = classify_position(breakpoint_minus, gene_minus, "minus")

    if classified_position_plus[1] and classified_position_minus[1]:
        data.append([breakpointNr, chromosome, "plus_and_minus", "plus"] + classified_position_plus)
        data.append([breakpointNr, chromosome, "plus_and_minus", "minus"] + classified_position_minus)
    elif classified_position_plus[1] and not classified_position_minus[1]: 
        data.append([breakpointNr, chromosome, "only_plus", "plus"] + classified_position_plus)
        data.append([breakpointNr, chromosome, "only_plus", "minus"] + classified_position_minus)
    elif not classified_position_plus[1] and classified_position_minus[1]:
        data.append([breakpointNr, chromosome, "only_minus", "plus"] + classified_position_plus)
        data.append([breakpointNr, chromosome, "only_minus", "minus"] + classified_position_minus)
    else: 
        data.append([breakpointNr, chromosome, "no_match", "plus"] + classified_position_plus)
        data.append([breakpointNr, chromosome, "no_match", "minus"] + classified_position_minus)

def filter_dataframe_by_position(dataframe, position):
    ##print(position)
    ##print(dataframe)
    return dataframe[(dataframe['start'] <= position) & (dataframe['end'] >= position)]

def parse_attribute(attribute):
    attribute_parts = attribute.split(';')
    return {part.split('=')[0]: part.split('=')[1] for part in attribute_parts}

def classify_regions_no_mrna(df, parent_df, position, list, matchLabel, noMatchLabel, isGene): 
    classified_list = list

    print(df)
    print(parent_df)
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
                #child_df = df[(df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
                 ##print(df)
                 ##print(parent_df)
                 ##print(parent_df['start'])
                 ##print(df['start'].min())
                if any((end1 < position) and (start2 > position) for end1, start2 in zip(df['end'].shift(0), df['start'].shift(-1))):
                    classified_list[5] = 'Intron'
                elif (parent_df['start'] < position) and (df['start'].min() > position):
                    classified_list[5] = '5UTR'
                elif (parent_df['end'] > position) and (df['end'].max() < position):
                    classified_list[5] = '3UTR'
                else:
                    classified_list[5] = 'Error'
            else: 
                classified_list[5] = noMatchLabel
    else:
        classified_list[5] = noMatchLabel

    return classified_list
    

def classify_position(position_to_check, gff_df, strand):
    ##print(strand, position_to_check)
    gene_mrna_df = gff_df[gff_df['feature'].isin(["gene", 'mRNA', 'transcript'])]

    matching_mrnas = filter_dataframe_by_position(gff_df[gff_df['feature'].isin(['mRNA', 'transcript'])], position_to_check)
    matching_genes = filter_dataframe_by_position(gff_df[gff_df['feature'] == 'gene'], position_to_check)

    result_list = [position_to_check, False, None, None, None, None, None]

    matching_feature = matching_mrnas.iloc[0] if not matching_mrnas.empty else matching_genes.iloc[0] if not matching_genes.empty else None
    ##print(matching_feature)
    if matching_feature is not None:
        attributes = parse_attribute(matching_feature['attribute'])
        result_list[1:5] = [True, matching_feature['feature'], attributes.get('ID', None), matching_feature['source']]
        tss = matching_feature['start']
        distance = tss - position_to_check
        result_list[6] = distance

        exons_df = gff_df[(gff_df['feature'] == 'exon') & (gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]

        if not exons_df.empty:
            is_inside_exon_df = filter_dataframe_by_position(exons_df, position_to_check)
            #print(is_inside_exon_df)
            #print(matching_genes)
            #print(matching_genes['start'])
            if not is_inside_exon_df.empty:
                result_list[5] = 'Exon'
                cds_df = gff_df[(gff_df['feature'] == 'CDS') & (gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
                if not filter_dataframe_by_position(cds_df, position_to_check).empty:
                    result_list[5] = 'ExonCDS'
            else:
                start_gene = matching_genes['start'].min()
                end_gene = matching_genes['end'].max()
                start_firstExon = exons_df['start'].min()
                end_lastExon = exons_df['end'].max()

                if any((end1 < position_to_check) and (start2 > position_to_check) for end1, start2 in zip(exons_df['end'].shift(0), exons_df['start'].shift(-1))):
                    result_list[5] = 'Intron'
                elif (start_gene < position_to_check) and (start_firstExon > position_to_check):
                    #print("555555555555555555555555555555555")
                    #print(f'({start_gene} < {position_to_check}) and ({start_firstExon} > {position_to_check})')
                    result_list[5] = '5UTR'
                elif (end_gene > position_to_check) and (end_lastExon < position_to_check):
                    #print("3333333333333333333333333")
                    #print(f'({end_gene} > {position_to_check}) and ({end_lastExon} < {position_to_check})')
                    result_list[5] = '3UTR'
                else:
                    #print('ERROR')
                    result_list[5] = 'Error'
        else:
             # ##print("hey")
             # ##print(matching_feature['attribute'])
             # ##print(attributes.get("ID", ""))
            feature_parent_df = gff_df[(gff_df['attribute'].str.contains(f'Parent={attributes.get("ID", "")};'))]
             # ##print(feature_parent_df)
            if feature_parent_df.empty: 
                patternID = f'ID={attributes.get("ID", "")};'
                feature_df = gff_df[(gff_df['attribute'].str.contains(patternID))]
                 # ##print(feature_df)
                result_list = classify_regions_no_mrna(feature_df, matching_feature, position_to_check, result_list, 'inGeneButNotInMrna', 'Error', True)
            else: 
                 # ##print("hmmm")
                 # ##print(position_to_check)
                 # ##print(matching_feature)
                result_list = classify_regions_no_mrna(feature_parent_df, matching_feature, position_to_check, result_list, 'inGeneButNotInMrna', 'Error', True)

    else: 
        other_regions_df = gff_df[~gff_df['feature'].isin(['gene', 'mRNA', 'exon', 'transcript'])]
        classified_region_list = classify_regions_no_mrna(other_regions_df, matching_feature, position_to_check, result_list, 'InOtherRegion', 'NoAnnotatedRegion', False)
        if classified_region_list[2] == "match": 
            classified_region_list[1] = False
            classified_region_list[5] = "NoAnnotatedRegion"

        result_list = classified_region_list
    
    return result_list

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
        chromosome = breakpoint['genomicID']

        matching_rows = annotation_df[annotation_df["seqname"] == chromosome].sort_values(by=["start"])
        matching_rows_region = matching_rows[matching_rows["feature"] == "region"]

        length_of_region = int(matching_rows_region["end"].iloc[0])
    
        matching_rows_gene = matching_rows[matching_rows["feature"] != "region"].sort_values(by=["start"])
        matching_rows_gene_plus = matching_rows_gene[matching_rows_gene["strand"] == "+"].sort_values(by=["start"])
        matching_rows_gene_minus = matching_rows_gene[matching_rows_gene["strand"] == "-"].sort_values(by=["start"])

        if not matching_rows_gene_plus.empty and not matching_rows_gene_minus.empty:
            process_row(index, breakpoint, matching_rows_gene_plus, matching_rows_gene_minus, length_of_region, data, chromosome)

    positions_filename = os.path.basename(args.positions_file_path)
    csv_file_path = os.path.join(args.output_path, os.path.splitext(positions_filename)[0] + '_features.csv')
    data_header = ["BreakpointNr", "Chromosome", "StrandMatch_info", "Strand", "Position", "AnnotationFound(True_False)", "Feature", "Feature_ID", "Fature_DB", "Classification", "DistanceToTSS"]

    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(data_header)
        csv_writer.writerows(data)
