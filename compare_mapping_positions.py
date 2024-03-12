import os
import argparse
import subprocess
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import logging

def concatenate_regions_with_gaps(bed_dataframe):
    """
    Concatenate genomic regions with gaps in between.

    Args:
        bed_dataframe (DataFrame): DataFrame containing genomic regions.

    Returns:
        DataFrame: Concatenated genomic regions with gaps removed.
    """    
    # Initialize variables to store the result
    result = []

    # Initialize variables to track the current chromosome and range
    current_chrom = None
    current_start = None
    current_end = None

    # Iterate through rows in the DataFrame
    for index, row in bed_dataframe.iterrows():
        chrom = row['chrom']
        start = row['pos_start']
        end = row['pos_stop']

        # If it's a new chromosome or there's a gap, add the previous range to the result
        if current_chrom is None or current_chrom != chrom or start > current_end + 0:
            if current_chrom is not None:
                result.append([current_chrom, current_start, current_end])
            current_chrom = chrom
            current_start = start
        current_end = end

    # Add the last range to the result
    if current_chrom is not None:
        result.append([current_chrom, current_start, current_end])

    # Create a DataFrame from the result
    result_df = pd.DataFrame(result, columns=['chrom', 'pos_start', 'pos_end'])
    
    return result_df

def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Compare BAM covered areas between two datasets.")
    parser.add_argument("dataset1_bam", help="Path to the first BAM file.")
    parser.add_argument("dataset2_bam", help="Path to the second BAM file.")
    parser.add_argument("output_path", help="Output directory path.")
    return parser.parse_args()

def main():
    try: 
        # Parse command-line arguments
        args = parse_arguments()
        dataset1_bam = args.dataset1_bam
        dataset2_bam = args.dataset2_bam
        output_path = args.output_path

        # Configure logging
        logging.basicConfig(level=logging.INFO)

        # Define output file names using the specified output path
        dataset1_bed = os.path.join(output_path, "dataset1.bed")
        dataset2_bed = os.path.join(output_path, "dataset2.bed")
        common_coverage_bed = os.path.join(output_path, "common_coverage.bed")

        # Convert pileup to BED using bedtools
        logging.info("Convert pileup to BED using bedtools")
        subprocess.call(f"bedtools genomecov -ibam {dataset1_bam} -bg > {dataset1_bed}", shell=True)
        subprocess.call(f"bedtools genomecov -ibam {dataset2_bam} -bg > {dataset2_bed}", shell=True)

        # Find common coverage using bedtools
        logging.info("Find common coverage using bedtools")
        subprocess.call(f"bedtools intersect -a {dataset1_bed} -b {dataset2_bed} -u > {common_coverage_bed}", shell=True)

        logging.info("Analyze")
        dataset1 = pd.read_csv(dataset1_bed, sep='\t', header=None, names=['chrom', 'pos_start', 'pos_stop', 'coverage'])
        dataset2 = pd.read_csv(dataset2_bed, sep='\t', header=None, names=['chrom', 'pos_start', 'pos_stop', 'coverage'])
        common_coverage = pd.read_csv(common_coverage_bed, sep='\t', header=None, names=['chrom', 'pos_start', 'pos_stop', 'coverage'])
        

        concat_dataset1 = concatenate_regions_with_gaps(dataset1)
        concat_dataset2 = concatenate_regions_with_gaps(dataset2)
        common_coverage_concat = concatenate_regions_with_gaps(common_coverage)
        

        concat_dataset1_filtered = concat_dataset1[(concat_dataset1['chrom'] != "PBGGPEx2-0p_TD4_BS14_LC_8417bp") & concat_dataset1['chrom'] != "PBGGPEx2-0m_TD4_BS14_HC_9418bp"]

        #print(concat_dataset1.to_string())
        print(concat_dataset1['chrom'].value_counts())
        print(concat_dataset2['chrom'].value_counts())

        total_positions_dataset1 = len(concat_dataset1)
        total_positions_dataset2 = len(concat_dataset2)
        common_positions = len(common_coverage_concat)

        unique_positions_dataset1 = total_positions_dataset1 - common_positions
        unique_positions_dataset2 = total_positions_dataset2 - common_positions

        logging.info("Total positions in Dataset 1: %d", total_positions_dataset1)
        logging.info("Total positions in Dataset 2: %d", total_positions_dataset2)
        logging.info("Common positions: %d", common_positions)
        logging.info("Unique positions in Dataset 1: %d", unique_positions_dataset1)
        logging.info("Unique positions in Dataset 2: %d", unique_positions_dataset2)

        # Create a Venn diagram to visualize the overlap and unique positions
        venn2(subsets=(unique_positions_dataset1, unique_positions_dataset2, common_positions), 
            set_labels=('no add T', 'add T'))
        plt.title("Positions with read coverage in whole genome (TTAA)")
        plt.show()

    except Exception as e:
        logging.error("An error occurred: %s", str(e))

if __name__ == "__main__":
    main()
