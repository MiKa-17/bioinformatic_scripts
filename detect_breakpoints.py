import argparse
import os
from collections import defaultdict
import pysam

def process_chromosome_ids(ids_file, sam_file):
    with open(ids_file, 'r') as f:
        chromosome_ids = [line.strip().split()[0] for line in f]

    print("ChromosomeID   NumberOfGroups NrOfGroupsMoreThan5Reads")

    results = {}
    for chromosome_id in chromosome_ids:
        positions = []
        samfile = pysam.AlignmentFile(sam_file, "r")

        for alignment in samfile:
            reference_id = samfile.get_reference_name(alignment.reference_id)
            if reference_id == chromosome_id:
                #handle matches on reverse strand 
                flag = alignment.flag
                alignment_start = alignment.reference_start
                if flag == "16" or flag == "272" or flag == "2064":
                    alignment_length = sum(length for op, length in alignment.cigartuples if op in (0, 2, 3))
                    positions.append(alignment_start+(alignment_length-1))
                else:
                    positions.append(alignment_start)
                           
        positions.sort()
        
        grouped_positions = group_positions(positions)
        results[chromosome_id] = grouped_positions

        num_groups = len(grouped_positions)
        num_groups_more_than_5 = 0

        for count in grouped_positions.values():
            if count > 5:
                num_groups_more_than_5 += 1
        print(f"{chromosome_id} \t {num_groups} \t {num_groups_more_than_5}")
    

    return results


def group_positions(positions):
    position_count = defaultdict(int)
    positions.sort()
    
    for position in positions:
        if position in position_count:
            position_count[position] += 1
        else:
            position_count[position] = 1
    
    # Add the last group
    #groups[(current_group_start, positions[-1])] += current_group_count
    
    return position_count

def save_results(results, ids_file):
    min_coverage = 10 
    directory = os.path.dirname(ids_file)
    output_file = os.path.join(directory, f"potentialBreakpoints_minCoverage{min_coverage}")

    with open(output_file, 'w') as f:
        for chromosome_id, grouped_positions in results.items():
            for position, count in grouped_positions.items():
                if count >= min_coverage:
                    f.write(f"{chromosome_id}\t{count}\t{position}\n")
def main():
    parser = argparse.ArgumentParser(description="Group start positions from SAM file based on chromosome IDs")
    parser.add_argument("ids_file", help="Path to file containing chromosome IDs")
    parser.add_argument("sam_file", help="Path to SAM file")
    args = parser.parse_args()

    results = process_chromosome_ids(args.ids_file, args.sam_file)
    
    save_results(results, args.ids_file)
    print("\nResults saved to 'potential_breakpoints' file.")

if __name__ == "__main__":
    main()
