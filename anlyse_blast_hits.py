import argparse

def get_matching_hits(input_file, endPositonOfInsert):
    matching_hits = []
    
    with open(input_file, "r") as blast_file:
        for line in blast_file:
            columns = line.strip().split('\t')
            
            if len(columns) >= 8 and columns[6] == '1' and int(columns[7]) == endPositonOfInsert:
                matching_hits.append(line)
    
    return matching_hits

def count_reads_with_hits(matching_hits):
    read_hits = {}
    
    for line in matching_hits:
        columns = line.strip().split('\t')
        read_id = columns[1]
        
        if read_id in read_hits:
            read_hits[read_id] += 1
        else:
            read_hits[read_id] = 1
    
    return read_hits

def group_reads_by_hit_count(read_hits):
    read_counts = {}
    
    for read_id, hit_count in read_hits.items():
        if hit_count in read_counts:
            read_counts[hit_count].append(read_id)
        else:
            read_counts[hit_count] = [read_id]
    
    return read_counts

def count_unique_read_ids(input_file):
    unique_read_ids = set()
    
    with open(input_file, "r") as blast_file:
        for line in blast_file:
            columns = line.strip().split('\t')
            
            if len(columns) >= 2:
                unique_read_ids.add(columns[1])
    
    return len(unique_read_ids)

def save_read_lists(read_counts, output_path):
    for hits, read_list in read_counts.items():
        with open(f"{output_path}/reads_with_{hits}_hits.txt", "w") as output_file:
            for read_id in read_list:
                output_file.write(f"{read_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Process BLAST results and group reads by hit count.")
    parser.add_argument("input_file", help="Path to the BLAST table file")
    parser.add_argument("output_path", help="Path to the output directory")
    parser.add_argument("--endPos", type=int, required=True, help="Specific number to match in column 8")
    args = parser.parse_args()

    unique_read_count = count_unique_read_ids(args.input_file)

    with open(args.input_file, "r") as blast_file:
        total_lines = sum(1 for _ in blast_file)  # Count total lines in the file

    matching_hits = get_matching_hits(args.input_file, args.endPos)
    
    read_hits = count_reads_with_hits(matching_hits)
    read_counts = group_reads_by_hit_count(read_hits)

    # Calculate total number of reads with hits
    total_reads_with_hits = sum(len(read_list) for read_list in read_counts.values())
    
    print(f"Number of blast hits including the whole insert including TTAA: {len(matching_hits)} of {total_lines} hits")
    print(f"Number of reads including the whole insert including TTAA: {total_reads_with_hits} of {unique_read_count} reads")

    print("Reads with Multiple Hits:")
    for hits, read_list in read_counts.items():
        print(f"{len(read_list)} reads have {hits} hits")
    
    # Save the read_lists to separate files
    save_read_lists(read_counts, args.output_path)

if __name__ == "__main__":
    main()
