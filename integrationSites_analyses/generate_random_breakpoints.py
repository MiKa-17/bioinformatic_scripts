import random
import pandas as pd
import argparse


# Function to generate random points on the human genome
def generate_random_points(num_points, region_lengths):
    points = []
    for _ in range(num_points):
        random_chromosome = region_lengths.sample()
        length = int(random_chromosome['length'].iloc[0])
        position = random.randint(1, length)
        region = random_chromosome['region'].iloc[0]  # Access the region value
        points.append((region, position))
    return points


# Function to save random points to a file
def save_points_to_file(points, file_path):
    with open(file_path, 'w') as file:
        for point in points:
            file.write(f"{point[0]}\t1\t{point[1]}\n")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process BLAST results and annotation data")
    parser.add_argument("num_points", help="Path to file with genomic positions")
    parser.add_argument("region_lengths", help="Path to the annotation file")
    parser.add_argument("output_path", help="Path to save the CSV and plots")
    args = parser.parse_args()

    region_header = ["region", "length"]
    region_lengths = pd.read_csv(args.region_lengths, sep=',', header=None, names=region_header)
    region_lengths['length'] = region_lengths['length'].astype(int)


    num_points = int(args.num_points)
    output_file = args.output_path

    random_points = generate_random_points(num_points, region_lengths)
    save_points_to_file(random_points, output_file)

    print(f"{num_points} random points saved to {output_file}.")


"""
    'NC_000002.12': 242193529,
    'NC_000003.12': 198295559,
    'NC_000004.12': 190214555,
    'NC_000005.10': 181538259,
    'NC_000006.12': 170805979,
    'NC_000007.14': 159345973,
    'NC_000008.11': 145138636,
    'NC_000009.12': 138394717,
    'NC_000010.11': 133797422,
    'NC_000011.10': 135086622,
    'NC_000012.12': 133275309,
    'NC_000013.11': 114364328,
    'NC_000014.9': 107043718,
    'NC_000015.10': 101991189,
    'NC_000016.10': 90338345,
    'NC_000017.11': 83257441,
    'NC_000018.10': 80373285,
    'NC_000019.10': 58617616,
    'NC_000020.11': 64444167,
    'NC_000021.9': 46709983,
    'NC_000022.11': 50818468,
    'NC_000023.11': 156040895,
    'NC_000024.10': 57227415
"""