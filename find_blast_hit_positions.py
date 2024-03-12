import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Process BLAST results and group reads by hit count.")
parser.add_argument("input_file", help="Path to the BLAST table file")
parser.add_argument("--startBackbone", type=int, required=True, help="Specific number to match in column 8")
args = parser.parse_args()



# Read the tab-separated table into a Pandas DataFrame
# Replace 'your_file.csv' with the actual filename
df = pd.read_csv(args.input_file, sep='\t', header=None)


# Select only the necessary columns
df = df[[1, 6, 7]]  # Columns 1, 6, and 7 correspond to ReadsID, StartPosition, and StopPosition

# Rename columns for clarity
df.columns = ['ReadsID', 'StartPosition', 'StopPosition']

# Group the reads by alignment position and count them
read_groups = df.groupby(['StartPosition', 'StopPosition'])['ReadsID'].count().reset_index()

# Iterate through the grouped data and print the results
for index, row in read_groups.iterrows():
    start = row['StartPosition']
    stop = row['StopPosition']
    count = row['ReadsID']
    #print(f"A group of {count} reads have a blast hit from position {start} to {stop}")


# Create a second group for reads with stop position > 960
second_group = read_groups[read_groups['StopPosition'] > args.startBackbone]

# Make ReadsID unique within the second group using .loc
second_group.loc[:, 'ReadsID'] = second_group.groupby(['StartPosition', 'StopPosition'])['ReadsID'].transform('nunique')

# Count the reads in the second group
count_in_second_group = second_group['ReadsID'].sum()

# Print the result for the second group
print(f"{count_in_second_group} reads include the backbone")



# Filter read_groups based on specific start and stop positions
start_position_min = 955
start_position_max = 970
stop_position_min = 1718
stop_position_max = 1728

filtered_groups = read_groups[
    (read_groups['StartPosition'] >= start_position_min) &
    (read_groups['StartPosition'] <= start_position_max) &
    (read_groups['StopPosition'] >= stop_position_min) &
    (read_groups['StopPosition'] <= stop_position_max)
]

# Make ReadsID unique within the filtered groups and count them
filtered_groups['ReadsID'] = filtered_groups.groupby(['StartPosition', 'StopPosition'])['ReadsID'].transform('nunique')

# Get the total count of reads within the filtered groups
total_count = filtered_groups['ReadsID'].sum()

print(f"{total_count} reads have unique IDs within the filtered groups.")