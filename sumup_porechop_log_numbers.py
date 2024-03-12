import os

# Initialize variables to store sums
start_sum = 0
end_sum = 0
split_sum = 0
total_reads = 0

dir_path = '/media/michele/MKY-Data/VVGT/lenti_insertion/fastq_pass/barcode11/porechop_logs/'

# Iterate over each file
for filename in os.listdir(dir_path):
    if filename.endswith('.log'):
        with open(os.path.join(dir_path, filename), 'r') as file:
            # Use list comprehension to filter lines starting with numbers followed by "adapters"
            for line in file:
                line = line.strip()
                if line and line[0].isdigit() and "adapters" in line:
                    parts = line.split()
                    count = int(parts[0].replace(',', ''))
                    

                    if "start" in line:
                        start_sum += count
                    elif "end" in line:
                        end_sum += count
                    elif "split" in line:
                        split_sum += count
            total_reads += int(parts[2].replace(',', ''))
            
# Output the sums
print(f"{start_sum} / {total_reads} reads had adapters trimmed from their start")
print(f"{end_sum} / {total_reads} reads had adapters trimmed from their end")
print(f"{split_sum} / {total_reads} reads were split based on middle adapters")