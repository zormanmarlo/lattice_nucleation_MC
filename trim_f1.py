import numpy as np
from sys import argv

# Path for the new output file
output_path = "f1"
size = int(argv[2])
file_path = str(argv[1])

# Read the input file, filter the lines, and write the valid lines to the output file
with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
    for line in infile:
        # Split the line into values based on whitespace
        parts = line.split()
        # Convert the first three columns to integers
        first_three_columns = [int(parts[i]) for i in range(3)]
        # Check if all values are 100 or less
        if all(x <= size for x in first_three_columns):
            # Write the valid line to the output file
            outfile.write(line)

