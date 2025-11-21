import numpy as np

import argparse
parser = argparse.ArgumentParser(description='calculate average')

# Arguments supported by the code.
parser.add_argument("--file_data",         default='data.dat', help='file data')

args      = parser.parse_args()
file_data = args.file_data

F = np.loadtxt(file_data)
X0 = F[:,0]
Y0 = F[:,1]
X0 = np.array(X0)
Y0 = np.array(Y0)

integral = sum(Y0)*(X0[1]-X0[0])

F[:,1] = F[:,1]/integral

# Specify the format: 8 characters wide, 2 decimal places
fmt = "{:8.2f}"

with open('normalized_' + file_data, 'w') as f:
    for row in F:
        # Format each number in the row
        formatted_row = ' '.join(fmt.format(num) for num in row)
        f.write(formatted_row + '\n')
