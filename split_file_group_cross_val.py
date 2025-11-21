import random
import os

import argparse


def split_file_randomly(input_file, num_parts):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return

    random.shuffle(lines)
    chunk_size = len(lines) // num_parts
    remainder = len(lines) % num_parts

    for i in range(num_parts):
        part_name = f"{os.path.splitext(input_file)[0]}_part{i+1}.txt"
        with open(part_name, 'w') as outfile:
            start = i * chunk_size + min(i, remainder)
            end = (i + 1) * chunk_size + min(i + 1, remainder)
            outfile.writelines(lines[start:end])
        print(f"Created: {part_name}")

        output_file = "file_split_" + str(i+1) + ".xyz"
        with open(output_file, 'w') as outfile:
            with open(part_name, "r") as file:
                for line in file:
                    filename= line.strip()
                    with open(filename, 'r') as infile:
                        outfile.write(infile.read())
                        outfile.write('\n')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
    # Arguments supported by the code.
    parser.add_argument("--file_list",             default='list_files.txt',  help='list of xyz files')
    parser.add_argument("--nfold",      type=int,  default=1,                 help='4, number of subfiles')

    args        = parser.parse_args()

    file_list = args.file_list
    nfold     = args.nfold
    split_file_randomly(file_list, nfold)
