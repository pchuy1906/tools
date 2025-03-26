from ase.io import read, write

input_file = "OUTCAR"
output_file = "output.extxyz"

# Read all frames from the OUTCAR file
atoms_list = read(input_file, index=":", format="vasp-out")

# Write all frames to a single extxyz file
write(output_file, atoms_list, format="extxyz")

print(f"Successfully converted {input_file} to {output_file}")
