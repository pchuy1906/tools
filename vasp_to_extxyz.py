from ase.io import read, write
import numpy as np

input_file = "OUTCAR"
output_file = "output.extxyz"

# Read all frames from the OUTCAR file
atoms_list = read(input_file, index=":", format="vasp-out")

#num_atoms = len(atoms_list[0])
#q = np.zeros(num_atoms)
q = np.loadtxt("qBader.dat")

atoms_list[0].new_array('hirshfeld', q)

# Write all frames to a single extxyz file
write(output_file, atoms_list, format="extxyz")

print(f"Successfully converted {input_file} to {output_file}")
