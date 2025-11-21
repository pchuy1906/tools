import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
from ase.io import read, write
from ase.geometry import  get_distances

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument('--atom_types', nargs='+')
parser.add_argument('--fscale', type=float, default=1.0)

args    = parser.parse_args()
atom_types = args.atom_types
fscale     = args.fscale

for i in range(len(atom_types)):
    atomic_number_1 = atomic_numbers[atom_types[i]]
    atomic_number_2 = atomic_numbers[atom_types[i]]
    atomic_radius_1 = covalent_radii[atomic_number_1]
    atomic_radius_2 = covalent_radii[atomic_number_2]
    bond_length = atomic_radius_1 + atomic_radius_2
    bond_length *= fscale
    print (f"{atom_types[i]} {atom_types[i]} {bond_length}")

for i in range(len(atom_types)):
    for j in range(i+1,len(atom_types)):
        atomic_number_1 = atomic_numbers[atom_types[i]]
        atomic_number_2 = atomic_numbers[atom_types[j]]
        atomic_radius_1 = covalent_radii[atomic_number_1]
        atomic_radius_2 = covalent_radii[atomic_number_2]
        bond_length = atomic_radius_1 + atomic_radius_2
        bond_length *= fscale
        print (f"{atom_types[i]} {atom_types[j]} {bond_length}")

