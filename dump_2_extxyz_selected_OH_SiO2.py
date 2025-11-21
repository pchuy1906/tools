from ase.io import read, write
from ase.geometry import  get_distances
import numpy as np


import argparse
parser = argparse.ArgumentParser(description='export LAMMPS calculations to extxyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",           default='dump_xyz_vxyz',     help='file LAMMPS dump')
parser.add_argument("--nwater",    type=int, default=100,                 help='file LAMMPS dump')
parser.add_argument('--element_list', nargs='+')

args        = parser.parse_args()
file_dump        = args.file_dump
nwater           = args.nwater
element_list     = args.element_list

#file_dump = "dump_xyz_vxyz"
#nwater = 56
#element_list = ['H', 'O', 'Si']

atoms_list = read(file_dump, format='lammps-dump-text', index=':') 

atoms = atoms_list[0]
unique_atomic_numbers = sorted(set(atom.number for atom in atoms))

for atoms in atoms_list:
    for atom in atoms:
        atom.symbol = element_list[atom.number-1]

atoms = atoms_list[0]
cell = atoms.get_cell()

atoms_H = atoms[[atom.symbol == 'H' for atom in atoms]]
surface_atoms_H = atoms_H[:-2*nwater]
#print (surface_atoms_H)

atoms_O = atoms[[atom.symbol == 'O' for atom in atoms]]
surface_atoms_O = atoms_O[:-nwater]
#print (surface_atoms_O)

ID_Os = []
for i in range(len(surface_atoms_H)):
    HO_distances = get_distances(surface_atoms_H.positions[i], surface_atoms_O.positions, pbc=True, cell=cell)[1][0]
    min_index = np.argmin(HO_distances)
    if HO_distances[min_index] < 1.1:
        ID_Os.append(min_index)

all_frames = []
for atoms in atoms_list:
    atoms_H = atoms[[atom.symbol == 'H' for atom in atoms]]
    surface_atoms_H = atoms_H[:-2*nwater]
    atoms_O = atoms[[atom.symbol == 'O' for atom in atoms]]
    surface_atoms_O = atoms_O[ID_Os]
    surface_atoms = surface_atoms_H + surface_atoms_O
    all_frames.append(surface_atoms)


write('output.extxyz', all_frames, format='extxyz')

print("LAMMPS dump file successfully converted to extxyz format.")
