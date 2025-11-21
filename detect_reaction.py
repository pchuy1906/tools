import numpy as np

from ase import Atoms
from ase.io import read, Trajectory
from ase.geometry import  get_distances

## input parameters
import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                         default='file.xyz',         help='file_XYZ format xyz')
parser.add_argument("--natom_per_molecule", type=int,     default=10,                 help='number of atom per molecule')
parser.add_argument("--rcut",               type=float,   default=1.5,                help='cutoff threshold')
args        = parser.parse_args()
file_xyz           = args.file_xyz
natom_per_molecule = args.natom_per_molecule
rcut               = args.rcut


cell_xyz = np.zeros(shape=(3,3))
with open(file_xyz) as fileVar:
    cell3 = fileVar.readlines()[1].split()
for k in range(3):
    cell_xyz[k,k] = float(cell3[k])

trajectory = read(file_xyz, index=':')
iframe = 1

stop_running = False

for atoms in trajectory:
    atomList = list(atoms.symbols)
    xyz = atoms.get_positions()
    natom = len(atomList)
    n_molecule = natom//natom_per_molecule
    for i in range(n_molecule):
        i1 = natom_per_molecule*i+0
        i2 = natom_per_molecule*i+natom_per_molecule
        sub_atomList = atomList[i1:i2]
        sub_xyz = xyz[i1:i2]
        sub_atoms = Atoms(symbols=sub_atomList, positions=sub_xyz, cell=cell_xyz, pbc=True)
        tmp = sub_atoms.get_positions()
        for j in range(natom_per_molecule):
            tmp_arr = get_distances(sub_atoms.positions[j], sub_atoms.positions,pbc=True, cell=cell_xyz)
            tmp_dist = tmp_arr[1][0]
            tmp_dist[tmp_dist<0.01] = 100.0
            if np.min(tmp_dist) > rcut:
                print (iframe)
                stop_running = True
        if stop_running:
            break
    if stop_running:
        break

    iframe += 1



