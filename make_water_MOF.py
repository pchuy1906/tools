import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import  get_distances
from ase.constraints import FixAtoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS

import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils


# read input 
import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument("--file_xyz_1",        default='xyz', help='file xyz')
parser.add_argument("--file_xyz_2",        default='xyz', help='file xyz')
parser.add_argument("--rcut",  type=float, default=1.5,   help='rcut')
parser.add_argument("--nwmax", type=int,   default=50,    help='nwmax')

args    = parser.parse_args()
file_xyz_1  = args.file_xyz_1
file_xyz_2  = args.file_xyz_2
rcut        = args.rcut
nwmax       = args.nwmax

cell_type = "cell_9"
natom, cell_3_3, atomList, xyz = utils.read_xyz(file_xyz_1, cell_type)
atoms1 = Atoms(symbols=atomList, positions=xyz, cell=cell_3_3, pbc=True)


cell_type = "cell_9"
natom2, cell_3_3_2, atomList2, xyz2 = utils.read_xyz(file_xyz_2, cell_type)

nwat = natom2//3
collect_xyz = []
collect_sym = []
all_atoms = atoms1

for i in range(nwat):

    dist_arr = np.array([])

    tmp_xyz = []

    i1 = 3*i
    TMP = get_distances(all_atoms.positions, xyz2[i1],pbc=True, cell=cell_3_3)
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(xyz2[i1].tolist())

    i1 = 3*i+1
    TMP = get_distances(all_atoms.positions, xyz2[i1],pbc=True, cell=cell_3_3)
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(xyz2[i1].tolist())

    i1 = 3*i+2
    TMP = get_distances(all_atoms.positions, xyz2[i1],pbc=True, cell=cell_3_3)
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(xyz2[i1].tolist())

    tmp_sym = ['H','H','O']

    if nwmax<1:
        break
    if (np.min(dist_arr) > rcut):
        collect_xyz.extend(tmp_xyz)
        collect_sym.extend(tmp_sym)
        selected_atoms = Atoms(symbols=collect_sym, positions=collect_xyz, cell=cell_3_3, pbc=True)
        all_atoms = atoms1 + selected_atoms
        if len(collect_sym)==3*nwmax:
            break

print ("number of adding water", len(collect_sym)//3)


cell_str = " ".join(f"{x:.6f}" for x in cell_3_3.flatten())

output_path = "xyz_all.xyz"
write(output_path, all_atoms, format='xyz', comment=cell_str)

output_path = "xyz_water.xyz"
write(output_path, selected_atoms, format='xyz', comment=cell_str)
