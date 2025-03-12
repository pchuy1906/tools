import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import  get_distances
from ase.constraints import FixAtoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
import xtb
from xtb.ase.calculator import XTB
#from ase.calculators.xtb import XTB

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR_1", default='POSCAR', help='file with format xyzf, xyzfe, xyzfes')
parser.add_argument("--file_POSCAR_2", default='POSCAR', help='file with format xyzf, xyzfe, xyzfes')

args    = parser.parse_args()
poscar1_path  = args.file_POSCAR_1
poscar2_path  = args.file_POSCAR_2

atoms1 = read(poscar1_path)
atoms2 = read(poscar2_path)
merged_atoms = atoms1 + atoms2

natom1 = len(atoms1)
natom2 = len(atoms2)

nwat = natom2//3

collect_xyz = []
collect_sym = []

cell = atoms1.get_cell()

for i in range(nwat):

    dist_arr = np.array([])

    tmp_xyz = []

    i1 = natom1 + 2*i
    TMP = get_distances(merged_atoms.positions[0:natom1], merged_atoms.positions[i1])
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(merged_atoms.positions[i1].tolist())

    i1 = natom1 + 2*i+1
    TMP = get_distances(merged_atoms.positions[0:natom1], merged_atoms.positions[i1])
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(merged_atoms.positions[i1].tolist())


    i1 = natom1 + 2*nwat + i
    TMP = get_distances(merged_atoms.positions[0:natom1], merged_atoms.positions[i1])
    dist_arr = np.append(dist_arr, np.array(TMP[1]))
    tmp_xyz.append(merged_atoms.positions[i1].tolist())

    tmp_sym = ['H','H','O']

    if (np.min(dist_arr) > 1.5):
        collect_xyz.extend(tmp_xyz)
        collect_sym.extend(tmp_sym)

print (collect_sym)
print (collect_xyz)

selected_atoms = Atoms(symbols=collect_sym, positions=collect_xyz, cell=cell, pbc=True)

all_atoms = atoms1 + selected_atoms
#output_path = "POSCAR_all"
#write(output_path, all_atoms, format='vasp')

# fix atoms
fixed_indices = [i for i in range(natom1)]
constraint = FixAtoms(indices=fixed_indices) 
all_atoms.set_constraint(constraint)

#all_atoms.calc = LennardJones()
all_atoms.calc = XTB(method='GFN2-xTB')

optimizer = BFGS(all_atoms)  
optimizer.run(fmax=2.961)  

output_path = "POSCAR_all"
write(output_path, all_atoms, format='vasp')

