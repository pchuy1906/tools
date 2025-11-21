import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np
from pyace import PyACECalculator
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.io import write
###############################################################################
import argparse
parser = argparse.ArgumentParser(description='vcopt using ase')
# Arguments supported by the code.
parser.add_argument("--file_parameter",   default="aaa.txt",     help='parameter file')
args = parser.parse_args()
file_parameter = args.file_parameter
###############################################################################
# read the configuration

from ase.io.vasp import read_vasp
atoms = read_vasp(file='POSCAR')
natom = len(atoms)
print (natom, "atoms in the cell")

print ()
print ()

cell_xyz = atoms.cell[:]
cell_lengths, cell_angles = utils.calculate_cell_parameters(cell_xyz)
print ("initial cell parameter")
print (cell_xyz)
print (cell_lengths)
print (cell_angles)
print ("cell_volume", atoms.cell.volume)
################################################################################
## Now let's create a calculator from built in models:
atoms.calc = PyACECalculator(file_parameter)
################################################################################
## Now run the dynamics:
ucf = UnitCellFilter(atoms)
opt = BFGS(ucf, trajectory='traj.ase')
opt.run(fmax=0.01, steps=1000)
###############################################################################
cell_xyz = atoms.cell[:]
cell_lengths, cell_angles = utils.calculate_cell_parameters(cell_xyz)
print ()
print ()

print ("final cell parameter")
print (cell_xyz)
print (cell_lengths)
print (cell_angles)
print ("cell_volume", atoms.cell.volume)

write("vcopt_POSCAR", atoms, format="vasp")
