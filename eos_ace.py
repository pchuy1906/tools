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
atoms = read_vasp(file='vcopt_POSCAR')
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
from ase.eos import EquationOfState

# List to store energies and volumes
volumes = []
energies = []

# Scale the cell to different volumes
scalings = np.linspace(0.96, 1.04, 7)  # 7 points around equilibrium

for scale in scalings:
    atoms_scaled = atoms.copy()
    atoms_scaled.set_cell(atoms.get_cell() * scale, scale_atoms=True)
    atoms_scaled.calc = PyACECalculator(file_parameter)
    volumes.append(atoms_scaled.get_volume())
    energies.append(atoms_scaled.get_potential_energy())

# Fit an equation of state
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()

print(f"Equilibrium volume: {v0:.2f} Å^3")
print(f"Equilibrium energy: {e0:.2f} eV")
print(f"Bulk modulus: {B:.2f} GPa")

# Optional: plot the EOS
eos.plot('eos.png')
