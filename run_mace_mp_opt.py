import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils
from utils import conv_eV_2_kcalmol, conv_eV_2_Ha, conv_Angstrom_2_Bohr


from mace.calculators import mace_mp

from ase import units
import numpy as np


from ase.optimize import BFGS

import argparse
parser = argparse.ArgumentParser(description='MD simulation using torchANI')
# Arguments supported by the code.
parser.add_argument('--T',           type=float, default=300.0, help="temperature (K)")
parser.add_argument('--dt',          type=float, default=0.1,   help="time step (fs) in a MD simulation")
parser.add_argument("--nstep",       type=int,   default=10,    help='number of steps in a MD simulation')
parser.add_argument("--iprint_ener", type=int,   default=5,     help='frequently print out energy')
parser.add_argument("--iprint_traj", type=int,   default=5,     help='frequently write traj file')
args = parser.parse_args()
T           = args.T
dt          = args.dt
nstep       = args.nstep
iprint_ener = args.iprint_ener
iprint_traj = args.iprint_traj

###############################################################################
# read the configuration

#from ase.lattice.cubic import Diamond
#atoms = Diamond(symbol="C", pbc=True)

from ase.io.vasp import read_vasp
atoms = read_vasp(file='POSCAR')
print (atoms)

natom = len(atoms)
print (natom, "atoms in the cell")

###############################################################################
# Now let's create a calculator from built in models:
#calculator = torchani.models.ANI1xnr().ase()
calculator = mace_mp() # return the default medium ASE calculator equivalent to mace_mp(model="medium") in MACE < 0.3.10 and mace_mp(model="medium-mpa-0") in MACE >= 0.3.10
#calculator = mace_mp(model="large") # return a larger model
#calculator = mace_mp(model="https://tinyurl.com/y7uhwpje") # downlaod the model at the given url
#calculator = mace_mp(dispersion=True) # return a model with D3 dispersion correction

atoms.calc = calculator
###############################################################################
opt = BFGS(atoms, trajectory='traj.ase')
opt.run(fmax=0.5)
###############################################################################
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = -1.0*atoms.get_stress(voigt=False)
# convert to ChIMES input
energy = energy * conv_eV_2_kcalmol
forces = forces * conv_eV_2_Ha/conv_Angstrom_2_Bohr
stress = stress / units.bar * 0.0001
print (energy)
print (forces)
print (stress)

fname = "input.xyzf"
cell_type = "NON_ORTHO"
export_stress = True
atomList = list(atoms.symbols)
cell_xyz = atoms.cell[:]
xyz = atoms.get_positions()

utils.write_xyzf(fname, natom, atomList, xyz, cell_xyz, cell_type, forces, energy, stress, export_stress)


