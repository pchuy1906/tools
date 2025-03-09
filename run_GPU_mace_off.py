from mace.calculators import mace_anicc
from mace.calculators import mace_off
from mace.calculators import MACECalculator

from ase import units
import numpy as np

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

# initial velocities
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
MaxwellBoltzmannDistribution(atoms, temperature_K=T)
#atoms.set_velocities(np.full((natom, 3), 1e-1))

###############################################################################
# Now let's create a calculator from built in models:
#calculator = torchani.models.ANI1xnr().ase()
calculator = mace_off(device='cuda')
#calculator = MACECalculator(model_paths="../../ani500k_large_CC.model", device='cuda')
atoms.set_calculator(calculator)

###############################################################################
# Now create a callback function that print interesting physical quantities:
def printenergy(a=atoms):
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    etot = epot + ekin
    temp = ekin / (1.5 * units.kB)
    stress = a.get_stress(include_ideal_gas=True)
    stress_GPa = stress / units.bar * 0.0001
    press = (stress_GPa[0] + stress_GPa[1] + stress_GPa[2])/3.0
    #print ("stress=", stress)
    #print ("stress_GPa=",stress_GPa)
    print('Energy per atom: T = %.3fK P = %.3fGPa Etot = %.3feV  Epot = %.3feV Ekin = %.3feV' 
       % (temp, press, etot, epot, ekin))


###############################################################################
from ase.md.langevin import Langevin
dyn = Langevin(atoms, timestep= dt * units.fs, temperature_K=T, friction=0.02)
dyn.attach(printenergy, interval=iprint_ener)
from ase.io.trajectory import Trajectory, TrajectoryReader
traj_file = "traj.ase"
traj = Trajectory(traj_file, 'w', atoms)
dyn.attach(traj.write, interval=iprint_traj)

###############################################################################
# Now run the dynamics:
print("Beginning dynamics...")
#printenergy()
dyn.run(nstep)
