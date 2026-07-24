import argparse
parser = argparse.ArgumentParser(description='NEB tools')
# Arguments supported by the code.
parser.add_argument("--nimage", type=int, default=11, help='number of images')
args = parser.parse_args()
nimage = args.nimage

import matplotlib.pyplot as plt
from ase.mep.neb import NEBTools
from ase.io import read, write

images = read(f'neb.traj@-{nimage}:')

nebtools = NEBTools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()
print(Ef,dE)

# Get the barrier without any interpolation between highest images.
Ef, dE = nebtools.get_barrier(fit=False)
print(Ef,dE)

# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()
print(max_force)

with open("energies.dat", "w") as f:
    for i, image in enumerate(images):
        energy = image.get_potential_energy()
        f.write(f"{i:3d}  {energy:.12f}\n")

write('traj_neb.xyz', images, format='extxyz')
