import argparse
parser = argparse.ArgumentParser(description='NEB calculations')
# Arguments supported by the code.
parser.add_argument("--path_1", default='.', help='path state 1')
parser.add_argument("--path_2", default='.', help='path state 2')
parser.add_argument("--nimage", type=int, default=11, help='number of images')
args = parser.parse_args()
path_1 = args.path_1
path_2 = args.path_2
nimage = args.nimage

from ase.calculators.chimes_ase import ChIMES 
param_file = "../ChIMES_params_2026_July_13.txt"
calc = ChIMES(param_file)
calc2 = ChIMES(param_file)

from ase import io
from ase.calculators.lammpslib import LAMMPSlib
from ase.mep import NEB
from ase.optimize import QuasiNewton, FIRE
from ase.io import read, write
from ase.io.trajectory import Trajectory
import sys

initial = read(f'{path_1}/POSCAR')
final = read(f'{path_2}/POSCAR')
initial.calc = calc
final.calc = calc2
images = [initial]
for i in range(nimage-2):
    image=initial.copy()
    image.calc=ChIMES(param_file)
    images.append(image)
images.append(final)
neb = NEB(images, method='improvedtangent', climb=True)
neb.interpolate('idpp')
qn = FIRE(neb, trajectory='neb.traj')
qn.run(fmax=0.01)
