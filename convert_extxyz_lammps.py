import os
import sys
import numpy as np
from ase.io import read, write
import subprocess
###############################################################################
import argparse
parser = argparse.ArgumentParser(description="Convert extxyz to pckl.gzip for PACE optimization.")
parser.add_argument('--file_xyz',         default="aaa.xyz", help='Input XYZ file')
args = parser.parse_args()
file_xyz       = args.file_xyz
###############################################################################
configs = read(file_xyz, index=':')

fname = "VALIDATION"
if os.path.isdir(fname):
    print(f"The folder '{fname}' exists. Check if you want to delete it")
    sys.exit(1)
else:
    subprocess.run(['mkdir', fname])

ncount = 0
for atoms in configs:
    ncount += 1
    unique_elements = sorted(set(atoms.get_chemical_symbols()))
    print(unique_elements)
    open('ELEMENTS', 'w').write(' '.join(sorted(set(atoms.get_chemical_symbols()))) + '\n')
    file_lammps_data = "data.lammps"
    write(file_lammps_data, atoms, format='lammps-data', atom_style='atomic', masses=True)
    sfold = fname+"/run-"+str(ncount)
    subprocess.run(['mkdir', sfold])
    subprocess.run(['mv', file_lammps_data, sfold])
    subprocess.run(['mv', 'ELEMENTS', sfold])



