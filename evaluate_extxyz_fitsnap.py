import sys
import numpy as np
#import pandas as pd
#from ase import Atoms
from ase.io import read, write
#from ase.io.extxyz import read_extxyz
#import re, sys
#from typing import Dict, List, Any, Optional, Union
import gzip
import pickle
#import warnings

import matplotlib.pyplot as plt

from ase.calculators.lammpsrun import LAMMPS


###############################################################################
import argparse
parser = argparse.ArgumentParser(description="Convert extxyz to pckl.gzip for PACE optimization.")
parser.add_argument('--file_xyz',         default="aaa.xyz", help='Input XYZ file')
parser.add_argument("--file_parameter",   default="aaa.txt", help='parameter file')
args = parser.parse_args()
file_xyz       = args.file_xyz
file_parameter = args.file_parameter
###############################################################################

configs = read(file_xyz, index=':')

ncount = 0
ref_forces = []
ace_forces = []
ref_energy = []
ace_energy = []
ref_stress = []
ace_stress = []
for atoms in configs:

    ncount += 1

    TMP_ref_forces = np.array(atoms.get_forces(), dtype = float)
    ref_forces.append(TMP_ref_forces)
    TMP_ref_energy = atoms.get_potential_energy()
    ref_energy.append(TMP_ref_energy)
    try:
        TMP_ref_stress = atoms.get_stress(voigt=False)
        tstress = True
        ref_stress.append(TMP_ref_stress)
    except NotImplementedError:
        tstress = False

    atoms.calc = LAMMPS(
        pair_style='hybrid/overlay zbl 4.0 4.8 snap',
        pair_coeff=[
            '* * zbl 94',
            '* * snap fitsnap_pot.snapcoeff fitsnap_pot.snapparam Pu'
        ],
        files=['fitsnap_pot.snapcoeff', 'fitsnap_pot.snapparam'],
        keep_tmp_files=True,
        command='/usr/workspace/pham20/codes/FITSNAP/lammps/build-fitsnap/lmp'
    )

    TMP_ace_forces = np.array(atoms.get_forces(), dtype = float)
    ace_forces.append(TMP_ace_forces)
    TMP_ace_energy = atoms.get_potential_energy()
    ace_energy.append(TMP_ace_energy)
    if tstress:
        TMP_ace_stress = -1.0*atoms.get_stress(voigt=False)
        ace_stress.append(TMP_ace_stress)

def write_comparison(ref_forces, ace_forces, fileoutput, fconversion):
    ref_forces = np.vstack(ref_forces)
    ace_forces = np.vstack(ace_forces)
    ref_forces_flat = np.concatenate(ref_forces).flatten()*fconversion
    ace_forces_flat = np.concatenate(ace_forces).flatten()*fconversion
    data = np.column_stack((ref_forces_flat, ace_forces_flat))
    np.savetxt(fileoutput, data, fmt='%.6f')

eV_A3_to_GPa = 160.2

write_comparison(ref_forces, ace_forces, fileoutput="compare_forces.dat",fconversion=1)
write_comparison(ref_energy, ace_energy, fileoutput="compare_energy.dat",fconversion=1)
write_comparison(ref_stress, ace_stress, fileoutput="compare_stress.dat",fconversion=eV_A3_to_GPa)
