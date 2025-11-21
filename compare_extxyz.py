import sys
import numpy as np
from ase.io import read, write

import matplotlib.pyplot as plt


###############################################################################
import argparse
parser = argparse.ArgumentParser(description="Convert extxyz to pckl.gzip for PACE optimization.")
parser.add_argument('--file_xyz_1',         default="aaa.xyz", help='Input XYZ file')
parser.add_argument('--file_xyz_2',         default="aaa.xyz", help='Input XYZ file')
args = parser.parse_args()
file_xyz_1       = args.file_xyz_1
file_xyz_2       = args.file_xyz_2
###############################################################################
def read_extxyz(file_xyz):
    configs = read(file_xyz, index=':')
    ncount = 0
    ref_forces = []
    ref_energy = []
    ref_stress = []
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
    ref_forces = np.vstack(ref_forces)
    ref_forces_flat = np.concatenate(ref_forces).flatten()
    ref_energy = np.vstack(ref_energy)
    ref_energy_flat = np.concatenate(ref_energy).flatten()
    ref_stress = np.vstack(ref_stress)
    ref_stress_flat = np.concatenate(ref_stress).flatten()
    return ref_forces_flat, ref_energy_flat, ref_stress_flat
###############################################################################
def write_comparison(ref_forces, ace_forces, fileoutput, fconversion):
    ref_forces = np.vstack(ref_forces)
    ace_forces = np.vstack(ace_forces)
    ref_forces_flat = np.concatenate(ref_forces).flatten()*fconversion
    ace_forces_flat = np.concatenate(ace_forces).flatten()*fconversion
    data = np.column_stack((ref_forces_flat, ace_forces_flat))
    np.savetxt(fileoutput, data, fmt='%.6f')

f1, e1, s1 = read_extxyz(file_xyz_1)
f2, e2, s2 = read_extxyz(file_xyz_2)

bar2GPa = 0.0001
write_comparison(f1, f2, fileoutput="compare_forces.dat",fconversion=1)
write_comparison(e1, e2, fileoutput="compare_energy.dat",fconversion=1)
write_comparison(s1, s2, fileoutput="compare_stress.dat",fconversion=bar2GPa)

