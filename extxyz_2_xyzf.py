import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils
###############################################################################
from ase.io import read, write
###############################################################################
from ase.units import kcal, mol, Hartree, Bohr, GPa, eV, Angstrom, bar
# others
conv_Angstrom_2_Bohr = Angstrom/Bohr


# energy
kcal_per_mol = kcal/mol
conv_Hartree_2_kcal_per_mol = Hartree/kcal_per_mol
conv_eV_2_kcal_per_mol = eV/kcal_per_mol
conv_eV_2_Hartreertree = eV/Hartree

# forces
Hartree_per_Bohr = Hartree/Bohr
eV_per_Angstrom = eV/Angstrom

# stress
eV_per_Angstrom3 = eV/Angstrom**3
con_eV_per_Angstrom3_to_bar = eV_per_Angstrom3/bar
con_eV_per_Angstrom3_to_GPa = eV_per_Angstrom3/GPa
###############################################################################
import argparse
parser = argparse.ArgumentParser(description="Convert extxyz to xyzf for ChIMES optimization.")
parser.add_argument('--file_xyz',         default="aaa.xyz", help='Input XYZ file')
parser.add_argument('--utype',            default="ASE", help='unit in the extxyz file')
parser.add_argument('--uenergy',          default="", help='unit of energy in the extxyz file')
parser.add_argument('--uforces',          default="", help='unit of forces in the extxyz file')
parser.add_argument('--ustress',          default="", help='unit of stress in the extxyz file')

args = parser.parse_args()
file_xyz       = args.file_xyz
utype          = args.utype
uenergy        = args.uenergy
uforces        = args.uforces
ustress        = args.ustress
###############################################################################
# Read all frames from the extxyz file
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

frames = read(file_xyz, format="extxyz", index=":")

for atoms in frames:
    # Get energy, forces, and stress if present
    try:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        stress = atoms.get_stress(voigt=True)
    except Exception:
        energy = None
        forces = None
        stress = None

    # convert units
    #                 ChIMES      this unit
    #    Energy      kcal/mol       eV
    #    Forces      Ha/Bohr        eV/A
    #    Stress      GPa            eV/A3

    if utype == "ASE":
        energy *= eV / kcal_per_mol
        forces *= eV_per_Angstrom / Hartree_per_Bohr
        stress *= eV_per_Angstrom3 / GPa
    elif utype == "LAMMPS_metal":
        energy *= eV / kcal_per_mol
        forces *= eV_per_Angstrom / Hartree_per_Bohr
        stress *= bar / GPa

    cell_xyz = atoms.cell.array
    atomlist = list(atoms.symbols)
    xyz = atoms.get_positions()
    fname = "corrected_" + file_xyz

    utils._chimes_write_xyzf(fname, atomlist, xyz, cell_xyz, forces, energy, stress)
