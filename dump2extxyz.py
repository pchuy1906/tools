from ase.units import kcal, mol, Hartree, Bohr, GPa, eV, Angstrom, bar
kcal_per_mol = kcal/mol

#conv_Hartree_2_kcal_per_mol = Hartree/kcal_per_mol
#conv_eV_2_kcal_per_mol = eV/kcal_per_mol
#conv_eV_2_Hartreertree = eV/Hartree
#conv_Angstrom_2_Bohr = Angstrom/Bohr
#eV_per_Angstrom3 = eV/Angstrom**3
#con_eV_per_Angstrom3_to_bar = eV_per_Angstrom3/bar
#con_eV_per_Angstrom3_to_GPa = eV_per_Angstrom3/GPa

conv_kcal_per_mol_2_eV = kcal_per_mol / eV
conv_atm_2_bar = 1.01325


import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

from ase.io import read
from ase import Atoms

import argparse
parser = argparse.ArgumentParser(description='export LAMMPS calculations to extxyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",           default='dump_xyz_vxyz',  help='file LAMMPS dump')
parser.add_argument("--file_elements",       default='ELEMENTS',       help='file ELEMENTS')
parser.add_argument("--units_input_LAMMPS",  default='metal',          help='unit used in LAMMPS')
parser.add_argument("--units_output_extxyz", default='metal',          help='unit used in extxyz')

args        = parser.parse_args()
file_dump        = args.file_dump
file_elements    = args.file_elements
units_input_LAMMPS  = args.units_input_LAMMPS
units_output_extxyz = args.units_output_extxyz

atoms = read(file_dump, format='lammps-dump-text', index=-1)


print(f"OLD symbols: {atoms.symbols}")
print(atoms.get_chemical_symbols())
print(atoms.numbers)

with open(file_elements, "r") as f:
    element_list = sorted(f.read().split())
old_atom_numbers = atoms.numbers
new_atom_symbols = [element_list[i-1] for i in old_atom_numbers]
atoms.symbols = new_atom_symbols

print(f"NEW symbols: {atoms.symbols}")
print(atoms.get_chemical_symbols())
print(atoms.numbers)

energy, stress = utils.read_energy_stress(file_LAMMPS_log="log.lammps")
forces = atoms.calc.results['forces']

if units_input_LAMMPS == "real":
    energy *= conv_kcal_per_mol_2_eV
    forces *= conv_kcal_per_mol_2_eV
    stress *= conv_atm_2_bar

from ase.calculators.singlepoint import SinglePointCalculator

selected_atoms = Atoms(
    symbols=atoms.get_chemical_symbols(),
    positions=atoms.get_positions(),
    cell=atoms.get_cell(),
    pbc=atoms.get_pbc(),
)

print (energy)
# Attach the calculator with forces
calc = SinglePointCalculator(
    selected_atoms,
    energy=energy,
    stress=stress,
    forces=forces,
)


selected_atoms.calc = calc
from ase.io import write
write('output.xyz', selected_atoms, format='extxyz')
