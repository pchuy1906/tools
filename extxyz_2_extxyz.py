
from ase.io import read, write
from ase.cell import Cell
from ase import Atoms
from ase.geometry import cell_to_cellpar


eVA3_to_bar = 1602176.6

###############################################################################
import argparse
parser = argparse.ArgumentParser(description="Convert extxyz to pckl.gzip for PACE optimization.")
parser.add_argument('--file_xyz',         default="aaa.xyz", help='Input XYZ file')
parser.add_argument("--good_cell",        default=False, action="store_true", help="rotate cell or not")
args = parser.parse_args()
file_xyz       = args.file_xyz
good_cell      = args.good_cell
###############################################################################
# Read all frames from the extxyz file
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

frames = read(file_xyz, format="extxyz", index=":")

for atoms in frames:

    oldcell9 = atoms.cell.array
    oldcell6 = cell_to_cellpar(atoms.cell)
    newcell9 = Cell.fromcellpar(oldcell6)

    # for convenient, we want :             newxyz[natomx3] = oldxyz[natomx3] @ R(3x3)
    # so the selected transformation is:    R = [oldcell]^{-1} @ newcell
    # or                                    newcell = oldcell @ R

    if good_cell:
        R = np.linalg.inv(oldcell9) @ newcell9
    else:
        R = np.eye(3)

    atomlist = list(atoms.symbols)
    oldxyz = atoms.get_positions()

    # change the cell_parameter, xyz, and stress
    newcell9 = oldcell9 @ R
    newxyz = oldxyz @ R

    old_stress = atoms.get_stress(voigt=False)
    R_T = R.T
    new_stress = R_T @ old_stress @ R  * eVA3_to_bar

    # change the forces
    try:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    except Exception:
        energy = None
        forces = None
    new_forces = forces @ R

    atoms.set_cell(newcell9)
    atoms.set_positions(newxyz)

    # Attach new calculator with updated stress
    atoms.calc = SinglePointCalculator(
        atoms,
        energy=energy,
        forces=new_forces,
        stress=new_stress
    )

output_file = "corrected_" + file_xyz
write(output_file, frames, format='extxyz')
