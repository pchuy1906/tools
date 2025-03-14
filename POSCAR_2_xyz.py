import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='POSCAR 2 xyz')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",    default="POSCAR", help=' POSCAR/CONTCAR')
parser.add_argument("--cell_type",      default="cell_3", help=' cell_3/cell_9/NON_ORTHO')

args = parser.parse_args()
file_POSCAR  = args.file_POSCAR
cell_type    = args.cell_type

natom, cell, atomList, xyz = utils.read_POSCAR(file_POSCAR)

fname = "POSCAR2xyz.xyz"
utils.write_xyz(fname, natom, cell_type, cell, atomList, xyz)
