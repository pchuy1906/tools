import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",  default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--cell_type", default='cell_3',   help='cell_3/cell_9/NON_ORTHO')

args        = parser.parse_args()
file_xyz    = args.file_xyz
cell_type   = args.cell_type

print ("")
print ("xyz_2_POSCAR")
print ("")
print ("read fileXYZ:", file_xyz)

natom, cell_3_3, atomList, xyz = utils.read_xyz(file_xyz, cell_type)

utils.write_POSCAR(natom, cell_3_3, atomList, xyz)
