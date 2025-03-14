import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np


import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",           default="POSCAR", help=' POSCAR/CONTCAR')

args = parser.parse_args()
file_POSCAR  = args.file_POSCAR

natom, unitcell, atomList, xyz = utils.read_POSCAR(file_POSCAR)

fname = "POSCAR2gen.gen"
utils.write_gen(fname, natom, unitcell, atomList, xyz)
