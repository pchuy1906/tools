import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='gen_2_xyz')
# Arguments supported by the code.
parser.add_argument("--file_dftb_gen",                         default='geo_end.gen',  help='file DFTB gen')
parser.add_argument("--cell_type",                             default='cell_3',       help='cell_3/cell_9/NON_ORTHO')

args = parser.parse_args()
file_dftb_gen    = args.file_dftb_gen
cell_type        = args.cell_type

natom, cell_xyz, AtomList, xyz = utils.read_gen(file_dftb_gen)

fname = "gen_2_xyz.xyz"
utils.write_xyz(fname, natom, cell_type, cell_xyz, AtomList, xyz)
