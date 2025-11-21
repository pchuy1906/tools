import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='DFTB_2_xyzf')
# Arguments supported by the code.
parser.add_argument("--file_SPARC_static",                     default='SPARC.static', help='file SPARC static')
parser.add_argument("--file_SPARC_out",                        default='SPARC.out',    help='file SPARC out')
parser.add_argument("--cell_type",                             default='cell_3',       help='cell_3/cell_9/NON_ORTHO')
parser.add_argument("--export_stress",    action='store_true',                         help='export stress')


args = parser.parse_args()
file_SPARC_static    = args.file_SPARC_static
file_SPARC_out = args.file_SPARC_out
cell_type        = args.cell_type
export_stress    = args.export_stress


AtomList, Cxyz, energy, fxyz, stress = utils.read_SPARC_static(file_SPARC_static)
cell_xyz = utils.read_SPARC_out(file_SPARC_out)
xyz = Cxyz @ cell_xyz
natom = len(AtomList)

print (stress)
fname = "input.xyzf"
export_stress = True
utils.write_xyzf(fname, natom, AtomList, xyz, cell_xyz, cell_type, fxyz, energy, stress, export_stress)

