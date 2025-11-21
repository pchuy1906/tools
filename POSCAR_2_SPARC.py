import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='POSCAR 2 SPARC')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",              default="POSCAR",                                             help=' POSCAR/CONTCAR')
parser.add_argument("--pathPP",                   default="/usr/workspace/pham20/codes/SPARC/SPARC-2.0.0/psps", help=' path to SPARC pseudopotential')
parser.add_argument("--MESH_SPACING", type=float, default=0.1,                                                  help=' MESH_SPACING')
parser.add_argument("--ECUT",         type=float, default=-100,                                                 help=' ECUT')
parser.add_argument("--DFT_METHOD",               default="GGA_PBE",                                            help=' DFT_METHOD')


args = parser.parse_args()
file_POSCAR  = args.file_POSCAR
pathPP       = args.pathPP
MESH_SPACING = args.MESH_SPACING
ECUT         = args.ECUT
DFT_METHOD   = args.DFT_METHOD

if (ECUT>0.0):
    MESH_SPACING = []
else:
    ECUT = []

natom, cell, atomList, xyz = utils.read_POSCAR_SPARC(file_POSCAR)

utils.write_ion_SPARC(natom, atomList, xyz, pathPP)
utils.write_input_SPARC(cell, MESH_SPACING, ECUT, DFT_METHOD)

