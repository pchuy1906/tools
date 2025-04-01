import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np


import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",              default="POSCAR",             help=' POSCAR')
parser.add_argument("--file_out_multicharge",     default="OUTPUT_MULTICHARGE", help=' output_multicharge')


args = parser.parse_args()
file_POSCAR          = args.file_POSCAR
file_out_multicharge = args.file_out_multicharge

natom, cell_xyz, atomList, xyz = utils.read_POSCAR(file_POSCAR)

q, energy, fxyz, stress = utils.read_multicharge(file_out_multicharge, natom)

with open('EEQ_q.dat', 'w') as f:
    np.savetxt(f, q, fmt='%12.4f')

print (stress)
fname = "EEQ.xyzf"
cell_type = "NON_ORTHO"
export_stress = True
utils.write_xyzf(fname, natom, atomList, xyz, cell_xyz, cell_type, fxyz, energy, stress, export_stress)
