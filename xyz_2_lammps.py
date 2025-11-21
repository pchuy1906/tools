import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                         default='file.xyz',                 help='file_XYZ format xyz')
parser.add_argument('--atom_type', nargs='+')
parser.add_argument('--atom_mass', nargs='+', type=float)
parser.add_argument("--cell_type",                        default='cell_3',                   help='cell_3, cell_9, NON_ORTHO')
parser.add_argument("--use_molecule_id",      type=int,   default=1,                          help='0, 1 (default 1 yes for ChIMES)')
parser.add_argument("--triclinic_cell",                   default=False, action="store_true", help="option for triclinic cell")

args        = parser.parse_args()
file_xyz        = args.file_xyz
atom_type       = args.atom_type
atom_mass       = args.atom_mass
cell_type       = args.cell_type
use_molecule_id = args.use_molecule_id
triclinic_cell  = args.triclinic_cell

print ("")
print ("convert xyz to lammps")
print ("")
print ("read fileXYZ:", file_xyz)

natom, cell_3_3, atomList, xyz = utils.read_xyz(file_xyz, cell_type)
print (cell_3_3)

syms, counts_syms = np.unique(atomList, return_counts=True)
ntype = len(syms)

atom_type = np.array(atom_type)

atom_type_used = atom_type
atom_mass_used = atom_mass
ntype = len(atom_type)

f2 = open("data.lammps_new", "w")
f2.write("%1s\n" %( "# Position data file" ))
f2.write("%1s\n" %( "" ))
f2.write("%1d %1s\n" %( natom, "atoms" ))
f2.write("%1d %1s\n" %( ntype, "atom types" ))
f2.write("%1s\n" %( "" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, cell_3_3[0,0], "xlo xhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, cell_3_3[1,1], "ylo yhi" ))
Lz = cell_3_3[2,2]
f2.write("%15.9f %15.9f %1s\n" %( 0.0, Lz, "zlo zhi" ))

if triclinic_cell:
    f2.write("%15.9f %15.9f %15.9f %1s\n" %( cell_3_3[1,0], cell_3_3[2,0], cell_3_3[2,1], "xy xz yz" ))

f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Masses" ))
f2.write("%1s\n" %( "" ))

f3 = open("atoms.dat", "w")
for k in range(len(atom_type_used)):
    print (k+1, atom_mass_used[k])
    f2.write("%1d %15.9f\n" %( k+1, atom_mass_used[k] ))
    f3.write("%s %s" %( atom_type_used[k], " "))

f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Atoms" ))
f2.write("%1s\n" %( "" ))

for k in range(0,natom):
    for isym in range(len(atom_type_used)):
        if (atomList[k] == atom_type_used[isym]):
            atype = isym + 1

    if use_molecule_id == 1:
        f2.write("%1d %1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, 0, xyz[k,0], xyz[k,1], xyz[k,2] ))
    else:
        f2.write("%1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, xyz[k,0], xyz[k,1], xyz[k,2] ))

f2.close()
