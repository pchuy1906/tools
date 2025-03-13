import numpy as np

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                         default='file.xyz', help='file_XYZ format xyz')
parser.add_argument('--atom_type', nargs='+')
parser.add_argument('--atom_mass', nargs='+', type=float)
parser.add_argument("--num_cell_input",       type=int,   default=3, help='3, 9, 10')
parser.add_argument("--use_molecule_id",      type=int,   default=1, help='0, 1(default, yes,ChIMES)')
parser.add_argument("--triclinic_cell",                   default=False, action="store_true", help="option for triclinic cell")
parser.add_argument("--dipole_sphere",                    default=False, action="store_true", help="option for atom_style dipole sphere")

args        = parser.parse_args()
file_xyz        = args.file_xyz
atom_type       = args.atom_type
atom_mass       = args.atom_mass
num_cell_input  = args.num_cell_input
use_molecule_id = args.use_molecule_id
triclinic_cell  = args.triclinic_cell
dipole_sphere   = args.dipole_sphere

print ("")
print ("Orthohombic cell case")
print ("")
print ("need to check variable atom_type and atom_mass")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

cell_3_3 = np.zeros(shape=(3,3))
tmp = f.readline().split()

if num_cell_input==3:
    for k in range(3):
        cell_3_3[k,k] = float(tmp[k])
if num_cell_input==9:
    cell_3_3[0,0] = float(tmp[0])
    cell_3_3[0,1] = float(tmp[1])
    cell_3_3[0,2] = float(tmp[2])

    cell_3_3[1,0] = float(tmp[3])
    cell_3_3[1,1] = float(tmp[4])
    cell_3_3[1,2] = float(tmp[5])

    cell_3_3[2,0] = float(tmp[6])
    cell_3_3[2,1] = float(tmp[7])
    cell_3_3[2,2] = float(tmp[8])

if num_cell_input==10:

    cell_3_3[0,0] = float(tmp[1])
    cell_3_3[0,1] = float(tmp[2])
    cell_3_3[0,2] = float(tmp[3])

    cell_3_3[1,0] = float(tmp[4])
    cell_3_3[1,1] = float(tmp[5])
    cell_3_3[1,2] = float(tmp[5])

    cell_3_3[2,0] = float(tmp[7])
    cell_3_3[2,1] = float(tmp[8])
    cell_3_3[2,2] = float(tmp[9])

myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

syms, counts_syms = np.unique(myList, return_counts=True)
ntype = len(syms)

atom_type = np.array(atom_type)

atom_type_used = atom_type
atom_mass_used = atom_mass
ntype = len(atom_type)

#atom_type_used = []
#atom_mass_used = []
#
#for i in range(ntype):
#    tsym = syms[i]
#    idx  = np.where(atom_type==tsym)[0][0]
#    atom_type_used.append(atom_type[idx])
#    atom_mass_used.append(atom_mass[idx])
#
#print (atom_type_used)
#print (atom_mass_used)

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
        if (myList[k] == atom_type_used[isym]):
            atype = isym + 1

    if not dipole_sphere:
        if use_molecule_id == 1:
            f2.write("%1d %1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, 0, xyz[k,0], xyz[k,1], xyz[k,2] ))
        else:
            f2.write("%1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, xyz[k,0], xyz[k,1], xyz[k,2] ))
    else:
        f2.write("%6d" %( k+1 ))
        f2.write("%6d" %( atype ))
        for i in range(3):
            f2.write("%15.9f" %( xyz[k,i]))
        #charge
        f2.write("%6d" %( 0 ))
        #muxyz
        for i in range(3):
            f2.write("%6d" %( 0 ))
        #diameter density
        f2.write("%6d" %( 1 ))
        f2.write("%6d" %( 1 ))
        f2.write("\n")

f2.close()
