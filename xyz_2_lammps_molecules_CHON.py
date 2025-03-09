import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_para", default='file.xyz', help='file parameter')
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--nmolecule", type=int, default=1, help='the number of molecules')

args        = parser.parse_args()
file_para    = args.file_para
file_xyz     = args.file_xyz
nmolecule    = args.nmolecule


print ("")
print ("Orthohombic cell case")
print ("")
print ("need to check variable data_sym and data_mass")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

natom_per_mol = natom / nmolecule


box = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
for k in range(3):
    box[k,k] = float(tmp[k])

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

f2 = open("data.lammps", "w")
f2.write("%1s\n" %( "# Position data file" ))
f2.write("%1s\n" %( "" ))
f2.write("%1d %1s\n" %( natom, "atoms" ))
f2.write("%1d %1s\n" %( ntype, "atom types" ))
f2.write("%1s\n" %( "" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[0,0], "xlo xhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[1,1], "ylo yhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[2,2], "zlo zhi" ))
f2.write("%1s\n" %( "" ))

data_sym = np.loadtxt(file_para, usecols=0, dtype=str)
print (data_sym)
data_mass = np.loadtxt(file_para, usecols=1, dtype=float)
print (data_mass)
data_sig = np.loadtxt(file_para, usecols=2, dtype=float)
print (data_sig)
data_eps = np.loadtxt(file_para, usecols=3, dtype=float)
print (data_eps)

f2.write("%1s\n" %( "Masses" ))
f2.write("%1s\n" %( "" ))
for k in range(len(data_sym)):
    print (k+1, data_mass[k])
    f2.write("%1d %15.9f\n" %( k+1, data_mass[k] ))
f2.write("%1s\n" %( "" ))

#f2.write("%1s\n" %( "Pair Coeffs" ) )
#f2.write("%1s\n" %( "" ))
#rcut = 8.0
#for k1 in range(len(data_sig)):
#    for k2 in range(k1,len(data_sig)):
#        eps = np.sqrt(data_eps[k1] * data_eps[k2])
#        sig =  0.5 * (data_sig[k1] + data_sig[k2])
#        f2.write("%d %d %15.9f %15.9f %15.9f\n" %( k1+1, k2+1, eps, sig, rcut ))
#f2.write("%1s\n" %( "" ))

f2.write("%1s\n" %( "Atoms" ))
f2.write("%1s\n" %( "" ))
for k in range(0,natom):
    for isym in range(len(data_sym)):
        if (myList[k] == data_sym[isym]):
            atype = isym + 1
    molecule_tag = k//natom_per_mol + 1
    #f2.write("%1d %4d %1d %1d %15.9f %15.9f %15.9f\n" %( k+1, molecule_tag, atype, 0, xyz[k,0], xyz[k,1], xyz[k,2] ))
    f2.write("%d %4d %15.9f %15.9f %15.9f\n" %( k+1, atype, xyz[k,0], xyz[k,1], xyz[k,2] ))

f2.close()


f4 = open("in.lammps", "w")
nmol = natom / natom_per_mol

f4.write("%s\n" %("units	real"))
f4.write("%s\n" %("atom_style	atomic"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("read_data    data.lammps"))
f4.write("%s\n" %(""))


f4.write("%s\n" %("pair_style   lj/cut 8.0"))
f4.write("%1s\n" %( "" ))
rcut = 8.0
for k1 in range(len(data_sig)):
    for k2 in range(k1,len(data_sig)):
        eps = np.sqrt(data_eps[k1] * data_eps[k2])
        sig =  0.5 * (data_sig[k1] + data_sig[k2])
        f4.write("%s %d %d %15.9f %15.9f %15.9f\n" %( "pair_coeff", k1+1, k2+1, eps, sig, rcut ))
f4.write("%1s\n" %( "" ))


f4.write("%s\n" %("velocity 	all create 300.0 4928459"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("# unconnected bodies"))
f4.write("%s\n" %(""))

sym="fix 1 all rigid/nvt group " + str(nmol)
for k in range(nmol):
    id_begin = 1 + k*natom_per_mol
    id_end   = (1 + k)*natom_per_mol
    f4.write("%s %s %d %d\n" %("group	clump" + str(k+1)," id <> ", id_begin, id_end))
    sym = sym + " clump" + str(k+1)
sym = sym + " &"

f4.write("%s\n" %(""))
f4.write("%s\n" %(sym))
f4.write("%s\n" %("                      temp 300.0 300.0 5.0 reinit no"))

f4.write("%s\n" %(""))
for k in range(nmol):
    sym = "neigh_modify exclude group clump"+ str(k+1)+ " clump" + str(k+1)
    f4.write("%s\n" %(sym))

f4.write("%s\n" %(""))
f4.write("%s\n" %("thermo		5000"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("dump             dump_1 all custom 5000 dump_xyz_vxyz id type x y z vx vy vz"))
f4.write("%s\n" %(""))
f4.write("%s\n" %(""))
f4.write("%s\n" %("timestep 	0.1"))
f4.write("%s\n" %("thermo		5000"))
f4.write("%s\n" %("run		10000000"))

f4.close()

