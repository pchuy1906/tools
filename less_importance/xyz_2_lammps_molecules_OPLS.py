import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_LAMMPS_data", default='data.lammps', help='file LAMMPS data')
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--nmolecule", type=int, default=1, help='the number of molecules')
parser.add_argument("--scale_q", type=float, default=1.0, help='scale partial charges')
parser.add_argument("--scale_eps", type=float, default=1.0, help='scale epsilon LJ')
parser.add_argument("--ishock", type=int, default=0, help='0-no, 1-shock')

args                = parser.parse_args()
file_LAMMPS_data    = args.file_LAMMPS_data
file_xyz            = args.file_xyz
nmolecule           = args.nmolecule
scale_q             = args.scale_q
scale_eps           = args.scale_eps
ishock              = args.ishock

print ("Generate the LAMMPS data file")
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

zmin = np.min(xyz[:,2])
zmax = np.max(xyz[:,2])
zboxmaxmin = zmax-zmin+0.02



f  = open(file_LAMMPS_data ,"r")
f2 = open("NEW_data.lammps", "w")
f2.write("%s\n" %( "Generate data LAMMPS" ))
f2.write("%s\n" %( "" ))


while True:
    line = f.readline().split()

    #if (len(line)==1) and (line[0] == 'Impropers') : break

    if (len(line)==2) and (line[1] == 'atoms'):
        natoms = int(line[0])
        f2.write("%1d %1s\n" %( natoms * nmolecule, " atoms" ))

    if (len(line)==2) and (line[1] == 'bonds'):
        nbonds = int(line[0])
        f2.write("%1d %1s\n" %( nbonds * nmolecule, " bonds" ))

    if (len(line)==2) and (line[1] == 'angles'):
        nangles = int(line[0])
        f2.write("%1d %1s\n" %( nangles * nmolecule, " angles" ))

    if (len(line)==2) and (line[1] == 'dihedrals'):
        ndihedrals = int(line[0])
        f2.write("%1d %1s\n" %( ndihedrals * nmolecule, " dihedrals" ))

    if (len(line)==2) and (line[1] == 'impropers'):
        nimpropers = int(line[0])
        f2.write("%1d %1s\n" %( nimpropers * nmolecule, " impropers" ))
        f2.write("%s\n" %( "" ))


    ################################################################

    if (len(line)==3) and (line[2] == 'types') and (line[1] == 'atom'):
        ntypes_atoms = int(line[0])
        f2.write("%1d %1s\n" %( ntypes_atoms , " atom types" ))

    if (len(line)==3) and (line[2] == 'types') and (line[1] == 'bond'):
        ntypes_bonds = int(line[0])
        f2.write("%1d %1s\n" %( ntypes_bonds , " bond types" ))

    if (len(line)==3) and (line[2] == 'types') and (line[1] == 'angle'):
        ntypes_angles = int(line[0])
        f2.write("%1d %1s\n" %( ntypes_angles , " angle types" ))

    if (len(line)==3) and (line[2] == 'types') and (line[1] == 'dihedral'):
        ntypes_dihedrals = int(line[0])
        f2.write("%1d %1s\n" %( ntypes_dihedrals , " dihedral types" ))

    if (len(line)==3) and (line[2] == 'types') and (line[1] == 'improper'):
        ntypes_impropers = int(line[0])
        f2.write("%1d %1s\n" %( ntypes_impropers , " improper types" ))
        f2.write("%s\n" %( "" ))
        f2.write("%15.9f %15.9f %s\n" %( 0.0,  box[0,0], " xlo xhi" ))
        f2.write("%15.9f %15.9f %s\n" %( 0.0,  box[1,1], " ylo yhi" ))
        if (ishock==0):
            f2.write("%15.9f %15.9f %s\n" %( 0.0,  box[2,2], " zlo zhi" ))
        if (ishock==1):
            f2.write("%15.9f %15.9f %s\n" %( 0.0,  zboxmaxmin, " zlo zhi" ))

        f2.write("%s\n" %( "" ))



    if (len(line)==1) and (line[0] == 'Masses') :
        line = f.readline()
        f2.write("%s\n" %( "Masses" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_atoms):
            line = f.readline()
            f2.write("%s" %( line ))
        f2.write("%s\n" %( "" ))

    if (len(line)==2) and (line[0] == 'Pair') and (line[1] == 'Coeffs') :
        line = f.readline()
        f2.write("%s\n" %( "Pair Coeffs" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_atoms):
            line = f.readline().split()
            f2.write("%s %12.6f %s\n" %( line[0], float(line[1]) * scale_eps, line[2] ))
        f2.write("%s\n" %( "" ))

    if (len(line)==2) and (line[0] == 'Bond') and (line[1] == 'Coeffs') :
        line = f.readline()
        f2.write("%s\n" %( "Bond Coeffs" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_bonds):
            line = f.readline()
            f2.write("%s" %( line ))
        f2.write("%s\n" %( "" ))

    if (len(line)==2) and (line[0] == 'Angle') and (line[1] == 'Coeffs') :
        line = f.readline()
        f2.write("%s\n" %( "Angle Coeffs" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_angles):
            line = f.readline()
            f2.write("%s" %( line ))
        f2.write("%s\n" %( "" ))

    if (len(line)==2) and (line[0] == 'Dihedral') and (line[1] == 'Coeffs') :
        line = f.readline()
        f2.write("%s\n" %( "Dihedral Coeffs" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_dihedrals):
            line = f.readline()
            f2.write("%s" %( line ))
        f2.write("%s\n" %( "" ))

    if (len(line)==2) and (line[0] == 'Improper') and (line[1] == 'Coeffs') :
        line = f.readline()
        f2.write("%s\n" %( "Improper Coeffs" ))
        f2.write("%s\n" %( "" ))
        for k in range(ntypes_impropers):
            line = f.readline()
            f2.write("%s" %( line ))
        f2.write("%s\n" %( "" ))



    if (len(line)==1) and (line[0] == 'Atoms') :
        line = f.readline()
        f2.write("%s\n" %( "Atoms" ))
        f2.write("%s\n" %( "" ))

        a_type = []
        a_q = []
        for k in range(ntypes_atoms):
            line = f.readline().split()
            a_type.append(int(line[2]))
            a_q.append(float(line[3]))
        print ("number of charges:", len(a_q))
        print ("sum charges:", sum(a_q))
        a_q = a_q - sum(a_q)/len(a_q)
        print ("number of charges:", len(a_q))
        print ("sum charges:", sum(a_q))

        for k in range(natom):
            idx = k%natom_per_mol
            idm = k//natom_per_mol+1
            if (ishock==0):
                f2.write("%5d %5d %5d %15.12f %15.9f %15.9f %15.9f\n" %( k+1, idm,  a_type[idx], a_q[idx] * scale_q, xyz[k,0], xyz[k,1], xyz[k,2]  ))
            if (ishock==1):
                f2.write("%5d %5d %5d %15.12f %15.9f %15.9f %15.9f\n" %( k+1, idm,  a_type[idx], a_q[idx] * scale_q, xyz[k,0], xyz[k,1], xyz[k,2]-zmin+0.01  ))
        f2.write("%s\n" %( "" ))



    if (len(line)==1) and (line[0] == 'Bonds') :
        line = f.readline()
        f2.write("%s\n" %( "Bonds" ))
        f2.write("%s\n" %( "" ))

        _type = []
        i_1 = []
        i_2 = []
        for k in range(ntypes_bonds):
            line = f.readline().split()
            _type.append(int(line[1]))
            i_1.append(int(line[2]))
            i_2.append(int(line[3]))

        ncount = 0
        for imol in range(nmolecule):
            base = imol * natom_per_mol
            for k in range(ntypes_bonds):
                ncount += 1
                f2.write("%5d %5d %5d %5d \n" %( ncount, _type[k], i_1[k]+base, i_2[k]+base ))
        f2.write("%s\n" %( "" ))



    if (len(line)==1) and (line[0] == 'Angles') :
        line = f.readline()
        f2.write("%s\n" %( "Angles" ))
        f2.write("%s\n" %( "" ))

        _type = []
        i_1 = []
        i_2 = []
        i_3 = []
        for k in range(ntypes_angles):
            line = f.readline().split()
            _type.append(int(line[1]))
            i_1.append(int(line[2]))
            i_2.append(int(line[3]))
            i_3.append(int(line[4]))

        ncount = 0
        for imol in range(nmolecule):
            base = imol * natom_per_mol
            for k in range(ntypes_angles):
                ncount += 1
                f2.write("%5d %5d %5d %5d %5d \n" %( ncount, _type[k], i_1[k]+base, i_2[k]+base, i_3[k]+base ))
        f2.write("%s\n" %( "" ))


    if (len(line)==1) and (line[0] == 'Dihedrals') :
        line = f.readline()
        f2.write("%s\n" %( "Dihedrals" ))
        f2.write("%s\n" %( "" ))

        _type = []
        i_1 = []
        i_2 = []
        i_3 = []
        i_4 = []
        for k in range(ntypes_dihedrals):
            line = f.readline().split()
            _type.append(int(line[1]))
            i_1.append(int(line[2]))
            i_2.append(int(line[3]))
            i_3.append(int(line[4]))
            i_4.append(int(line[5]))

        ncount = 0
        for imol in range(nmolecule):
            base = imol * natom_per_mol
            for k in range(ntypes_dihedrals):
                ncount += 1
                f2.write("%5d %5d %5d %5d %5d %5d \n" %( ncount, _type[k], i_1[k]+base, i_2[k]+base, i_3[k]+base, i_4[k]+base ))
        f2.write("%s\n" %( "" ))


    if (len(line)==1) and (line[0] == 'Impropers') :
        line = f.readline()
        f2.write("%s\n" %( "Impropers" ))
        f2.write("%s\n" %( "" ))

        _type = []
        i_1 = []
        i_2 = []
        i_3 = []
        i_4 = []
        for k in range(ntypes_impropers):
            line = f.readline().split()
            _type.append(int(line[1]))
            i_1.append(int(line[2]))
            i_2.append(int(line[3]))
            i_3.append(int(line[4]))
            i_4.append(int(line[5]))

        ncount = 0
        for imol in range(nmolecule):
            base = imol * natom_per_mol
            for k in range(ntypes_impropers):
                ncount += 1
                f2.write("%5d %5d %5d %5d %5d %5d \n" %( ncount, _type[k], i_1[k]+base, i_2[k]+base, i_3[k]+base, i_4[k]+base ))
        f2.write("%s\n" %( "" ))
        break

