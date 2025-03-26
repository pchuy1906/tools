import numpy as np
import os

import argparse
parser = argparse.ArgumentParser(description='setup DFTB input for file_xyz')
# Arguments supported by the code.
parser.add_argument("--file1_xyz",                  default='file1.xyz', help='file_xyz format xyz')
parser.add_argument("--file2_xyz",                  default='file2.xyz', help='file_xyz format xyz')
parser.add_argument("--cell_option_input",          default='cell_3',    help='cell_option_input  = "NON_ORTHO", "cell_3", or "cell_9" ')
parser.add_argument("--cell_option_output",         default='cell_3',    help='cell_option_output = "NON_ORTHO", "cell_3", or "cell_9" ')
parser.add_argument("--export_quan",                default='xyzfes',    help='export xyzfe/xyzfes')


args        = parser.parse_args()
file1_xyz            = args.file1_xyz
file2_xyz            = args.file2_xyz
cell_option_input    = args.cell_option_input
cell_option_output   = args.cell_option_output
export_quan          = args.export_quan


Ha2kcalmol = 627.503
Bohr2Ang   = 0.529177

f1  = open(file1_xyz ,"rt")
f2  = open(file2_xyz ,"rt")
f3 = open("ndelta.xyz", "w")
f4 = open("data_file1.dat", "w")
f5 = open("data_file2.dat", "w")
f6 = open("data_natom.dat", "w")

nframe = 0
while True:
    #print ("frame", nframe)
    tmp1  = f1.readline()
    line = tmp1.strip()
    if line == '': break
    #v, e = [float(x) for x in line.split()[:2]]

    # line-1: number of atoms
    natom = int(tmp1)
    tmp2  = f2.readline()
    f3.write("%-d\n" %( natom))
    
    # line-2: cell_parameter, stress tensor and energy
    tmp1 = f1.readline().split()
    tmp2 = f2.readline().split()

    if (cell_option_input == "NON_ORTHO"):
        if (export_quan == "xyzfes"):
            cell1_9 = [float(x) for x in tmp1[1:10]]
            stress1_6 = [float(x) for x in tmp1[10:16]]
            energy1 = float(tmp1[16])
            stress2_6 = [float(x) for x in tmp2[10:16]]
            energy2 = float(tmp2[16])
    elif (cell_option_input == "cell_3"):
        if (export_quan == "xyzfes"):
            cell1_9 = [0 for i in range(9)]
            cell1_9[0] = float(tmp1[0])
            cell1_9[4] = float(tmp1[1])
            cell1_9[8] = float(tmp1[2])
            energy1 = float(tmp1[9])
            energy2 = float(tmp2[9])
            stress1_6 = [float(tmp1[i]) for i in range(3,9)]
            stress2_6 = [float(tmp2[i]) for i in range(3,9)]
        elif (export_quan == "xyzfe"):
            cell1_9 = [0 for i in range(9)]
            cell1_9[0] = float(tmp1[0])
            cell1_9[4] = float(tmp1[1])
            cell1_9[8] = float(tmp1[2])
            energy1 = float(tmp1[3])
            energy2 = float(tmp2[3])
            stress1_6 = [0 for i in range(6)]
            stress2_6 = [0 for i in range(6)]
        else:
            print ("Not implemented yet!")
            exit()
    else:
        print ("Not implemented yet!")
        exit()

    dstress = np.subtract(stress1_6, stress2_6)
    denergy = energy1 - energy2

    if (cell_option_output == "cell_3"):
        if (export_quan == "xyzfes"):
            f3.write("%12.6f" %( cell1_9[0] ))
            f3.write("%12.6f" %( cell1_9[4] ))
            f3.write("%12.6f" %( cell1_9[8] ))
            for i in range(6):
                f3.write("%12.6f" %( dstress[i] ))
            f3.write("%15.6f\n" %( denergy ))
        elif (export_quan == "xyzfe"):
            f3.write("%12.6f" %( cell1_9[0] ))
            f3.write("%12.6f" %( cell1_9[4] ))
            f3.write("%12.6f" %( cell1_9[8] ))
            f3.write("%15.6f\n" %( denergy ))
        else:
            print ("Not implemented yet!")
            exit()
    elif (cell_option_output == "NON_ORTHO"):
        if (export_quan == "xyzfes"):
            f3.write("%s" %( "NON_ORTHO" ))
            for i in range(9):
                f3.write("%12.6f" %( cell1_9[i] ))
            for i in range(6):
                f3.write("%12.6f" %( dstress[i] ))
            f3.write("%15.6f\n" %( denergy ))
        elif (export_quan == "xyzfe"):
            f3.write("%s" %( "NON_ORTHO" ))
            for i in range(9):
                f3.write("%12.6f" %( cell1_9[i] ))
            f3.write("%15.6f\n" %( denergy ))
    else:
        print ("Not implemented yet!")
        exit()

    # line-next: xyz coordinate
    atomList1 = []
    for k in range(natom):
        tmp1 = f1.readline().split()
        tmp2 = f2.readline().split()
        atomList1.append(tmp1[0])
        xyz1 = [float(x) for x in tmp1[1:4]]
        fxyz1 = [float(x) for x in tmp1[4:7]]
        fxyz2 = [float(x) for x in tmp2[4:7]]
        #print (len(fxyz1), len(fxyz2))
        #print (fxyz1)
        #print (fxyz2)
        dfxyz = np.subtract(fxyz1, fxyz2)
        f3.write("%s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" %(tmp1[0], xyz1[0], xyz1[1], xyz1[2], dfxyz[0], dfxyz[1], dfxyz[2] ))
        for ixyz in range(3):
            f4.write("%12.6f %s\n" %( fxyz1[ixyz] * Ha2kcalmol/Bohr2Ang, "FORCE_"+str(nframe)+"_"+tmp1[0] ))
            f5.write("%12.6f %s\n" %( fxyz2[ixyz] * Ha2kcalmol/Bohr2Ang, "FORCE_"+str(nframe)+"_"+tmp1[0] ))

    f6.write("%d \n" %( natom ))
    f4.write("%12.6f %s\n" %( energy1, "ENERGY_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( energy2, "ENERGY_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[0], "STRESS_XX_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[1], "STRESS_YY_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[2], "STRESS_ZZ_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[3], "STRESS_XY_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[4], "STRESS_XZ_"+str(nframe) ))
    f4.write("%12.6f %s\n" %( stress1_6[5], "STRESS_YZ_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[0], "STRESS_XX_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[1], "STRESS_YY_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[2], "STRESS_ZZ_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[3], "STRESS_XY_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[4], "STRESS_XZ_"+str(nframe) ))
    f5.write("%12.6f %s\n" %( stress2_6[5], "STRESS_YZ_"+str(nframe) ))

    nframe += 1
