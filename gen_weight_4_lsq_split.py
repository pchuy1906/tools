import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate weights')
# Arguments supported by the code.
parser.add_argument("--nAtomType",  type=int,   default=4)
parser.add_argument("--nCondensed", type=int,   default=10)
parser.add_argument("--wE",         type=float, default=200)
parser.add_argument("--wS",         type=float, default=500)
parser.add_argument("--file_nxyz",  type=str,   default="NEW.xyz")

args    = parser.parse_args()
nAtomType  = args.nAtomType
nCondensed = args.nCondensed
wE         = args.wE
wS         = args.wS
file_nxyz  = args.file_nxyz

import os.path
check_file = os.path.isfile(file_nxyz)

wS_vary = False

nconf = 0
if (check_file):
    wS_vary = True
    stress_collect = []
    f4 = open(file_nxyz, "rt")
    while True:
        tmp  = f4.readline()
        # if EOF, stop 
        line = tmp.strip()
        if line == '': break
        # otherwise, continue reading file
        tmp  = tmp.split()
        tmp_natom = int(tmp[0])

        tmp  = f4.readline().split()

        ncomment = len(tmp)
        if (ncomment==17):
            # format: NON_ORTHO cell(1,1:3), cell(2,1:3), cell(3,1:3), stress(1:6), energy
            cell_9 = [float(tmp[i]) for i in range(1,10)]
            stress = [float(tmp[i]) for i in range(10,16)]
        elif (ncomment==11):
            # format: NON_ORTHO cell(1,1:3), cell(2,1:3), cell(3,1:3), energy
            cell_9 = [float(tmp[i]) for i in range(1,10)]
            stress = []
        elif (ncomment==10):
            # format: cell(1), cell(2), cell(3), stress(1:6), energy
            cell_9 = [0.0 for i in range(1,10)]
            cell_9[0] = float(tmp[0])
            cell_9[4] = float(tmp[1])
            cell_9[8] = float(tmp[2])
            stress = [float(tmp[i]) for i in range(3,9)]
        elif (ncomment==4):
            # format: cell(1), cell(2), cell(3), energy
            cell_9 = [0.0 for i in range(1,10)]
            cell_9[0] = float(tmp[0])
            cell_9[4] = float(tmp[1])
            cell_9[8] = float(tmp[2])
            stress = []
        else:
            print ("unknown option")
            print ("ncomment=",ncomment)
            print (tmp)
            exit()

        stress_collect.append(stress)

        for i in range(tmp_natom):
            tmp  = f4.readline()

        nconf += 1

#print ("len(stress_collect)=", len(stress_collect))

f2 = open('label.txt', "w")
f3 = open('new_weight.dat', "w")

def gen_molecule_name(AType, ANumber):
    mole_name = ""
    nAtomType = len(AType)
    for i in range(nAtomType):
        mole_name = mole_name + AType[i] + str(ANumber[i])
    return mole_name





tmp_file = "frames.all.log"

f  = open(tmp_file ,"rt")
while True:

    tmp  = f.readline()
    # if EOF, stop 
    line = tmp.strip()
    if line == '': break
    # otherwise, continue reading file
    tmp  = tmp.split()
    id_frame = int(tmp[2])

    tmp  = f.readline().split()
    AType = []
    ANumber = []
    for j in range(nAtomType):
        AType.append(tmp[2*j])
        ANumber.append(int(tmp[2*j+1]))
    mole_name = gen_molecule_name(AType, ANumber)

    tmp  = f.readline()

    natom = sum(ANumber)

    for j in range(natom):
        f2.write("force_"+mole_name+"\n")
        f2.write("force_"+mole_name+"\n")
        f2.write("force_"+mole_name+"\n")
        f3.write("1.0\n")
        f3.write("1.0\n")
        f3.write("1.0\n")

    if (id_frame < nCondensed):
        f2.write("dia_stress_" +mole_name+"\n")
        f2.write("stress_xy_"  +mole_name+"\n")
        f2.write("stress_xz_"  +mole_name+"\n")
        f2.write("stress_xy_"  +mole_name+"\n")
        f2.write("dia_stress_" +mole_name+"\n")
        f2.write("stress_yz_"  +mole_name+"\n")
        f2.write("stress_xz_"  +mole_name+"\n")
        f2.write("stress_yz_"  +mole_name+"\n")
        f2.write("dia_stress_" +mole_name+"\n")
        if check_file:
            stress_values = stress_collect[id_frame]
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[0]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[3]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[4]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[3]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[1]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[5]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[4]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[5]/50)**2 )) )
            f3.write('%15.6f \n' %(wS*1.0/(1.0+(stress_values[2]/50)**2 )) )
        else:
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )
            f3.write('%15.6f \n' %(wS) )

    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f3.write('%15.6f \n' %(wE/float(natom)))
    f3.write('%15.6f \n' %(wE/float(natom)))
    f3.write('%15.6f \n' %(wE/float(natom)))

