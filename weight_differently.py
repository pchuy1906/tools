import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate weights')
# Arguments supported by the code.
parser.add_argument("--nAtomType",     type=int,   default=4)
parser.add_argument("--nCondensed",    type=int,   default=10)
parser.add_argument("--wE",            type=float, default=200)
parser.add_argument("--wS",            type=float, default=500)
parser.add_argument('--natom_special', nargs='+', type=int)

args    = parser.parse_args()
nAtomType  = args.nAtomType
nCondensed = args.nCondensed
wE         = args.wE
wS         = args.wS
natom_special = args.natom_special

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

        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )


    fE_extra = 1.0
    if (natom in natom_special):
        fE_extra = 5.0

    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))

