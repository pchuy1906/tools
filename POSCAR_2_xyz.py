import numpy as np

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",           default="POSCAR", help=' POSCAR/CONTCAR')
parser.add_argument("--option_cell",  type=int, default=2,        help=' 1-(box3)  2-(box9) 3-(box9:NON_ORTHO)')

args = parser.parse_args()
file_POSCAR  = args.file_POSCAR
option_cell   = args.option_cell

f  = open(file_POSCAR ,"r")
tmp = f.readline()
tmp = f.readline()
fscale = float(tmp.split()[0])
print ("scale-factor = ", fscale)

unitcell = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
unitcell[0,:] = tmp
tmp = f.readline()
tmp = tmp.split()
unitcell[1,:] = tmp
tmp = f.readline()
tmp = tmp.split()
unitcell[2,:] = tmp

unitcell = unitcell*fscale
print ("unitcell=")
print (unitcell)

tmp = f.readline()
atomNameList = tmp.split()
print ("atomNameList=",atomNameList)
tmp = f.readline().split()
atomNumList = [int(x) for x in tmp]
print ("atomNumList=",atomNumList)

natom = sum(atomNumList)
print ("natom=", natom)
 
tmp = f.readline().split()
if tmp[0]=="Cartesian":
    xyz_type = "C"
if tmp[0]=="Direct":
    xyz_type = "D"

xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    xyz[k,:] =  tmp[:3]
f.close()

if xyz_type == "D":
    #print (xyz)
    xyz = np.dot(xyz, unitcell)

f2 = open('POSCAR2xyz.xyz', "w")
f2.write("%-d\n" %(natom))

if (option_cell==1):
    print (unitcell)
    for i in range(3):
        f2.write("  %15.9f" % unitcell[i,i])
    f2.write("\n")
elif (option_cell==2):
    unitcell = unitcell.reshape(9)
    print (unitcell)
    for icellxyz in unitcell:
        f2.write("  %15.9f" % icellxyz)
    f2.write("\n")
elif (option_cell==3):
    unitcell = unitcell.reshape(9)
    print (unitcell)
    for icellxyz in unitcell:
        f2.write("  %15.9f" % icellxyz)
    f2.write("\n")



ixyz = 0
for i in range(len(atomNameList)):
    for j in range(atomNumList[i]):
        f2.write("%-s %15.9f %15.9f %15.9f\n" %(atomNameList[i], xyz[ixyz,0], xyz[ixyz,1], xyz[ixyz,2] ))
        ixyz += 1
f2.close()
