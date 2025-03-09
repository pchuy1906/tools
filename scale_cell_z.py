import numpy as np

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",   default='file.xyz',     help='file_XYZ format xyz')
parser.add_argument("--file_ncell", default='new_cell.dat', help='file new cell parameters')

print ("This work only for orthorgonal cell")
print ("Requirements:")
print ("1. Second line of file_xyz contains cell lengths: a b c")
print ("2. File new_cell.dat contains cell lengths: a b c_new")

args        = parser.parse_args()
file_xyz     = args.file_xyz
file_ncell   = args.file_ncell

print ("")
print ("Orthohombic cell case")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

cell = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
for k in range(3):
    cell[k,k] = float(tmp[k])

myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

print ("old cell")
print (cell)
newcell = np.zeros(shape=(3,3))
ncell = np.loadtxt(file_ncell)
newcell[0,0] = ncell[0]
newcell[1,1] = ncell[1]
newcell[2,2] = ncell[2]

print ("new cell")
print (newcell)

Cxyz = np.dot(xyz, np.linalg.inv(cell))
Axyz = np.dot(Cxyz, newcell)

f2 = open("deform.xyz", "w")
f2.write("%4d\n" %( natom))
f2.write("%15.9f %15.9f %15.9f\n" %( newcell[0,0], newcell[1,1], newcell[2,2] ))
for k4 in range(0,natom):
    tmpp  = Axyz[k4,:]
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k4], tmpp[0], tmpp[1], tmpp[2] ))
f2.close()
