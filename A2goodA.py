import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np
from math import floor

import argparse
parser = argparse.ArgumentParser(description='xyz_2_good_xyz where all atoms are in the box')
# Arguments supported by the code.
parser.add_argument("--file_xyz",    default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--cell_type",   default='cell_3',     help='cell_3, cell_9, NON_ORTHO')

args        = parser.parse_args()
file_xyz        = args.file_xyz
cell_type  = args.cell_type

natom, cell_3_3, atomList, xyz = utils.read_xyz(file_xyz, cell_type)

numrows = xyz.shape[0]
numcols = xyz.shape[1]
print ("size of xyz:", numrows, numcols)

print ("")
print ("Angstrom to Crystal calculation...")
Cxyz = np.dot(xyz, np.linalg.inv(cell_3_3))

for k1 in range(0,numrows):
    for k2 in range(0,numcols):
        tmp = Cxyz[k1][k2]-floor(Cxyz[k1][k2])
        Cxyz[k1][k2] = tmp

print ("Crystal to Angstrom calculation...")
Axyz = np.dot(Cxyz, cell_3_3)

print ("")
print ("Output: <Axyz.xyz>")
f2 = open('Axyz.xyz', "w")

f2.write("%-d\n" %(natom))

print (cell_3_3)
for ixyz in range(3):
    f2.write("  %15.9f" % cell_3_3[ixyz,ixyz])
f2.write("\n")

for k in range(0,natom):
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( atomList[k], Axyz[k,0], Axyz[k,1], Axyz[k,2] ))
f2.close
