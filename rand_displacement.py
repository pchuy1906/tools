import numpy as np
import random

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",         default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--dxyz", type=float, default=0.1,        help='delta XYZ for random displacement (in Angstrom)')
parser.add_argument("--cell_option",      default='cell_3',   help='cell_3/cell_9/NON_ORTHO')


args        = parser.parse_args()
file_xyz     = args.file_xyz
dxyz         = args.dxyz
cell_option  = args.cell_option

print ("")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

box = np.zeros(shape=(3,3))
tmp = f.readline().split()
try:
    # format: box[1,:] box[2,:] box[3,:]
    k = 0
    for k1 in range(3):
        for k2 in range(3):
            box[k1,k2] = float(tmp[k])
            k += 1
except:
    try:
        # format: NON_ORTHOR box[1,:] box[2,:] box[3,:]
        k = 1
        for k1 in range(3):
            for k2 in range(3):
                box[k1,k2] = float(tmp[k])
                k += 1
    except:
        box = np.zeros(shape=(3,3))
        # Orthor cell a, b, c
        for k1 in range(3):
            box[k1,k1] = float(tmp[k1])

print (box)
atomList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
    tmp = f.readline()
    tmp = tmp.split()
    atomList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

asym = np.array(atomList)
asym_unique = np.unique(asym)
asym_unique = np.sort(asym_unique)
asym_list = ' '.join(asym_unique)

f2 = open( "file_out.xyz", "w")
f2.write("%4d\n" %( natom))

if (cell_option=="cell_3"):
    f2.write("%15.9f %15.9f %15.9f  %d %s\n" %( box[0,0], box[1,1], box[2,2] ))
elif (cell_option=="cell_9"):
    f2.write("%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %( box[0,0], box[0,1], box[0,2], box[1,0], box[1,1], box[1,2], box[2,0], box[2,1], box[2,2] ))
else:
    print ("unknown cell_option ", cell_option)
    exit()


for k in range(0,natom):
    tmpp = np.zeros(shape=(3))
    for k1 in range(3):
        rdn = random.random()
        #print (rdn)
        tmpp[k1]  = xyz[k,k1] + dxyz *(2.0*rdn-1.0)
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( atomList[k], tmpp[0], tmpp[1], tmpp[2] ))
f2.close()


