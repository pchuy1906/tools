import numpy as np

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",              default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--option_cell", type=int, default=1,          help='1-(box3)  2-(box9) 3-(box9:NON_ORTHO)')

args        = parser.parse_args()
file_xyz      = args.file_xyz
option_cell   = args.option_cell

print ("")
print ("xyz_2_POSCAR")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

box = np.zeros(shape=(3,3))
tmp = f.readline().split()
if option_cell==1:
    for k in range(3):
        box[k,k] = float(tmp[k])
elif option_cell==2:
    tmp_cell = [float(tmp[i]) for i in range(9)]
    tmp_cell = np.array(tmp_cell)
    box = tmp_cell.reshape((3, 3))
elif option_cell==3:
    tmp_cell = [float(tmp[i]) for i in range(1,10)]
    tmp_cell = np.array(tmp_cell)
    box = tmp_cell.reshape((3, 3))
else:
    print ("Not implemented yet!", option_cell)
    exit()

    
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

f2 = open("POSCAR", "w")
f2.write("%1s\n" %( "COMMENT" ))
f2.write("%15.9f\n" %( 1.0 ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[0,0], box[0,1], box[0,2] ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[1,0], box[1,1], box[1,2] ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[2,0], box[2,1], box[2,2] ))

syms, counts_syms = np.unique(atomList, return_counts=True)

print (syms)
print (counts_syms)

f3 = open("ntype.dat", "w")
for item in syms:
    f2.write("%s %s" %(" ", item))
    f3.write("%s\n" %( item))

f2.write("\n")
for item in counts_syms:
    f2.write("%s %s" %(" ", item))
f2.write("\n")
f2.write("%1s\n" %( "Direct" ))


xyz = np.dot(xyz, np.linalg.inv(box))

for item in syms:
    for k in range(0,natom):
        if (item == atomList[k]):
            f2.write("%20.15f %20.15f %20.15f\n" %( xyz[k,0], xyz[k,1], xyz[k,2] ))

f2.close()
