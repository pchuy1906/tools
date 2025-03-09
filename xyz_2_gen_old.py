import numpy as np

import argparse
parser = argparse.ArgumentParser(description='xyz_2_gen')
# Arguments supported by the code.
parser.add_argument("--file_input", default='file.xyz', help='xyz_2_gen')
parser.add_argument("--num_cell_input", type=int, default=3, help='3, 9, 10')

args        = parser.parse_args()
file_input     = args.file_input
num_cell_input = args.num_cell_input

print ("read fileXYZ:", file_input)
f  = open(file_input ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

unitcell = np.zeros(shape=(3,3))
tmp = f.readline().split()
tmp = np.array(tmp)

if num_cell_input==3:
    for k in range(3):
        unitcell[k,k] = float(tmp[k])
if num_cell_input==9:
    unitcell[0,:] = [float(x) for x in tmp[0:3]]
    unitcell[1,:] = [float(x) for x in tmp[3:6]]
    unitcell[2,:] = [float(x) for x in tmp[6:9]]
if num_cell_input==10:
    unitcell[0,:] = [float(x) for x in tmp[1:4]]
    unitcell[1,:] = [float(x) for x in tmp[4:7]]
    unitcell[2,:] = [float(x) for x in tmp[7:10]]



print (unitcell)

atomList = []
xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    atomList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

syms, counts_syms = np.unique(atomList, return_counts=True)

print (syms)
asym_list = ' '.join(syms)
print (asym_list)
print (counts_syms)

f2 = open('file.gen', "w")
f2.write("%-d %4s\n" %(natom, "S"))
f2.write("%-s \n" %(asym_list))

for k in range(natom):
    ind = np.where(syms == atomList[k])[0][0]+1
    f2.write("%-5d %5d %15.9f %15.9f %15.9f\n" %(k+1, ind, xyz[k,0], xyz[k,1], xyz[k,2] ))

f2.write("%15.9f%15.9f%15.9f\n" % (0.0,0.0,0.0))
for ixyz in range(3):
    for icellxyz in unitcell[ixyz,:]:
        f2.write("%15.9f" % icellxyz)
    f2.write("\n")

f2.close()

