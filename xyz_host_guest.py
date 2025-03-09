import numpy as np

import argparse
parser = argparse.ArgumentParser(description='export host and guest seperately')
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

unitcell = unitcell.reshape(9)
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

sym_join = []
for i in range(natom-6):
    tsym = []
    for j in range(6):
        tsym.append(atomList[i+j])
    tsym_join = ''.join(tsym)
    sym_join.append(tsym_join)


sym_join = np.array(sym_join)
select_indices = np.where(sym_join=="HHOHHO")
print  (select_indices[0][0])
idcut = select_indices[0][0]

#syms, counts_syms = np.unique(atomList, return_counts=True)
#print (syms)
#asym_list = ' '.join(syms)
#print (asym_list)
#print (counts_syms)

f2 = open('host.xyz', "w")
f2.write("%-d\n" %(idcut))
for i in range(9):
    f2.write("%15.9f" %(unitcell[i]))
f2.write("\n" )
for k in range(idcut):
    f2.write("%-5s %15.9f %15.9f %15.9f\n" %(atomList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
f2.close()


f2 = open('guest.xyz', "w")
f2.write("%-d\n" %(natom-idcut))
for i in range(9):
    f2.write("%15.9f" %(unitcell[i]))
f2.write("\n" )
for k in range(idcut,natom):
    f2.write("%-5s %15.9f %15.9f %15.9f\n" %(atomList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
f2.close()

