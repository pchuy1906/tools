import os
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='correct stress tensor for file xyzf')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz           = args.file_xyz

print ("Calculate the atomic energies")
print ("")
print ("read fileXYZ:", file_xyz)

f  = open(file_xyz ,"rt")
f2 = open("corrected_" + file_xyz, "w")

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp.split()[0])
    f2.write("%s" %( tmp ))

    tmp  = f.readline().split()
    ## print the cell length first
    #f2.write(tmp[0]+ " " + tmp[1] + " " + tmp[2] + " ")
    ## then the stress tensor
    #f2.write(tmp[3]+ " " + tmp[4] + " " + tmp[5] + " " + tmp[6] + " " + tmp[8] + " " + tmp[7] + " ")
    ## then the energy
    #f2.write(tmp[9])
    #f2.write('\n')

    # print the cell first
    f2.write('NON_ORTHO ')
    for i in range(1,10):
        f2.write("%s " %( tmp[i] ))
    # then the stress tensor
    f2.write("%s " %( tmp[10] ))
    f2.write("%s " %( tmp[11] ))
    f2.write("%s " %( tmp[12] ))
    f2.write("%s " %( tmp[13] ))
    f2.write("%s " %( tmp[15] ))
    f2.write("%s " %( tmp[14] ))
    # then the energy
    f2.write(tmp[16])
    f2.write('\n')


    for iatom in range(natom):
        tmp  = f.readline()
        f2.write("%s" %( tmp ))

f.close
f2.close

