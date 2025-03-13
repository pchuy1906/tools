import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='different formats for file xyzf')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--in_format", default='non_orthor-stress6-energy', help='input format')
parser.add_argument("--out_format", default='cell3-energy', help='output format')

args        = parser.parse_args()
file_xyz           = args.file_xyz
in_format          = args.in_format
out_format         = args.out_format

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
    if (in_format == "non_orthor-stress6-energy"):
        cell9   = [float(x) for x in tmp[1:10]]
        stress6 = [float(x) for x in tmp[10:16]]
        energy  = float(tmp[16])
        cell3   = [cell9[0],cell9[4],cell9[8]]
        stress3 = [stress6[0],stress6[1],stress6[2]]
    if (in_format == "cell3-stress6-energy"):
        cell3   = [float(x) for x in tmp[0:3]]
        stress6 = [float(x) for x in tmp[3:9]]
        energy  = float(tmp[9])

    if (out_format=="cell3-energy"):
        # print the cell length first
        f2.write(`cell3[0]`+ " " + `cell3[1]` + " " + `cell3[2]` + " ")
        # then the stress tensor
        #f2.write(tmp[3]+ " " + tmp[4] + " " + tmp[5] + " " + tmp[6] + " " + tmp[8] + " " + tmp[7] + " ")
        # then the energy
        f2.write(`energy`)
        f2.write('\n')

    for iatom in range(natom):
        tmp  = f.readline()
        f2.write("%s" %( tmp ))

f.close
f2.close

