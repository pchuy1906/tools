#import os
import numpy as np
#from numpy import *
#from numpy import matrix
#from numpy import linalg
#import commands

import argparse
parser = argparse.ArgumentParser(description='convert to good xyz')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz           = args.file_xyz

print ("-----------------------------------------------------------")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")
f2 = open("new_" + file_xyz, "w")

istruc=0
while True:
    tmp  = f.readline()
    f2.write("%1s" %( tmp ))
    line = tmp.strip()
    if line == '': break
    natom = int(tmp)
    
    tmp = f.readline()
    f2.write("%1s" %( tmp ))
    tmp_ = tmp.split()

    for k in range(0,natom):
        tmp = f.readline().split()
        xyz = [float(tmp[i]) for i in range(1,4)]
        for i in range(3):
            if (xyz[i]>150):
                xyz[i] += -200
        print (xyz)
        f2.write("%4s %15.9f %15.9f %15.9f\n" %( tmp[0], xyz[0], xyz[1], xyz[2] ))
    istruc += 1
f.close

print 
print ("the number of structures:", istruc)
print 


