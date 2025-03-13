import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='equally spaces file xyz')
# Arguments supported by the code.
parser.add_argument("--file_xyz",          default='file.xyz', help='file with format xyzf, xyzfe, xyzfes')
parser.add_argument("--nevery", type=int,  default=100,        help='export file every this number of step')
args    = parser.parse_args()
file_xyz          = args.file_xyz
nevery            = args.nevery

f  = open( file_xyz ,"rt")
f2 = open( "nevery.xyz", "w")

nconf = 0

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)

    if (nconf%nevery == 0): f2.write("%-s" %(tmp))
    
    tmp = f.readline()
    if (nconf%nevery == 0): f2.write("%-s" %(tmp))

    for i in range(natom):
        tmp = f.readline()
        if (nconf%nevery == 0): f2.write("%-s" %(tmp))
    nconf += 1


f2.close()
