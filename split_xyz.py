import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file with format xyzf, xyzfe, xyzfes')
args    = parser.parse_args()
file_xyz          = args.file_xyz

f  = open(file_xyz ,"rt")
nconf = 1

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)

    fname = "conf_" + str(nconf) + ".xyz"
    f2 = open( fname, "w")
    f2.write("%-s" %(tmp))
    
    tmp = f.readline()
    f2.write("%-s" %(tmp))

    for i in range(natom):
        tmp = f.readline()
        f2.write("%-s" %(tmp))

    f2.close()
    nconf += 1



