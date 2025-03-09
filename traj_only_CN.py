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

fname = "CN_" + file_xyz
f2 = open( fname, "w")

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)

    f2.write("%-s\n" %("@natom@"))
    
    tmp = f.readline()
    f2.write("%-s" %(tmp))

    natom_true = 0

    for i in range(natom):
        tmp = f.readline()
        if ("C" in tmp) or ("N" in tmp):
            natom_true += 1
            f2.write("%-s" %(tmp))

    nconf += 1


print ("number of configuration", nconf)
print ("number of atom", natom_true)
