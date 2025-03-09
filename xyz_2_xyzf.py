import numpy as np
import re

import argparse
parser = argparse.ArgumentParser(description='export FORCES, ENERGIES, and STRESS of file xyzfes')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz           = args.file_xyz

print ("")
print ("read fileXYZ:", file_xyz)

istruc = 0
fout="new_"+file_xyz
f2 = open(fout, "w")

f  = open(file_xyz ,"rt")
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    natom = int(tmp)
    f2.write("%-d\n" %(natom))
    
    tmp = f.readline()
    tmp2= re.findall(r'"([^"]*)"', tmp)
    tmp3= np.array(tmp2[0].split())
    cell_9 = tmp3.astype(np.float)
    #print ("AAA")
    #print (cell_9)
    cell_9 = cell_9.reshape((3, 3))
    vec1 = cell_9[0,:]
    vec2 = cell_9[1,:]
    vec3 = cell_9[2,:]
    if (vec2[2]<0.0):
        print (vec2)
        vec2 = vec2 + vec3
        print (vec2)
    x = np.vstack((vec1,vec2,vec3))
    x = x.reshape(9)
    cell_9 = np.array(["%.6f" % w for w in x.reshape(x.size)])
    cell_9 = cell_9.reshape(x.shape)
    #print (cell_9)
    cell_9_str = " ".join(cell_9)
    #print (cell_9_str)
    f2.write("%-s\n" %(cell_9_str))

    for k in range(natom):
        tmp = f.readline().split()
        f2.write("%-s %s %s %s %s %s %s\n" %(tmp[0], tmp[1],tmp[2],tmp[3],  tmp[6],tmp[7],tmp[8] ))
    istruc += 1
f.close

print 
print ("the number of structures:", istruc)
print 


