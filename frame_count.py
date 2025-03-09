#! /usr/bin/python3
# Counts the number of frames in an xyz file

import sys
if ( len(sys.argv) != 2 ):
    print("Counts the number of frames in an xyz file") 
    print("Usage: frame_count.py <xyz file>")
    exit(1)
    
argfile = sys.argv[1]

f = open(argfile, 'r')
atom_count = [] ;
c_count = []
h_count = []
n_count = []
o_count = []
cluster_count = 0
for line in f:
    natoms = int(line)
    c_count.append(0)
    h_count.append(0)
    n_count.append(0)
    o_count.append(0)
    
    if natoms < 0:
        print("Error reading xyz file at: ", line)
    atom_count.append(natoms) 
    for j in range(natoms+1):
        line = f.readline()
        if "C" in line:
            c_count[cluster_count] += 1
        elif "H" in line:
            h_count[cluster_count] += 1
        elif "N" in line:
            n_count[cluster_count] += 1
        elif "O" in line:
            o_count[cluster_count] += 1
            
    cluster_count = cluster_count + 1

print("Number of frames = ", str(cluster_count))
print("Frame #atoms #C  #H  #N  #O")
for j in range(cluster_count) :
    print("%5d %4d C %4s H %4s N %4s O %4s" % (j+1, atom_count[j], c_count[j], h_count[j], n_count[j], o_count[j]))
