import os
import numpy as np
import commands

import argparse

parser = argparse.ArgumentParser(description='QEq')

# Arguments supported by the code.
parser.add_argument("--nframe", type=int, default=100, help='file_input')

args              = parser.parse_args()
nframe            = args.nframe

print ("")
print ("read file A.txt")
A = np.loadtxt('A.txt')

nline = np.shape(A)[0]
ncolumn = np.shape(A)[1]

print ("")
print ("the number of line, column", nline, ncolumn)

print ("")
print ("assuming the number of frame", nframe)


natom = (nline/nframe-3)/3
print ("")
print ("Found the number of atoms", natom)


Aextra = []
for k in range(nframe):
    for iatom in range(natom):
        Aextra.append(0.0)
        Aextra.append(0.0)
        Aextra.append(0.0)
    Aextra.append(0)
    Aextra.append(0)
    Aextra.append(0)


Aextra = np.asarray(Aextra)
Aextra = Aextra.reshape(len(Aextra),1)
print (type(Aextra), type (A))
print (np.shape(Aextra), np.shape(A))

Anew = np.concatenate((A,Aextra),1)

np.savetxt('newA.txt', Anew)

#  tmp = f.readline()
#  tmp = tmp.split()
#  myList.append(tmp[0])
#  xyz[k,0] =  float(tmp[1])
#  xyz[k,1] =  float(tmp[2])
#  xyz[k,2] =  float(tmp[3])
#f.close
#
#numrows = xyz.shape[0]
#numcols = xyz.shape[1]
#print ("size of xyz:", numrows, numcols)
#
#print ("")
#print "Angstrom to Crystal calculation..."
#Cxyz = np.dot(xyz, np.linalg.inv(box))
#
##for k1 in range(0,numrows):
##  for k2 in range(0,numcols):
##    tmp = Cxyz[k1][k2]-floor(Cxyz[k1][k2])
##    Cxyz[k1][k2] = tmp
#
#print ("")
#print "Output: <Cxyz.out>"
#f2 = open('Cxyz.out', "w")
#for k in range(0,natom):
#  f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], Cxyz[k,0], Cxyz[k,1], Cxyz[k,2] ))
#f2.close
