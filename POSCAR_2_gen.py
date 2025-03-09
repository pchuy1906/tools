import os
import numpy as np
from numpy import *
import commands

f_scale=1.0

f  = open('POSCAR' ,"r")
tmp = f.readline()
tmp = f.readline()

unitcell = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
unitcell[0,:] = tmp
tmp = f.readline()
tmp = tmp.split()
unitcell[1,:] = tmp
tmp = f.readline()
tmp = tmp.split()
unitcell[2,:] = tmp

print ("unitcell=",unitcell)
unitcell = unitcell * f_scale

tmp = f.readline()
atomNameList = tmp.split()
print ("atomNameList=",atomNameList)
asym_list = ' '.join(atomNameList)

tmp = f.readline()
atomNumList = tmp.split()
atomNumList = map(int, atomNumList)
print ("atomNumList=",atomNumList)

natom = sum(atomNumList)
print ("natom=", natom)
 
xyztype = f.readline().split()[0]
print (xyztype)

xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    xyz[k,:] =  tmp[:3]
f.close()

#print (xyz)
if xyztype=="direct" or xyztype=="Direct":
    print ("BBB")
    Axyz = np.dot(xyz, unitcell)
if xyztype=="Cartesian":
    print ("AAA")
    Axyz = xyz

f2 = open('POSCAR2gen.gen', "w")
f2.write("%-d %4s\n" %(natom, "S"))
f2.write("%-s \n" %(asym_list))

#atomNumList
#atomNameList
ixyz = 0
for i in range(len(atomNameList)):
    for j in range(atomNumList[i]):
        f2.write("%-5d %5d %15.9f %15.9f %15.9f\n" %(ixyz+1, i+1, Axyz[ixyz,0], Axyz[ixyz,1], Axyz[ixyz,2] ))
        ixyz += 1

f2.write("%15.9f%15.9f%15.9f\n" % (0.0,0.0,0.0))
for ixyz in range(3):
    for icellxyz in unitcell[ixyz,:]:
        f2.write("%15.9f" % icellxyz)
    f2.write("\n")

f2.close()
