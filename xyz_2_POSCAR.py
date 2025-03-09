#import os
import numpy as np
#from numpy import *
#from numpy import matrix
#from numpy import linalg
#import commands

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--option_cell", default=1, help='0-(box,xyz)   1-(box3)  2-(box9) 3-(box9:NON_ORTHO)')

args        = parser.parse_args()
file_xyz     = args.file_xyz
option_cell   = int(args.option_cell)

def Cell_XYZ_ABC(cellXYZ):
    a1= cellXYZ[0,:]
    a2= cellXYZ[1,:]
    a3= cellXYZ[2,:]
    a = math.sqrt(dot(a1, a1))
    b = math.sqrt(dot(a2, a2))
    c = math.sqrt(dot(a3, a3))
    alp = math.acos(dot(a2, a3)/(b*c))*180.0/pi
    bet = math.acos(dot(a1, a3)/(a*c))*180.0/pi
    gam = math.acos(dot(a1, a2)/(a*b))*180.0/pi
    return np.array([a,b,c,alp,bet,gam])

def Cell_ABC_XYZ(cellABC):
    tmp = np.zeros(shape=(3,3))
    tmp[0,0] = cellABC[0]
    tmp[1,0] = cellABC[1]*cos(cellABC[5])
    tmp[1,1] = cellABC[1]*sin(cellABC[5])
    tmp[2,0] = cellABC[2]*cos(cellABC[4])
    tmp[2,1] = cellABC[2]*cos(cellABC[3])*sin(cellABC[5])-((cellABC[2]*cos(cellABC[4])-cellABC[2]*cos(cellABC[3])*cos(cellABC[5]))/tan(cellABC[5]))
    tmp[2,2] = sqrt((cellABC[2])**2 -(tmp[2,0])**2 - (tmp[2,1])**2)
    return tmp

if (option_cell==0):
    print ("")
    print ("read 2 files: box and xyz")
    print ("")
    print ("read Cell Parameters: file <box>")
    box = np.loadtxt('box')
    print ("")
    print ("read Atomic Coordinates: file <xyz>")
    status, natom = commands.getstatusoutput(" wc xyz | awk '{print $1}' ")
    natom = int(natom)
    
    f  = open('xyz' ,"r")
    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
      tmp = f.readline()
      tmp = tmp.split()
      myList.append(tmp[0])
      xyz[k,0] =  float(tmp[1])
      xyz[k,1] =  float(tmp[2])
      xyz[k,2] =  float(tmp[3])
    f.close
elif (option_cell==1):
    print ("")
    print ("Orthohombic cell case")
    print ("")
    print ("read fileXYZ:", file_xyz)
    f  = open(file_xyz ,"r")
    
    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)
    
    box = np.zeros(shape=(3,3))
    tmp = f.readline()
    tmp = tmp.split()
    for k in range(3):
        box[k,k] = float(tmp[k])
    
    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        myList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
elif (option_cell==2) or (option_cell==3):
    print ("")
    print ("cell_9")
    print ("")
    print ("read fileXYZ:", file_xyz)
    f  = open(file_xyz ,"r")

    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)

    box = np.zeros(shape=(3,3))
    tmp = f.readline()
    tmp = tmp.split()
    tmp = np.array(tmp)

    #box = [float(x) for x in tmp]
    #box = np.array(box)
    #box = box.reshape((3, 3))
    if (option_cell==2):
        box[0,:] = [float(x) for x in tmp[0:3]]
        box[1,:] = [float(x) for x in tmp[3:6]]
        box[2,:] = [float(x) for x in tmp[6:9]]
    elif (option_cell==3):
        box[0,:] = [float(x) for x in tmp[1:4]]
        box[1,:] = [float(x) for x in tmp[4:7]]
        box[2,:] = [float(x) for x in tmp[7:10]]
    else:
        print ("Not implemented yet!")
        exit()

    print (box)

    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        myList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
else:
    print ("Not implemented yet!")
    exit()



f2 = open("POSCAR", "w")
f2.write("%1s\n" %( "COMMENT" ))
f2.write("%15.9f\n" %( 1.0 ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[0,0], box[0,1], box[0,2] ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[1,0], box[1,1], box[1,2] ))
f2.write("%20.15f %20.15f %20.15f\n" %( box[2,0], box[2,1], box[2,2] ))

syms, counts_syms = np.unique(myList, return_counts=True)

print (syms)
print (counts_syms)

f3 = open("ntype.dat", "w")
for item in syms:
    f2.write("%s %s" %(" ", item))
    f3.write("%s\n" %( item))

f2.write("\n")
for item in counts_syms:
    f2.write("%s %s" %(" ", item))
f2.write("\n")
f2.write("%1s\n" %( "Direct" ))


xyz = np.dot(xyz, np.linalg.inv(box))

for item in syms:
    for k in range(0,natom):
        if (item == myList[k]):
            f2.write("%20.15f %20.15f %20.15f\n" %( xyz[k,0], xyz[k,1], xyz[k,2] ))

f2.close()
