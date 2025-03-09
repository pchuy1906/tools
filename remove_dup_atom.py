import os
import math
import numpy as np
import commands
from numpy import *
 
## read Cell parameters --> cell_ABC
def readCell(CIFfile):
    cmd = "grep _cell_length_a "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >  boxABC"
    os.system(cmd)
    cmd = "grep _cell_length_b "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >> boxABC"
    os.system(cmd)
    cmd = "grep _cell_length_c "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >> boxABC"
    os.system(cmd)
    cmd = "grep _cell_angle_alpha "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >> boxABC"
    os.system(cmd)
    cmd = "grep _cell_angle_beta  "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >> boxABC"
    os.system(cmd)
    cmd = "grep _cell_angle_gamma "+ CIFfile + " | awk '{print $2}' | sed -e 's/([^()]*)//g' >> boxABC"
    os.system(cmd)
    cellABC = loadtxt('boxABC')
    return cellABC

## convert cell_ABC --> cell_XYZ
def Cell_ABC_XYZ(cellABC):
    tmp_cellABC = cellABC
    tmp_cellABC[3:6] = tmp_cellABC[3:6]*pi/180.0
    tmp = np.zeros(shape=(3,3))
    tmp[0,0] = tmp_cellABC[0]
    tmp[1,0] = tmp_cellABC[1]*cos(tmp_cellABC[5])
    tmp[1,1] = tmp_cellABC[1]*sin(tmp_cellABC[5])
    tmp[2,0] = tmp_cellABC[2]*cos(tmp_cellABC[4])
    tmp[2,1] = tmp_cellABC[2]*cos(tmp_cellABC[3])*sin(tmp_cellABC[5])-((tmp_cellABC[2]*cos(tmp_cellABC[4])-tmp_cellABC[2]*cos(tmp_cellABC[3])*cos(tmp_cellABC[5]))/tan(tmp_cellABC[5]))
    tmp[2,2] = sqrt((tmp_cellABC[2])**2 -(tmp[2,0])**2 - (tmp[2,1])**2)
    os.system("rm -rf BOX")
    np.savetxt('BOX', tmp, fmt='%15.9f %15.9f %15.9f')
    return tmp  

## find in an array: smallest number, that larger than specified number
def goodnumber(specified_number,an_array):
    condition = 1
    count = 0
    while (condition==1):
      count = count + 1
      if (specified_number < an_array[count]):
        goodline = int(an_array[count])
        condition = 0
    return goodline

## read xyz from CIFfile
def readXYZ(CIFfile):
    cmd = "sed -n '/xyz_begin/=' " + CIFfile
    status, line_begin = commands.getstatusoutput(cmd)
    line_begin = int(line_begin)
    print "XYZ, line begins: ", line_begin
    cmd = "sed -n '/xyz_end/=' " + CIFfile
    status, line_end = commands.getstatusoutput(cmd)
    line_end = int(line_end)
    print "XYZ, line ends: ", line_end
    #cmd = "sed -n "+ str(line_begin+1) + "," + str(line_end-1)+"p "+ CIFfile +" | awk '{printf \"%6s %15.9f %15.9f %15.9f \\n\", $2, $3+0, $4+0, $5+0}'  > xyz.1"
    cmd = "sed -n "+ str(line_begin+1) + "," + str(line_end-1)+"p "+ CIFfile +" | awk '{printf \"%6s %15.9f %15.9f %15.9f \\n\", $1, $2+0, $3+0, $4+0}'  > xyz.1"
    print cmd
    os.system(cmd)

## read symmetry operators
def readSYM(CIFfile):
    cmd = "sed -n '/sym_begin/=' " + CIFfile
    status, line_begin = commands.getstatusoutput(cmd)
    line_begin = int(line_begin)
    print "symmetry, line begins: ", line_begin
    cmd = "sed -n '/sym_end/=' " + CIFfile
    status, line_end = commands.getstatusoutput(cmd)
    line_end = int(line_end)
    print "symmetry, line ends: ", line_end
    cmd = "sed -n "+ str(line_begin+1) + "," + str(line_end-1)+"p "+ CIFfile + " > operator"
    os.system(cmd)
    os.system("sed -i 's/x/$2/g' operator")
    os.system("sed -i 's/y/$3/g' operator")
    os.system("sed -i 's/z/$4/g' operator")
    os.system("sed -i 's/X/$2/g' operator")
    os.system("sed -i 's/Y/$3/g' operator")
    os.system("sed -i 's/Z/$4/g' operator")
    os.system("sed -i 's/,/,  /g' operator")
    ncolumn = np.loadtxt('operator', dtype='str').shape[1]
    print (ncolumn)
    if (ncolumn==4):
      os.system("awk '{print $2, $3, $4}' operator > tmp; mv tmp operator")
    cmd = "sed -i 's/'" + "\"'\"" + "'/ /g' operator"
    os.system(cmd)

## from symmetry operator and xyz, creat whole_xyz
def wholeXYZ():
    ## Input file: xyz.1 and operator
    os.system("rm -rf XYZ")
    for line in open('operator'):
      line = line.rstrip('\n')
      cmd = "awk '{printf \"%6s %15.9f %15.9f %15.9f \\n\",$1, " + line + "}' xyz.1 >> XYZ"
      os.system(cmd)

## check whether two point are the same or not
def checkDuplicate(xyz1,xyz2):
    res = 0 ## 2 points are different
    dxyz = np.subtract(xyz1, xyz2)
    nxyz = np.around(dxyz,decimals=0)
    tmp  = np.subtract(dxyz, nxyz)
    maxval = max(abs(i) for i in tmp)
    eps = 0.0001
    if (maxval < eps):
      res = 1 ## 2 points are the same
    return res

## write the Duplicate points and write in the file
def writeDupp():
    os.system("rm -rf fileDuplicate")
    os.system("awk '{print $2, $3, $4}' XYZ > tmp")
    XYZ = loadtxt('tmp')
    os.system("rm -rf tmp")
    N0 = XYZ.shape[0]
    for k1 in range(0,N0-1):
      for k2 in range(k1+1,N0):
        test = checkDuplicate(XYZ[k1],XYZ[k2])
        if (test==1):
          cmd = "echo " +str(k2+1) + "  " + str(k1+1) + " >> fileDuplicate"
          os.system(cmd)

def removeDupp():
    os.system("awk '{print $1}' fileDuplicate > tmp1")
    tmp1 = loadtxt('tmp1')
    os.system("rm -rf tmp1")
    #print tmp1
    s = []
    for i in tmp1:
      if i not in s:
        s.append(i)
    os.system("rm -rf nXYZ")
    k = 0
    for line in open('XYZ'):
      tmp = line.strip()
      if len(tmp) == 1: break
      #print len(tmp)
      #print (len(line))
      line = line.rstrip('\n')
      #print  len(line)
      k = k + 1
      count = 0
      for k2 in range(0,len(s)):
        if (k==s[k2]):
           count = count + 1
      if (count==0):
         cmd = "echo \"" + line + "\" >> nXYZ"
         os.system(cmd)

#####################################################################################################

import argparse

parser = argparse.ArgumentParser(description='remove duplicated atoms')

# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file.xyz')
parser.add_argument("--cell_option", default='cell_3', help='cell_3/cell_9')
args        = parser.parse_args()
file_xyz   = args.file_xyz
cell_option     = args.cell_option


print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

cell_3_3 = np.zeros(shape=(3,3))
tmp = f.readline().split()

if cell_option == "cell_3":
    for k in range(3):
        cell_3_3[k,k] = float(tmp[k])

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

Cxyz = np.dot(xyz, np.linalg.inv(cell_3_3))
print "Output: XYZ"
f2 = open('XYZ', "w")
for k in range(natom):
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], Cxyz[k,0], Cxyz[k,1], Cxyz[k,2] ))
f2.close



## Read cell from file and convert angles
####cellABC = np.zeros(shape=(6))
####cellABC = readCell(CIFfile)
####print "cellABC", cellABC
#####cellABC[3:6] = cellABC[3:6]*pi/180.0
####print
###### convert cell_ABC to cell_XYZ
####cellXYZ = np.zeros(shape=(3,3))
####cellXYZ = Cell_ABC_XYZ(cellABC)
####
####print "cellXYZ"
####print cellXYZ
####print
###### read xyz from CIFfile
####readXYZ(CIFfile)
####print "read XYZ ok"
####print
####
###### read symmetry operators
####readSYM(CIFfile)
####print "read symmetry ok"
####print 
###### create wholeXYZ
####wholeXYZ()
####
###### write the duplicated into file
writeDupp()
####
###### remove duplicated points
removeDupp()
####
status, natom = commands.getstatusoutput(" wc nXYZ | awk '{print $1}' ")
natom = int(natom)
####
f  = open('nXYZ' ,"r")
myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
  tmp = f.readline()
  tmp = tmp.split()
  myList.append(tmp[0])
  xyz[k,0] =  float(tmp[1])
  xyz[k,1] =  float(tmp[2])
  xyz[k,2] =  float(tmp[3])
f.close()
####
numrows = xyz.shape[0]
####numcols = xyz.shape[1]
####print "size of XYZ:", numrows, numcols
####
#####print "C2goodC calculation..."
#####for k1 in range(0,numrows):
#####  for k2 in range(0,numcols):
#####    tmp = xyz[k1][k2]-floor(xyz[k1][k2])
#####    xyz[k1][k2] = tmp
####
####f2 = open('XYZ_CRYSTAL', "w")
####for k in range(0,natom):
####    f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
####f2.close()
###### finalize
####os.system("cat tmp_in BOX > qe_CRYSTAL.in")
####os.system("echo ATOMIC_POSITIONS {crystal} >> qe_CRYSTAL.in")
####os.system("cat XYZ_CRYSTAL >> qe_CRYSTAL.in")
####cmd = "sed -i 's/nat=tmp/nat = "+str(int(natom))+"/g' qe_CRYSTAL.in"
####os.system(cmd)
####
####
cellXYZ = cell_3_3
xyz = np.dot(xyz, cellXYZ)
f4 = open('cell_111.xyz', "w")
f4.write("%d\n" %( numrows ))


cellXYZ = cellXYZ.reshape(9)

print (cellXYZ)
for icellxyz in cellXYZ:
    f4.write("  %15.9f" % icellxyz)
f4.write("\n")

for k in range(0,natom):
    f4.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
f4.close()
###### finalize
####os.system("cat tmp_in BOX > qe_ANGSTROM.in")
####os.system("echo ATOMIC_POSITIONS {angstrom} >> qe_ANGSTROM.in")
####os.system("cat XYZ_ANGSTROM >> qe_ANGSTROM.in")
####cmd = "sed -i 's/nat=tmp/nat = "+str(int(natom))+"/g' qe_ANGSTROM.in"
####os.system(cmd)
####
os.system("rm -rf nXYZ Check_by_hand tmp1 XYZ ")
####os.system("mv XYZ_ANGSTROM XYZ_CRYSTAL BOX fileDuplicate boxABC operator XYZ xyz.1 Check_by_hand/")
