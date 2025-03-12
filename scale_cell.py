# Dictionary of all elements matched with their atomic masses.
elements_dict_mass = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'Fe' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}




import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--option_cell", type=int,   default=1,          help='1-(box3:orthohombic cell)  2-(box9:general) 3-NONORTHO')
parser.add_argument("--f_scale",     type=float, default=1.0,        help='scale the cell parameters')
parser.add_argument("--ixyz",        type=int,   default=0,          help='0-xyz 1-x 2-y 3-z')

args        = parser.parse_args()
file_xyz     = args.file_xyz
option_cell  = args.option_cell
f_scale      = args.f_scale
ixyz         = args.ixyz

if (option_cell==1):
    print ("\n")
    print ("Orthohombic cell case")
    print ("read fileXYZ: %s\n" % file_xyz)
    f  = open(file_xyz ,"r")
    
    natom = f.readline()
    natom = int(natom)
    print ("the number of atom: %d\n" % natom)
    
    box = np.zeros(shape=(3,3))
    tmp = f.readline()
    tmp = tmp.split()
    for k in range(3):
        box[k,k] = float(tmp[k])
    
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
elif (option_cell==2):
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
    box = [float(x) for x in tmp]
    box = np.array(box)
    box = box.reshape((3, 3))

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
elif (option_cell==3):
    print ("")
    print ("NON_ORTHO")
    print ("")
    print ("read fileXYZ:", file_xyz)
    f  = open(file_xyz ,"r")

    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)

    box = np.zeros(shape=(3,3))
    tmp = f.readline()
    tmp = tmp.split()[1:]
    tmp = np.array(tmp)
    box = [float(x) for x in tmp]
    box = np.array(box)
    box = box.reshape((3, 3))

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
else:
    print ("Not implemented yet!")
    exit()


values, counts = np.unique(atomList, return_counts=True)
print ("system detail:")
mass = []
for i in range(len(values)):
    mass.append( elements_dict_mass.get(values[i]) )
    print ("%4s %10d %12.3f" % (values[i], counts[i],  mass[i] ) )

print ("\n")
Navo = 0.6022140857


mass_total = np.dot(counts, mass)
volume0 = abs(np.linalg.det(box))
rho0 = mass_total / (Navo*volume0)

print ("initial volume %f" % volume0)
print ("initial density %f" % rho0)

newbox = np.zeros(shape=(3,3))
for i in range(3):
    newbox[i,:] = box[i,:]

if ixyz==0:
    newbox[0,:] = box[0,:]*f_scale
    newbox[1,:] = box[1,:]*f_scale
    newbox[2,:] = box[2,:]*f_scale
elif ixyz==1:
    newbox[0,:] = box[0,:]*f_scale
elif ixyz==2:
    newbox[1,:] = box[1,:]*f_scale
elif ixyz==3:
    newbox[2,:] = box[2,:]*f_scale



volume_ = abs(np.linalg.det(newbox))
rho_ = mass_total/ (Navo*volume_)

print ("Final volume %f" % volume_)
print ("Final density %f" % rho_)

#print (xyz[0,:])
Cxyz = np.dot(xyz, np.linalg.inv(box))
#print (Cxyz[0,:])
Axyz = np.dot(Cxyz, newbox)
#print (Axyz[0,:])


f2 = open("scale_cell_"+str(f_scale)+".xyz", "w")
f2.write("%4d\n" %( natom))
if (option_cell==1):
    f2.write("%15.9f %15.9f %15.9f\n" %( newbox[0,0], newbox[1,1], newbox[2,2] ))
else:
    f2.write("%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %( newbox[0,0], newbox[0,1], newbox[0,2], newbox[1,0], newbox[1,1], newbox[1,2], newbox[2,0], newbox[2,1], newbox[2,2] ))

for k4 in range(0,natom):
    tmpp  = Axyz[k4,:]
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( atomList[k4], tmpp[0], tmpp[1], tmpp[2] ))
f2.close()
