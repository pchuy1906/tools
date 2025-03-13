import numpy as np
import os
import time
#import math

start = time.time()

import argparse
parser = argparse.ArgumentParser(description='Generate the twister bilayer graphene')
# Arguments supported by the code.
parser.add_argument("--int_m", type=int, default=1, help='parameter m')
parser.add_argument("--int_n", type=int, default=1, help='parameter n')
parser.add_argument("--d_inter", type=float, default=3.4, help='interlayer distance')
parser.add_argument("--d_intra", type=float, default=1.42, help='intralayer distance')
parser.add_argument("--types_layer1", default="CC", help='layer1 atom types')
parser.add_argument("--types_layer2", default="CC", help='layer2 atom types')

args        = parser.parse_args()
int_m   = args.int_m
int_n   = args.int_n
d_intra = args.d_intra
d_inter = args.d_inter
types_layer1 = args.types_layer1
types_layer2 = args.types_layer2


def Cell_XYZ_ABC(boxXYZ):
    a1= boxXYZ[0,:]
    a2= boxXYZ[1,:]
    a3= boxXYZ[2,:]
    a = np.sqrt(np.dot(a1, a1))
    b = np.sqrt(np.dot(a2, a2))
    c = np.sqrt(np.dot(a3, a3))
    alp = np.arccos(np.dot(a2, a3)/(b*c))
    bet = np.arccos(np.dot(a1, a3)/(a*c))
    gam = np.arccos(np.dot(a1, a2)/(a*b))
    return np.array([a,b,c,alp,bet,gam])

def Cell_ABC_XYZ(boxABC):
    tmp = np.zeros(shape=(3,3))
    tmp[0,0] = boxABC[0]
    tmp[1,0] = boxABC[1]*np.cos(boxABC[5])
    tmp[1,1] = boxABC[1]*np.sin(boxABC[5])
    tmp[2,0] = boxABC[2]*np.cos(boxABC[4])
    tmp[2,1] = boxABC[2]*np.cos(boxABC[3])*np.sin(boxABC[5])-((boxABC[2]*np.cos(boxABC[4])-boxABC[2]*np.cos(boxABC[3])*np.cos(boxABC[5]))/np.tan(boxABC[5]))
    tmp[2,2] = np.sqrt((boxABC[2])**2 -(tmp[2,0])**2 - (tmp[2,1])**2)
    return tmp

def unit_cell_t(m, n):
    if (np.gcd(n,3)==1):
        t1x = m
        t1y = m+n
        t2x = -m-n
        t2y = n+2*m
        natom = 4*( (n+m)*(n+m)+m*(n+2*m) )
        tan_ang = float(n)/( float(n+2*m) *np.sqrt(3.0) )
    if (np.gcd(n,3)==3):
        t1x = n/3+m
        t1y = n/3
        t2x = -n/3
        t2y = 2*n/3+m
        natom = 4*( m*m+m*n+n*n/3 )
        tan_ang = -np.sqrt(3.0)*float(m)/float(2*n+3*m)
    return natom, (t1x, t1y), (t2x, t2y), np.arctan(tan_ang)

def vector_a_b(phi, phi_TL):
    a1x =  np.sqrt(3.0)/2.0*np.cos(phi) -3.0/2.0*np.sin(phi)
    a1y = -np.sqrt(3.0)/2.0*np.sin(phi) -3.0/2.0*np.cos(phi)
    a2x =  np.sqrt(3.0)*np.cos(phi)
    a2y = -np.sqrt(3.0)*np.sin(phi)
    a1  = np.array([a1x,a1y,0.0])
    a2  = np.array([a2x,a2y,0.0])
    b1  = (np.cos(phi_TL)-np.sin(phi_TL)/np.sqrt(3.0))*a1 + 2.0*np.sin(phi_TL)/np.sqrt(3.0)*a2
    b2  = -2.0*np.sin(phi_TL)/np.sqrt(3.0)*a1 + (np.cos(phi_TL)+np.sin(phi_TL)/np.sqrt(3.0))*a2
    return np.array([a1,a2]), np.array([b1,b2])

def twisted_angle(m, n):
    cos_phi = ( 0.5*n*n + 3.0*m*n + 3.0*m*m )/ ( n*n + 3.0*m*n + 3.0*m*m )
    return np.arccos(cos_phi)

def r_from_a(cell_a):
    r3 = -1.0/3.0*(cell_a[0,:]+cell_a[1,:])
    r1 = r3 + cell_a[0,:]
    r2 = r3 + cell_a[1,:]
    return r1, r2, r3

def xyz_2_atom_from_r(r1, r2, r3, zcoor):
    xyz1 =  r2
    xyz2 = -r3
    xyz1[2] = zcoor
    xyz2[2] = zcoor
    return np.array( [xyz1,xyz2] )

def write_xyz(file_xyz, xyz_2_atoms, type_2_atoms, cell_a, nx_minmax, ny_minmax, cell_xyz):
    cell_abc = Cell_XYZ_ABC(cell_xyz)
    gcell_xyz = Cell_ABC_XYZ(cell_abc)
    f3 = open("cell.dat", "w")
    f3.write("%15.9f %15.9f %15.9f\n" %( gcell_xyz[0,0], gcell_xyz[0,1], gcell_xyz[0,2] ))
    f3.write("%15.9f %15.9f %15.9f\n" %( gcell_xyz[1,0], gcell_xyz[1,1], gcell_xyz[1,2] ))
    f3.write("%15.9f %15.9f %15.9f\n" %( gcell_xyz[2,0], gcell_xyz[2,1], gcell_xyz[2,2] ))
    f3.close()
    f2 = open(file_xyz, "w")
    for ix in range(nx_minmax[0],nx_minmax[1]):
        for iy in range(ny_minmax[0],ny_minmax[1]):
            tran = float(ix)*cell_a[0,:] + float(iy)*cell_a[1,:]
            R = tran + xyz_2_atoms[0,:]
            R = np.dot(R, np.linalg.inv(cell_xyz))
            for ixyz in range(3):
                R[ixyz] = R[ixyz]-np.floor(R[ixyz])
            f2.write("%4s %15.9f %15.9f %15.9f\n" %( type_2_atoms[0], R[0], R[1], R[2] ))
            R = tran + xyz_2_atoms[1,:]
            R = np.dot(R, np.linalg.inv(cell_xyz))
            for ixyz in range(3):
                R[ixyz] = R[ixyz]-np.floor(R[ixyz])
            f2.write("%4s %15.9f %15.9f %15.9f\n" %( type_2_atoms[1], R[0], R[1], R[2] ))

    f2.close()

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
def writeDupp(file_xyz, fileDuplicate):
    cmd = "rm -rf " + fileDuplicate
    os.system(cmd)
    cmd = "awk '{print $2, $3, $4}' " + file_xyz + " > tmp"
    os.system(cmd)
    XYZ = np.loadtxt('tmp')
    os.system("rm -rf tmp")
    #XYZ = np.loadtxt('_tmp_xyz',usecols=[1,2,3])
    N0 = XYZ.shape[0]
    for k1 in range(0,N0-1):
        for k2 in range(k1+1,N0):
            test = checkDuplicate(XYZ[k1],XYZ[k2])
            if (test==1):
                cmd = "echo " +str(k2+1) + "  " + str(k1+1) + " >> " + fileDuplicate
                os.system(cmd)

def removeDupp(file_xyz, fileDuplicate, file_out_xyz):
    cmd = "awk '{print $1}' " + fileDuplicate + " > tmp1"
    os.system(cmd)
    tmp1 = np.loadtxt('tmp1')
    os.system("rm -rf tmp1")
    s = []
    for i in tmp1:
        if i not in s:
            s.append(i)
    cmd = "rm -rf " + file_out_xyz
    os.system(cmd)
    k = 0
    for line in open(file_xyz):
        line = line.rstrip('\n')
        k = k + 1
        count = 0
        for k2 in range(0,len(s)):
            if (k==s[k2]):
                count = count + 1
        if (count==0):
            cmd = "echo \"" + line + "\" >> " + file_out_xyz
            os.system(cmd)

a0 = d_intra
m = int_m
n = int_n
dz = 5
zcoor = d_inter

phi_TL = twisted_angle(m,n)

natom, t1, t2, phi = unit_cell_t(m, n)
print (t1)
print (t2)
print ("The number of atom is %d" %natom)
cell_a, cell_b = vector_a_b(phi, phi_TL)
cell_a = cell_a*a0
cell_b = cell_b*a0

r1,r2,r3 = r_from_a(cell_a)
xyz_2_atoms_layer1 = xyz_2_atom_from_r(r1, r2, r3, dz)
r1,r2,r3 = r_from_a(cell_b)
xyz_2_atoms_layer2 = xyz_2_atom_from_r(r1, r2, r3, dz+zcoor)


vec_t1 = t1[0]*cell_a[0,:] + t1[1]*cell_a[1,:]
vec_t2 = t2[0]*cell_a[0,:] + t2[1]*cell_a[1,:]
vec_t3 = [0.0, 0.0, 25.0]
cell_xyz = np.array( [vec_t1,vec_t2,vec_t3] )

#print (vec_t1)
#print (vec_t2)

tmp_min = min(t1[0], t2[0])
tmp_max = max(t1[0], t2[0])+1
if (tmp_min > 0): tmp_min = 0
nx_minmax = (tmp_min, tmp_max)

tmp_min = min(t1[1], t2[1])
tmp_max = max(t1[1], t2[1])+1
if (tmp_min > 0): tmp_min = 0
ny_minmax = (tmp_min, tmp_max)

file_xyz_layer1 = "_tmp_xyz_layer1"
#type_2_atoms_layer1 = ["B","N"]
#write_xyz(file_xyz_layer1, xyz_2_atoms_layer1, type_2_atoms_layer1, cell_a, nx_minmax, ny_minmax, cell_xyz)
write_xyz(file_xyz_layer1, xyz_2_atoms_layer1, types_layer1, cell_a, nx_minmax, ny_minmax, cell_xyz)

file_xyz_layer2 = "_tmp_xyz_layer2"
#type_2_atoms_layer2 = ["C","C"]
#write_xyz(file_xyz_layer2, xyz_2_atoms_layer2, type_2_atoms_layer2, cell_b, nx_minmax, ny_minmax, cell_xyz)
write_xyz(file_xyz_layer2, xyz_2_atoms_layer2, types_layer2, cell_b, nx_minmax, ny_minmax, cell_xyz)


end1 = time.time()
print(end1 - start)

fileDuplicate1 = "Duplicate_layer1.dat"
writeDupp(file_xyz_layer1, fileDuplicate1)
fileDuplicate2 = "Duplicate_layer2.dat"
writeDupp(file_xyz_layer2, fileDuplicate2)

end2 = time.time()
print(end2 - end1)

file_out_xyz_layer1 = "xyz.layer1"
removeDupp(file_xyz_layer1, fileDuplicate1, file_out_xyz_layer1)
file_out_xyz_layer2 = "xyz.layer2"
removeDupp(file_xyz_layer2, fileDuplicate2, file_out_xyz_layer2)

end3 = time.time()
print(end3 - end2)

