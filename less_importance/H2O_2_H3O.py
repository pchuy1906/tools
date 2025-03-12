import numpy as np

data = np.loadtxt('xyz.last_water', usecols = (1,2,3))

xyz_O = data[2,:]
xyz_H1= data[0,:]
xyz_H2= data[1,:]

#print (xyz_O)
xyz_H3= 3*xyz_O-xyz_H1-xyz_H2
#print(xyz_H3)


f2 = open("H3O.xyz", "w")
f2.write("%d\n" % 4)
f2.write("Comment\n")

tmp_xyz = xyz_O
f2.write("O ")
for ixyz in range(3):
    f2.write("%15.9f" % tmp_xyz[ixyz])
f2.write("\n")

tmp_xyz = xyz_H1
f2.write("H ")
for ixyz in range(3):
    f2.write("%15.9f" % tmp_xyz[ixyz])
f2.write("\n")

tmp_xyz = xyz_H2
f2.write("H ")
for ixyz in range(3):
    f2.write("%15.9f" % tmp_xyz[ixyz])
f2.write("\n")

tmp_xyz = xyz_H3
f2.write("H ")
for ixyz in range(3):
    f2.write("%15.9f" % tmp_xyz[ixyz])
f2.write("\n")

