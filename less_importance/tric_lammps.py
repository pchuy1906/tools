import numpy as np

file = 'boxABC'
data = np.loadtxt(file)

a = data[0]
b = data[1]
c = data[2]
alp = data[3] * np.pi / 180.0 
bet = data[4] * np.pi / 180.0 
gam = data[5] * np.pi / 180.0 

lx = a
xy = b * np.cos(gam)
xz = c * np.cos(bet)
ly = np.sqrt(b*b - xy*xy)
yz = (b*c*np.cos(alp)-xy*xz)/ly
lz = np.sqrt(c*c-xz*xz-yz*yz)


print ("%15.4f %15.4f xlo xhi" %(0.0, lx))
print ("%15.4f %15.4f ylo yhi" %(0.0, ly))
print ("%15.4f %15.4f zlo zhi" %(0.0, lz))
print ("%15.4f %15.4f %15.4f xy xz yz" %(xy, xz, yz))




#unix(['echo 0.000   ' num2str(Lattice(1,1),9) ' xlo xhi >  box.lammps']);
#unix(['echo 0.000   ' num2str(Lattice(2,2),9) ' ylo yhi >> box.lammps']);
#unix(['echo 0.000   ' num2str(Lattice(3,3),9) ' zlo zhi >> box.lammps']);
#unix(['echo ' num2str(Lattice(2,1),9) ' ' num2str(Lattice(3,1),9) ' ' num2str(Lattice(3,2),9) ' xy xz yz >> box.lammps']);
#
