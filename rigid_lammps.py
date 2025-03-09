f4 = open("in.lammps", "w")
natom = 260
natom_per_mol = 26
nmol = natom / natom_per_mol

f4.write("%s\n" %("units	real"))
f4.write("%s\n" %("atom_style	atomic"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("pair_style   lj/cut 8.0"))
f4.write("%s\n" %("read_data	data.lammps"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("velocity 	all create 300.0 4928459"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("# unconnected bodies"))
f4.write("%s\n" %(""))

sym="fix 1 all rigid/nvt group " + str(nmol)
for k in range(nmol):
    id_begin = 1 + k*natom_per_mol
    id_end   = (1 + k)*natom_per_mol
    f4.write("%s %s %d %d\n" %("group	clump" + str(k+1)," id <> ", id_begin, id_end))
    sym = sym + " clump" + str(k+1)
sym = sym + " &"

f4.write("%s\n" %(""))
f4.write("%s\n" %(sym))
f4.write("%s\n" %("                      temp 300.0 300.0 5.0 reinit no"))

f4.write("%s\n" %(""))
for k in range(nmol):
    sym = "neigh_modify exclude group clump"+ str(k+1)+ " clump" + str(k+1)
    f4.write("%s\n" %(sym))

f4.write("%s\n" %(""))
f4.write("%s\n" %("thermo		100"))
f4.write("%s\n" %(""))
f4.write("%s\n" %("dump             dump_1 all custom 50 dump_xyz_vxyz id type x y z vx vy vz"))
f4.write("%s\n" %(""))
f4.write("%s\n" %(""))
f4.write("%s\n" %("timestep 	0.0001"))
f4.write("%s\n" %("thermo		50"))
f4.write("%s\n" %("run		10000"))

f4.close()

