import heapq
import numpy as np

nmin = 10
file_rmin = "../1-fit/rmin.dat"
file_xyz  = "../0-data/subtract_AE_all.xyz"




rmins = np.loadtxt(file_rmin)

print (f"read file rmin {file_rmin}")
nrow, ncol = rmins.shape
print (f"the number of row and column is {nrow} {ncol}")

ilocs_min = []
for icol in range(ncol):
    rmin_col = rmins[:,icol]
    ilocs = heapq.nsmallest(nmin, range(len(rmin_col)), key=rmin_col.__getitem__)
    ilocs_min.extend(ilocs)

ilocs_min_unique = list(dict.fromkeys(ilocs_min))
ilocs_min_unique.sort(reverse=False)
print (ilocs_min_unique)




f  = open(file_xyz ,"rt")
f2 = open("rmins.xyzf", "w")
nconf = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    write_to_file = False
    if nconf in ilocs_min_unique:
        write_to_file = True
        f2.write("%-s" %(tmp))

    tmp = f.readline()
    if write_to_file:
        f2.write("%-s" %(tmp))

    for i in range(natom):
        tmp = f.readline()
        if write_to_file:
            f2.write("%-s" %(tmp))

    nconf += 1

