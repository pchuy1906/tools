import numpy as np

import argparse
parser = argparse.ArgumentParser(description='export dump frame with selective density')
# Arguments supported by the code.
parser.add_argument("--file_AAA_density",                     default='AAA_Density.dat', help='file AAA_Density.dat')
parser.add_argument("--file_lmp_dump",                        default='dump_xyz_vxyz',   help='file dump_xyz_vxyz')
parser.add_argument("--target_density",   type=float,         default=1.0,               help='desired density')



args        = parser.parse_args()
file_AAA_density = args.file_AAA_density
file_lmp_dump    = args.file_lmp_dump
target_density   = args.target_density


print ("reading file_AAA_density ", file_AAA_density)
data = np.loadtxt(file_AAA_density)

#print (data)
arr_density = data[:,1]
min_dens = arr_density[-1]
max_dens = arr_density[0]

print (f"density should be in range {min_dens}, {max_dens}")
if min_dens > target_density or max_dens < target_density:
    print (f"target density {target_density} is not in the range [{min_dens}, {max_dens}]")
    exit()
else:
    # since the array is in descending order, find the first postion that is smaller than targer value
    indices = np.where(arr_density < target_density)
    find_id = indices[0][0]
    print (data[find_id,:])
    print (data[find_id-1,:])
    frame_num = int(data[find_id,0])
    f2 = open('density_'+str(target_density)+'.dump','w')

f = open(file_lmp_dump, 'rt')
while True:
    line = f.readline()
    if line == '': break

    keywords = "TIMESTEP"
    if keywords in line:
        line2 = f.readline()
        frame_id = int(line2.split()[0])
        if frame_id == frame_num:
            print (frame_id)

            f2.write("%s" %( line ))

            f2.write("%s" %( line2 ))

            line3 = f.readline()
            f2.write("%s" %( line3 ))

            line3 = f.readline()
            f2.write("%s" %( line3 ))
            natom = int(line3.split()[0])

            for i in range(5):
                line3 = f.readline()
                f2.write("%s" %( line3 ))

            for i in range(natom):
                line3 = f.readline()
                f2.write("%s" %( line3 ))

