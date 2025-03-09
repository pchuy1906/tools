import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file with format xyzf, xyzfe, xyzfes')
args    = parser.parse_args()
file_xyz          = args.file_xyz

id_rmin = np.loadtxt('id_rmin.dat')
#print (id_rmin)
id_rmin = np.array(id_rmin, dtype=int)
#print (id_rmin)
lenid = len(id_rmin)

f2 = open('configs_rmin.xyz', "w")


f  = open(file_xyz ,"rt")
iconf = 0
id_search_match = 0
id_match = id_rmin[id_search_match]
while iconf < id_rmin[-1]+1:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    
    tmp2 = f.readline()
    if (iconf == id_match):
        f2.write("%-s" %(tmp))
        f2.write("%-s" %(tmp2))

    for i in range(natom):
        tmp3 = f.readline()
        if (iconf == id_match):
            f2.write("%-s" %(tmp3))

    if (iconf == id_match):
        id_search_match += 1
        print (id_search_match, id_match)
        if (id_search_match < lenid):
            id_match = id_rmin[id_search_match]

    iconf += 1


