import numpy as np

Ha2kcal = 627.503

import argparse
parser = argparse.ArgumentParser(description='read Gaussian MP2 calculation')

# Arguments supported by the code.
parser.add_argument("--file_gauss_output", default='file.out', help='file gaussian output')
parser.add_argument("--file_xyz_output", default='file.xyz', help='file xyz output')
parser.add_argument("--method", default='DFT', help='calculation method (DFT/MP2)')


args        = parser.parse_args()

file_gauss_output   = args.file_gauss_output
file_xyz_output     = args.file_xyz_output
method              = args.method


f  = open(file_gauss_output,"r")
f2 = open(file_xyz_output,"w")


still_read_xyz  = True
still_read_fxyz = True
found_energy    = False

def get_asym(anum):
    if (str(anum)=="1"):
        return "H"
    if (str(anum)=="6"):
        return "C"
    if (str(anum)=="7"):
        return "N"
    if (str(anum)=="8"):
        return "O"

nconfig = 0

while True:

    line = f.readline().split()

    if (len(line)>3) and (line[0] == 'Job') and (line[1] == 'cpu') and (line[2] == 'time:') : break

    if (len(line)==2) and (line[0] == 'Z-Matrix') and (line[1] == 'orientation:'):
        line = f.readline().split()
        line = f.readline().split()
        line = f.readline().split()
        line = f.readline().split()
        still_read_xyz = True
        x = []
        y = []
        z = []
        atomList = []
        while (still_read_xyz):
            line = f.readline().split()
            if (len(line)==6):
                #print (line)
                if (int(line[1])>0):
                    atomList.append(get_asym(line[1]))
                    x.append(float(line[3]))
                    y.append(float(line[4]))
                    z.append(float(line[5]))
            else:
                still_read_xyz = False

    if (len(line)>3) and (line[0] == 'Center') and (line[1] == 'Atomic') and (line[2] == 'Forces'):
        line = f.readline().split()
        line = f.readline().split()
        still_read_fxyz = True
        fx = []
        fy = []
        fz = []
        while (still_read_fxyz):
            line = f.readline().split()
            if (len(line)==5):
                #print (line)
                fx.append(float(line[2]))
                fy.append(float(line[3]))
                fz.append(float(line[4]))
            else:
                still_read_fxyz = False

    if (method=="DFT"):
       if (len(line)>3) and (line[0] == 'SCF') and (line[1] == 'Done:') and (not found_energy):
           energy = float(line[4]) * Ha2kcal
           found_energy = True

    if (method=="MP2"):
       if (len(line)>3) and (line[0] == 'E2') and (line[3] == 'EUMP2') and (not found_energy):
           energy = float(line[5].replace('D','E')) * Ha2kcal
           found_energy = True

    if found_energy and (not still_read_xyz) and (not still_read_fxyz):
        nconfig += 1
        found_energy    = False
        still_read_xyz  = True
        still_read_fxyz = True

        asym = np.array(atomList)
        asym_unique = np.unique(asym)
        asym_unique = np.sort(asym_unique)
        asym_list = ' '.join(asym_unique)

        f2.write("%d\n" %(len(atomList)))
        f2.write("%d %d %d %15.9f %d %s\n" %( 200, 200, 200, energy, len(asym_unique), asym_list ))
        for k in range(len(atomList)):
            f2.write("%s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %( atomList[k], x[k], y[k], z[k], fx[k], fy[k], fz[k] ))

f2.close()
