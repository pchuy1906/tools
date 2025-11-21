import numpy as np

import argparse
parser = argparse.ArgumentParser(description='compute error dlars')
# Arguments supported by the code.
parser.add_argument('--file_traj',  default='traj.txt')
parser.add_argument('--file_label', default='label.txt')
parser.add_argument('--file_b',     default='b.txt')

args    = parser.parse_args()
file_traj     = args.file_traj
file_label    = args.file_label
file_b        = args.file_b


f  = open(file_traj,"rt")


b = np.loadtxt(file_b)
label = np.loadtxt(file_label,dtype=str)

id_forces = [i for i, s in enumerate(label) if "force"  in s]
id_energy = [i for i, s in enumerate(label) if "energy" in s]
id_stress = [i for i, s in enumerate(label) if "stress" in s]

ref_forces = b[id_forces]
ref_energy = b[id_energy]
ref_stress = b[id_stress]


print (label)

iterations = []
numvars = []
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    keywords = "Iteration"
    if keywords in line:
        i = line.split()[-1]
        iterations.append(int(i))

    keywords = "Number of vars"
    if keywords in line:
        i = line.split()[-1]
        numvars.append(int(i))


f2 = open("final_data.dat", "w")

for i in range(len(iterations)):
    if iterations[i]%100==1:
        file_x = f"x.{iterations[i]:05d}.txt"
        file_Ax = f"Ax.{iterations[i]:05d}.txt"
        ChIMES_values = np.loadtxt(file_Ax)
        ChIMES_forces = ChIMES_values[id_forces]
        ChIMES_energy = ChIMES_values[id_energy]
        ChIMES_stress = ChIMES_values[id_stress]
        RMSE_forces = np.sqrt(np.mean((ref_forces - ChIMES_forces)**2))
        RMSE_energy = np.sqrt(np.mean((ref_energy - ChIMES_energy)**2))
        RMSE_stress = np.sqrt(np.mean((ref_stress - ChIMES_stress)**2))
        f2.write("%d %15.9f %15.9f %15.9f %d\n" %( numvars[i], RMSE_forces, RMSE_energy, RMSE_stress, iterations[i] ))




#    for i in range(5):
#        tmp = f.readline().split()
#        #print (tmp)
#        if len(tmp) == 0:
#            exit()
#        x.append(float(tmp[0]))
#        y.append(float(tmp[1]))
#
#    x = np.array(x)
#    y = np.array(y)
#
#    #print (x, y)
#    x0 = np.mean(x)
#    dx = 0.5*(x[2]-x[1])
#    y0 = np.mean(y)
#    dy = 0.5*(y[3]-y[4])
#    print (x0, y0, dx, dy)


#print (iterations)
#print (numvars)
