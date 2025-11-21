import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans']})

import argparse
parser = argparse.ArgumentParser(description='plot contour from a 2d array')
# Arguments supported by the code.
parser.add_argument("--file_input",              default='file.dat', help='file input 2d array')
parser.add_argument("--contour_min", type=float, default=0.0,        help='min contour value')
parser.add_argument("--contour_max", type=float, default=8000.0,     help='max contour value')
parser.add_argument("--fscale",      type=float, default=1.0,        help='scaling factor for the quantity')
parser.add_argument("--vcmap",                   default='cm_jet',   help='cmap for contour plot')

args        = parser.parse_args()
file_input   = args.file_input
contour_min  = args.contour_min
contour_max  = args.contour_max
fscale       = args.fscale
vcmap        = args.vcmap


def print_guide():
    guide = """
####################################################################################
old name: ~/tools/others/plot_contour_plt_SiO2.py 

if we want to read maximum values of the contour from file, set "--contour_max -10"
####################################################################################
"""
    print (guide)

print_guide()

print ("")
print ("read file input:", file_input)
print ("")

f  = open(file_input ,"r")
tmp = f.readline()
tmp = tmp.split()
Ny = int(tmp[0])
Nz = int(tmp[1])
print (Ny, Nz)

y = np.zeros(shape=(Ny))
z = np.zeros(shape=(Nz))
fyz = np.zeros(shape=(Ny,Nz))

for iy in range(Ny):
    for iz in range(Nz):
        tmp = f.readline()
        tmp = tmp.split()
        y[iy] = float(tmp[0])
        z[iz] = float(tmp[1])
        fyz[iy,iz] = float(tmp[2]) * fscale
f.close
if (contour_max<0.01):
    contour_max = np.amax(fyz)

ave_z = np.zeros(shape=(Nz))
collect_z = []
collect_ave = []
for iz in range(Nz):
    #ave_z[iz] = sum(fyz[:,iz])/float(Ny)
    #ave_z[iz] = sum( fyz[:,iz] * y[:] )/sum(fyz[:,iz])
    if (sum(fyz[:,iz])>0.001):
        collect_z.append(z[iz])
        collect_ave.append(sum( fyz[:,iz] * y[:] )/sum(fyz[:,iz]))

res = np.vstack((collect_z, collect_ave)).T
newfile = "ave_"+file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f')

import matplotlib.pyplot as plt
#plt.style.use('seaborn-white')

Y, Z = np.meshgrid(z, y)

v  = np.linspace(contour_min, contour_max, 2000, endpoint=True)
v2 = np.linspace(contour_min, contour_max, 5,   endpoint=True)

from matplotlib import cm
if (vcmap=="cm_jet"):
    vcmap=cm.jet
plt.contourf(Y, Z, fyz, 80, cmap=vcmap, levels=v)
#plt.contourf(Y, Z, fyz, 80, cmap=cm.jet, levels=v)

plt.colorbar(shrink=0.75,ticks=v2)

plt.xlim(z[0], z[-1])
plt.ylim(y[0], y[-1])

#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tick_params( axis='x', which='both', bottom=False, top=False, labelbottom=False)
#plt.tick_params( axis='y', which='both', right=False, left=False, labelbottom=False)
#plt.axis('off')

#plt.xlabel('Lz',fontsize=24)
#plt.ylabel('Ly',fontsize=24)

#plt.axis('scaled')

file_input = file_input.split('.')
pre_out = file_input[0]

#plt.savefig(pre_out + '.pdf')
plt.savefig(pre_out + '.png')


