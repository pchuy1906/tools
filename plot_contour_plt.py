import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

import argparse
parser = argparse.ArgumentParser(description='plot contour from a 2d array')
# Arguments supported by the code.
parser.add_argument("--file_input",              default='file.dat', help='file input 2d array')
parser.add_argument("--contour_min", type=float, default=0.0,        help='min contour value')
parser.add_argument("--contour_max", type=float, default=8000.0,     help='max contour value')
parser.add_argument("--plot_STIB",   type=int,   default=0,          help='plot STIBs or not')
parser.add_argument("--STIB_min",    type=float, default=0.0,        help='plot line for STIBs')
parser.add_argument("--STIB_max",    type=float, default=100.0,      help='plot line for STIBs')
parser.add_argument("--fscale",      type=float, default=1.0,        help='scaling factor for the quantity')
parser.add_argument("--vcmap",                   default='cm_jet',   help='cmap for contour plot')

args        = parser.parse_args()
file_input   = args.file_input
contour_min  = args.contour_min
contour_max  = args.contour_max
STIB_min     = args.STIB_min
STIB_max     = args.STIB_max
plot_STIB    = args.plot_STIB
fscale       = args.fscale
vcmap        = args.vcmap

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

ave_z = np.zeros(shape=(Nz))
for iz in range(Nz):
    ave_z[iz] = sum(fyz[:,iz])/float(Ny)
res = np.vstack((z, ave_z)).T
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

if (plot_STIB != 0):
    x_stib = y
    y1_stib = np.linspace(STIB_min, STIB_min, len(x_stib),   endpoint=True)
    y2_stib = np.linspace(STIB_max, STIB_max, len(x_stib),   endpoint=True)
    plt.plot(y1_stib, x_stib, 'r')
    plt.plot(y2_stib, x_stib, 'r')

plt.colorbar(shrink=0.75,ticks=v2)

plt.xlim(0, z[-1])
plt.ylim(0, y[-1])

#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tick_params( axis='x', which='both', bottom=False, top=False, labelbottom=False)
#plt.tick_params( axis='y', which='both', right=False, left=False, labelbottom=False)
plt.axis('off')

#plt.xlabel('Lz',fontsize=24)
#plt.ylabel('Ly',fontsize=24)

plt.axis('scaled')

file_input = file_input.split('.')
pre_out = file_input[0]

#plt.savefig(pre_out + '.pdf')
plt.savefig(pre_out + '.png')


