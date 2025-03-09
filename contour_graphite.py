import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

import argparse
parser = argparse.ArgumentParser(description='plot contour from a 2d array')
# Arguments supported by the code.
parser.add_argument("--file_input", default='file.dat', help='file input 2d array')
parser.add_argument("--contour_min", type=float, default=0.0, help='min contour value')
parser.add_argument("--contour_max", type=float, default=8000.0, help='max contour value')


args        = parser.parse_args()
file_input   = args.file_input
contour_min  = args.contour_min
contour_max  = args.contour_max

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
        fyz[iy,iz] = float(tmp[2])
f.close


vmin = fyz.min()
for iy in range(Ny):
    for iz in range(Nz):
        fyz[iy,iz] = (fyz[iy,iz]-vmin)
        if (fyz[iy,iz] > contour_max):
            fyz[iy,iz] = contour_max


import matplotlib.pyplot as plt

plt.style.use('seaborn-white')

Y, Z = np.meshgrid(z, y)

v  = np.linspace(contour_min, contour_max, 2000, endpoint=True)
v2 = np.linspace(contour_min, contour_max, 5,   endpoint=True)


plt.plot(3.356, 2.463, 'wx')

from matplotlib import cm
plt.contourf(Y, Z, fyz, 150, cmap=cm.jet, levels=v)
cbar = plt.colorbar(shrink=0.75,ticks=v2)
cbar.ax.tick_params(labelsize=20) 

plt.gcf().subplots_adjust(left=0.15)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plt.axis('off')
plt.xlabel('c/2 ($\mathrm{\AA}$)',fontsize=24)
plt.ylabel('a ($\mathrm{\AA}$)',fontsize=24)

plt.axis('scaled')

ax = plt.gca()
ax.set_ylim(y[0], y[-1])

v_xticks = np.arange(z[0],z[-1]+0.001,0.2)
v_yticks = np.arange(y[0],y[-1]+0.001,0.1)
ax.set_xticks(v_xticks)
ax.tick_params(direction='out', length=6, width=2, colors='k', grid_color='k', grid_alpha=0.5)

file_input = file_input.split('.')
pre_out = file_input[0]
plt.savefig(pre_out + '.eps')


