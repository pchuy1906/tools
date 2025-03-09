import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands
import math

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

import argparse
parser = argparse.ArgumentParser(description='plot contour from a 2d array')
# Arguments supported by the code.
parser.add_argument("--file_input", default='file.dat', help='file input 2d array')
parser.add_argument("--contour_min", type=float, default=0.0, help='min contour value')
parser.add_argument("--contour_max", type=float, default=10000.0, help='max contour value')
parser.add_argument("--vcmap", default='cm_jet', help='cmap for contour plot')
parser.add_argument("--Nx", type=int, default=11, help='number of point x-direction')
parser.add_argument("--Ny", type=int, default=11, help='number of point y-direction')
parser.add_argument("--xmin", type=float, default=0.0, help='x-min')
parser.add_argument("--xmax", type=float, default=10.0, help='x-max')
parser.add_argument("--ymin", type=float, default=0.0, help='y-min')
parser.add_argument("--ymax", type=float, default=10.0, help='y-max')


args        = parser.parse_args()
file_input   = args.file_input
contour_min  = args.contour_min
contour_max  = args.contour_max
vcmap        = args.vcmap
Nx           = args.Nx
Ny           = args.Ny
xmin         = args.xmin
xmax         = args.xmax
ymin         = args.ymin
ymax         = args.ymax

print ("")
print ("read file input:", file_input)
print ("")

Fdata = np.loadtxt(file_input)
ndata = Fdata.shape[0]
if (Nx*Ny != ndata):
    exit(1)

x = np.zeros(shape=(Nx))
for ix in range(Nx):
    x[ix] = xmin + (xmax-xmin)/float(Nx-1)*float(ix)
y = np.zeros(shape=(Ny))
for iy in range(Ny):
    y[iy] = ymin + (ymax-ymin)/float(Ny-1)*float(iy)

print (x)
print (y)
F = Fdata[:,1]
Fmin = min(F)
F = F-Fmin
for k in range(len(F)):
    if (F[k] >= contour_max): F[k] = contour_max

F = F.reshape((Nx,Ny))
print (F)

import matplotlib.pyplot as plt
plt.style.use('seaborn-white')

X, Y = np.meshgrid(y, x)

v  = np.linspace(contour_min, contour_max, 2000, endpoint=True)
v2 = np.linspace(contour_min, contour_max, 5,   endpoint=True)

from matplotlib import cm
if (vcmap=="cm_jet"):
    vcmap=cm.jet

plt.contourf(X, Y, F, 150, cmap=vcmap, levels=v)

cbar = plt.colorbar(shrink=0.75,ticks=v2)
cbar.ax.tick_params(labelsize=20) 

plt.gcf().subplots_adjust(left=0.15, bottom=0.15)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plt.axis('off')
plt.xlabel('$\gamma (\mathrm{^{\circ}})$',fontsize=24)
plt.ylabel('a ($\mathrm{\AA}$)',fontsize=24)

#plt.axis('scaled')

ax = plt.gca()
#ax.set_ylim(2.33, y[-1])

#v_xticks = np.arange(x[0],x[-1],2)
#v_yticks = np.arange(y[0],y[-1],0.1)
#ax.set_xticks(v_xticks)
#ax.tick_params(direction='out', length=6, width=2, colors='k', grid_color='k', grid_alpha=0.5)


#file_input = file_input.split('.')
#pre_out = file_input[0]
#plt.savefig(pre_out + '.pdf')
#plt.savefig(pre_out + '.png')

plt.savefig('output.pdf')
