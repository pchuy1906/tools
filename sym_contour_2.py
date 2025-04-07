import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

import argparse

parser = argparse.ArgumentParser(description='calculate average')


# Arguments supported by the code.
parser.add_argument("--file_input_density", default='density.dat', help='file_input_density')
parser.add_argument("--file_input",         default='contour.dat', help='file_input_contour')
parser.add_argument("--isym", type=int,     default=0,             help='1')
parser.add_argument('--xlim', nargs='+', type=float)

args      = parser.parse_args()
file_input_density   = args.file_input_density
file_input           = args.file_input
isym                 = args.isym
xlim                 = args.xlim

F = np.loadtxt(file_input_density)
X0 = F[:,0]
Y0 = F[:,1]
X0 = np.array(X0)
Y0 = np.array(Y0)

Ndata = F.shape[0]
print ("shape of the file:")
print (F.shape)

indices = np.where(Y0 > 0)[0]
idmin = indices[0]
idmax = indices[-1]
print (idmin, idmax)

idmiddle = (idmin+idmax)//2

f  = open(file_input ,"r")
tmp = f.readline().split()
Ny = int(tmp[0])
Nx = int(tmp[1])
print (Ny, Nx)

y = np.zeros(shape=(Ny))
x = np.zeros(shape=(Nx))
fyx = np.zeros(shape=(Ny,Nx))

for iy in range(Ny):
    for ix in range(Nx):
        tmp = f.readline().split()
        y[iy] = float(tmp[0])
        x[ix] = float(tmp[1])
        fyx[iy,ix] = float(tmp[2])
f.close

xshift = 0.5*(x[idmin]+x[idmax])
x = x - xshift

if isym != 0:
    for iy in range(Ny):
        for i in range(idmiddle-idmin+1):
            if (isym==-1):
                f1 = 0.5*(fyx[iy,idmin+i]+fyx[Ny-1-iy,idmax-i])
                fyx[iy,     idmin+i] = f1
                fyx[Ny-1-iy,idmax-i] = f1
            if (isym==1):
                f1 = 0.5*(fyx[iy,idmin+i]+fyx[iy,idmax-i])
                fyx[iy,idmin+i] = f1
                fyx[iy,idmax-i] = f1

Y, X = np.meshgrid(x, y)

contour_max = np.amax(fyx)
contour_min = 0.0
v  = np.linspace(contour_min, contour_max, 2000, endpoint=True)
v2 = np.linspace(contour_min, contour_max, 5,   endpoint=True)
v2 = np.linspace(0.0, 0.4, 5,   endpoint=True)

from matplotlib import cm
vcmap=cm.jet
plt.contourf(Y, X, fyx, 80, cmap=vcmap, levels=v)

plt.colorbar(shrink=0.75,ticks=v2)

plt.xlim(xlim[0], xlim[1] )
plt.ylim(y[0], y[-1])

plt.xlabel(r'$z$ $(\AA)$')
plt.ylabel(r'cos $\mathrm{\theta_{OH}}$')

file_input = file_input.split('.')
if isym==0:
    pre_out = file_input[0]
else:
    pre_out = "sym_" + file_input[0]

plt.savefig(pre_out + '.png', bbox_inches='tight')

