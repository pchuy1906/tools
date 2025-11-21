import numpy as np

import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='calculate average')


# Arguments supported by the code.
parser.add_argument("--file_input_density", default='density.dat', help='file_input_density')
parser.add_argument("--file_input",         default='contour.dat', help='file_input_contour')
parser.add_argument("--isym", type=int,     default=0,             help='1')
args      = parser.parse_args()
file_input_density   = args.file_input_density
file_input           = args.file_input
isym                 = args.isym

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

idmiddle_real = 0.5 * float(idmin+idmax)
print (idmiddle_real)
idmiddle_int = int(idmiddle_real)

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

xshift = 0.5*(x[idmin]+x[idmax-1])
x = x[idmin:idmax] - xshift
fyx = fyx[:,idmin:idmax]

Ny = fyx.shape[0]
Nx = fyx.shape[1]
Nxhalf = int(Nx/2)
print (Ny, Nx)

if isym != 0:
    for iy in range(Ny):
        for ix in range(Nxhalf):
            #print (iy, ix, Ny-1-iy, Nx-1-ix)
            if (isym==-1):
                f1 = 0.5*(fyx[iy,ix]+fyx[Ny-1-iy,Nx-1-ix])
                fyx[iy,ix] = f1
                fyx[Ny-1-iy,Nx-1-ix] = f1
            if (isym==1):
                f1 = 0.5*(fyx[iy,ix]+fyx[iy,Nx-1-ix])
                fyx[iy,ix] = f1
                fyx[iy,Nx-1-ix] = f1

Y, X = np.meshgrid(x, y)

contour_max = np.amax(fyx)
contour_min = 0.0
v  = np.linspace(contour_min, contour_max, 2000, endpoint=True)
v2 = np.linspace(contour_min, contour_max, 5,   endpoint=True)

from matplotlib import cm
vcmap=cm.jet
plt.contourf(Y, X, fyx, 80, cmap=vcmap, levels=v)

plt.colorbar(shrink=0.75,ticks=v2)

dx = 1.0
plt.xlim(x[0]-dx, x[-1]+dx)
plt.ylim(y[0], y[-1])

file_input = file_input.split('.')
if isym==0:
    pre_out = file_input[0]
else:
    pre_out = "sym_" + file_input[0]

plt.savefig(pre_out + '.png')

