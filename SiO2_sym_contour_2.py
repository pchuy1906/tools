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

print ("symmetry, isym=", isym)

## first read file density,  
F = np.loadtxt(file_input_density)
X0 = F[:,0]
Y0 = F[:,1]
X0 = np.array(X0)
Y0 = np.array(Y0)
Ndata = F.shape[0]
## array Y0 has form: 0,0,0,[],0,0,0
## Y0 is symmetry with the center of []
indices = np.where(Y0 > 0)[0]
idmin = indices[0]
idmax = indices[-1]
idmiddle = (idmin+idmax)//2
xshift = 0.5*(X0[idmin]+X0[idmax])
print ("density profile is symmetry with the center z=", xshift)



f  = open(file_input ,"r")
tmp = f.readline().split()
Ny = int(tmp[0])
Nx = int(tmp[1])

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

def compute_2D_distribution(x, y, fyx, xlim, file_input):
    #print (x[0], x[-1], len(x))
    #print (y[0], y[-1], len(y))

    collect_ave = []

    for i in range(len(y)):
        a_tmp = 0.0
        for j in range(len(x)):
            if x[j]>= xlim[0] and x[j]<= xlim[-1]:
                a_tmp += fyx[i,j]
        #print (y[i], sum(fyx[i,:]))
        #collect_ave.append(sum(fyx[i,:]))
        collect_ave.append(a_tmp)

    res = np.vstack((y, collect_ave)).T
    newfile = "DISTRIBUTION_" + file_input
    np.savetxt(newfile, res, fmt='%15.9f %15.9f')


compute_2D_distribution(x, y, fyx, xlim, file_input)


plt.colorbar(shrink=0.75,ticks=v2)

plt.xlim(xlim[0], xlim[1] )
plt.ylim(y[0], y[-1])

plt.xlabel(r'$z$ $(\AA)$')
#plt.ylabel(r'cos $\theta_{\mathrm{OH}}$')
plt.ylabel(r'$q$')


file_input = file_input.split('.')
if isym==0:
    pre_out = file_input[0]
else:
    pre_out = "sym_" + file_input[0]

plt.savefig(pre_out + '.png', bbox_inches='tight')

