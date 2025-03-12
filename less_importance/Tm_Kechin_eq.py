import argparse
parser = argparse.ArgumentParser(description='fit Kechin equation for melting temperature as a function of pressure')
# Arguments supported by the code.
parser.add_argument("--fname", default='Tm.dat', help='melting temperature as a function of pressure')

args        = parser.parse_args()
fname      = args.fname

import numpy as np
import sys, math
from scipy.optimize import leastsq

# Kechin equation
def kechin(params, press):
    'From '
    T0, a, b, c = params 
    Tm = T0 * np.power(1.0+press/a,b) * np.exp(-c*press)
    return Tm

print "Filename containing Tm vs pressure:    ", fname
f = open(fname, 'rt')

print "Data read from file:"
press = []
Tm = []
while True:
    line = f.readline().strip()
    if line == '': break
    #if line[0] == '#' or line[0] == '!': continue
    p, t = [float(x) for x in line.split()[:2]]
    press.append(p)
    Tm.append(t)
print
f.close()

# transform to np arrays
press = np.array(press)
Tm = np.array(Tm)

## inital guess
T0 = 782.0
a = 3.5
b = 0.3
c = 0.0016
x0 = [T0,a,b,c]

target = lambda params, y, x: y - kechin(params, x)
opt_params, ier = leastsq(target, x0, args=(Tm, press))
print (opt_params)
#
#
xmin = min(press)
xmax = max(press)
print "xmin, xmax =", xmin, xmax
npoints = 50000
vfit = np.linspace(xmin, xmax, npoints)


efit = kechin(opt_params,vfit)
res = np.vstack((vfit, efit)).T
newfile = "fit_eos_murnaghan.dat"
np.savetxt(newfile, res, fmt='%15.9f %15.9f')


