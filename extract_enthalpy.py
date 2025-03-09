import numpy as np

import argparse
parser = argparse.ArgumentParser(description='find the H(p) from the fit')
# Arguments supported by the code.
parser.add_argument("--fname",             default='fit.dat', help='filename from the fit')
parser.add_argument("--pmin",  type=float, default=0.0,       help='pressure min')
parser.add_argument("--pmax",  type=float, default=100.0,     help='pressure max')
parser.add_argument("--ndata",  type=int,  default=101,       help='number of data points')

args   = parser.parse_args()
fname = args.fname
pmin  = args.pmin
pmax  = args.pmax
ndata = args.ndata

print "Filename containing energy vs volume:    ", fname

print "Data read from file:"
data = np.loadtxt(fname)
P = data[:,2]
H = data[:,3]

f2 = open("enthalpy_"+fname, "w")
Pdata = np.linspace(pmin, pmax, ndata)
for Pvalue in Pdata:
    delta = np.absolute(P-Pvalue)
    i = delta.argmin()
    if (i != 0) and (i < len(P)-1):
        if (P[i]-Pvalue)*(P[i+1]-Pvalue) < 0.0:
            Hvalue = 0.5*(H[i]+H[i+1])
            f2.write("%15.9f %15.9f\n" %( Pvalue, Hvalue ))
        if (P[i]-Pvalue)*(P[i-1]-Pvalue) < 0.0:
            Hvalue = 0.5*(H[i]+H[i-1])
            f2.write("%15.9f %15.9f\n" %( Pvalue, Hvalue ))

f2.close()

