import numpy as np

import argparse
parser = argparse.ArgumentParser(description='find the H(p) from the fit')
# Arguments supported by the code.
parser.add_argument("--fname1",             default='fit.dat', help='filename from the fit')
parser.add_argument("--fname2",             default='fit.dat', help='filename from the fit')
parser.add_argument("--pmin",  type=float, default=0.0,       help='pressure min')
parser.add_argument("--pmax",  type=float, default=100.0,     help='pressure max')
parser.add_argument("--ndata",  type=int,  default=101,       help='number of data points')

args   = parser.parse_args()
fname1 = args.fname1
fname2 = args.fname2
pmin  = args.pmin
pmax  = args.pmax
ndata = args.ndata

data1 = np.loadtxt(fname1)
P1 = data1[:,0]
H1 = data1[:,1]

data2 = np.loadtxt(fname2)
P2 = data2[:,0]
H2 = data2[:,1]

f2 = open("delta_enthalpy.dat", "w")
Pdata = np.linspace(pmin, pmax, ndata)
for Pvalue in Pdata:
    delta1 = np.absolute(P1-Pvalue)
    i1 = delta1.argmin()
    delta2 = np.absolute(P2-Pvalue)
    i2 = delta2.argmin()

    if np.abs(P1[i1] - P2[i2]) < 0.05:
        if (np.abs(P1[i1])-Pvalue < 0.05):
            f2.write("%15.9f %15.9f %15.9f %15.9f\n" %( Pvalue, H2[i2]-H1[i1], H1[i1], H2[i2] ))

f2.close()

