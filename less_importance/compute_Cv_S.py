# some constants
################
# gas constant (J/mol-K)
R = 8.31446261815324
# to (kcal/mol-K)
R = R* 0.239001 / 1000.0
# Boltzman (kcal/mol-K)
kb = 1.987204259e-3
# Temperature (K)
T = 298.00
# Planck constant (kcal/mol / (cm-1))
hbar =  0.00285911

import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

print ("-----------------------------------------------------------")
print ("Calculate the HOF")
print ("")

import argparse
parser = argparse.ArgumentParser(description='compute HOF for CHNO molecules')
# Arguments supported by the code.
parser.add_argument("--file_freq",                              default='',                help='')
parser.add_argument("--nfu",                          type=int, default=1,                 help='')

args        = parser.parse_args()
file_freq                      = args.file_freq
nfu                            = args.nfu

freq = np.loadtxt( file_freq )
print (freq)

#ZPE  = sum(freq) * hbar * 0.5
#
Tmin = 100.0
Tmax = 1000.0
npoint=101
dT = (Tmax-Tmin)/(float(npoint-1))
f2 = open("Cv_T.dat", "w")

for i in range(npoint):
    T0 = Tmin + dT*(float(i-1))
    x = hbar*freq/(kb*T0)
    y = kb*x*x*np.exp(x)/((np.exp(x)-1.0)*((np.exp(x)-1.0)))
    Cv= 0.5*sum(y)/R
    f2.write("%15.9f %15.9f\n" %( T0, Cv ))
f2.close()

def _sinh(x):
    return 0.5*( np.exp(x)-np.exp(-x) )

def _coth(x):
    return ( np.exp(x)+np.exp(-x) )/ ( np.exp(x)-np.exp(-x) )

x = hbar*freq/(2.0*kb*T)
y1 = hbar*freq/(2.0*T)*_coth(x)
y2 = kb*np.log(2.0*_sinh(x))
S = sum(y1)-sum(y2)
#
S                *= 1.0/float(nfu)
print ("standard entropy (kcal/mol-K) %15.6f" %S)

kcal_2_kJ = 4.18400
S                *= kcal_2_kJ
print ("standard entropy (kJ/mol-K) is %15.6f" %S)


