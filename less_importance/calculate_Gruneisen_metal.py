import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='calculate Hugoniot point')


# Arguments supported by the code.
parser.add_argument("--file1", default='data1.dat', help='file_input1')
parser.add_argument("--file2", default='data2.dat', help='file_input2')
parser.add_argument("--unit_P", default='bars',      help='unit for pressure')
parser.add_argument("--unit_V", default='ang^3',    help='unit for pressure')
parser.add_argument("--unit_E", default='eV', help='unit for pressure')

args        = parser.parse_args()

file1   = args.file1
file2   = args.file2
unit_P  = args.unit_P
unit_V  = args.unit_V
unit_E  = args.unit_E

V1, p1, Etot1 = np.loadtxt(file1)
print "State with higher temperature:"
print "Volume, Press, Energy"
print V1, p1, Etot1
print

V2, p2, Etot2 = np.loadtxt(file2)
print "State with lower temperature:"
print "Volume, Press, Energy"
print V2, p2, Etot2
print

#bars_2_atm = 0.986923
#if (unit_P == 'bars'):
#    p1 = p1 * bars_2_atm
#    p2 = p2 * bars_2_atm
#eV_2_kcalmol = 23.0609
#if (unit_E == 'eV'):
#    Etot1 = Etot1 * eV_2_kcalmol
#    Etot2 = Etot2 * eV_2_kcalmol
#atm_2_GPa = 0.000101325
#GPa_A3_to_kcalmol = 6.0221367 / 41.84
#atm_A3_to_kcalmol = atm_2_GPa * GPa_A3_to_kcalmol
## unit atm*A3/kcalmol
#gruneisen = V2*(p2-p1)/(Etot2-Etot1) * atm_A3_to_kcalmol

eV_2_GPaA3 = 160.21766208
bar2GPa = 0.0001
gruneisen = V2*(p2-p1)/(Etot2-Etot1) * bar2GPa / eV_2_GPaA3 


print "gruneisen ", gruneisen
print
K = 2.0/gruneisen + 1.0
print "K ", K
print



