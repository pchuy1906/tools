import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='calculate Hugoniot point')

# Arguments supported by the code.
parser.add_argument("--file0", default='data0.dat', help='file_input0')
parser.add_argument("--file1", default='data1.dat', help='file_input1')
parser.add_argument("--file2", default='data2.dat', help='file_input2')
parser.add_argument("--unit_P", default='atm',      help='unit for pressure')
parser.add_argument("--unit_V", default='Ang3',     help='unit for volume')
parser.add_argument("--unit_T", default='K',        help='unit for temperature')
parser.add_argument("--unit_E", default='kcal/mol', help='unit for energy')

args        = parser.parse_args()

file0   = args.file0
file1   = args.file1
file2   = args.file2
unit_P  = args.unit_P
unit_V  = args.unit_V
unit_T  = args.unit_T
unit_E  = args.unit_E

V0, T0, p0, Etot0 = np.loadtxt(file0)
print "Initial state:"
print "Volume, Temp, Press, Energy"
print V0, T0, p0, Etot0
print

V1, T1, p1, Etot1 = np.loadtxt(file1)
print "State with higher temperature:"
print "Volume, Temp, Press, Energy"
print V1, T1, p1, Etot1
print "Higher P(unit)/T(K):", p1, T1
print

V2, T2, p2, Etot2 = np.loadtxt(file2)
print "State with lower temperature:"
print "Volume, Temp, Press, Energy"
print V2, T2, p2, Etot2
print "Lower P(unit)/T(K):", p2, T2
print

bars_2_atm = 0.986923
GPa_2_atm  = 9869.23
if (unit_P == 'bars'):
    p0 = p0 * bars_2_atm
    p1 = p1 * bars_2_atm
    p2 = p2 * bars_2_atm
elif (unit_P == 'GPa'):
    p0 = p0 * GPa_2_atm
    p1 = p1 * GPa_2_atm
    p2 = p2 * GPa_2_atm

eV_2_kcalmol = 23.0609
if (unit_E == 'eV'):
    Etot0 = Etot0 * eV_2_kcalmol
    Etot1 = Etot1 * eV_2_kcalmol
    Etot2 = Etot2 * eV_2_kcalmol

atm_2_GPa = 0.000101325
GPa_A3_to_kcalmol = 6.0221367 / 41.84
atm_A3_to_kcalmol = atm_2_GPa * GPa_A3_to_kcalmol

dE1  = Etot1 - Etot0
dPV1 = 0.5 * (p1+p0) * (V0-V1) * atm_A3_to_kcalmol
hug1 = dE1 - dPV1
print "Hugo1 -- should be positive"
print hug1

dE2  = Etot2 - Etot0
dPV2 = 0.5 * (p2+p0) * (V0-V2) * atm_A3_to_kcalmol
hug2 = dE2 - dPV2
print "Hugo2 -- should be negative"
print hug2

if (hug1*hug2 < 0.0):
    print "Changed sign"

print
hug_diff = hug2 - hug1
p_diff = p2 - p1
t_diff = T2 - T1
dpdh = p_diff/hug_diff
dtdh = t_diff/hug_diff
phug = p1 - hug1*dpdh
print
thug = T1 - hug1*dtdh
print "Hugoniot P(GPa)/T(K):", phug*atm_2_GPa, thug
print


