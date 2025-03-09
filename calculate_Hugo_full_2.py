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

[V0,dV0], [T0,dT0], [p0,dp0], [Etot0,dEtot0] = np.loadtxt(file0)
print "Initial state:"
print "Volume, Temp, Press, Energy"
print V0, T0, p0, Etot0
print

[V1,dV1], [T1,dT1], [p1,dp1], [Etot1,dEtot1] = np.loadtxt(file1)
print "State with higher temperature:"
print "Volume, Temp, Press, Energy"
print V1, T1, p1, Etot1
print "Higher P(unit)/T(K):", p1, T1
print

[V2,dV2], [T2,dT2], [p2,dp2], [Etot2,dEtot2] = np.loadtxt(file2)
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
    dp0 = dp0 * bars_2_atm
    dp1 = dp1 * bars_2_atm
    dp2 = dp2 * bars_2_atm
elif (unit_P == 'GPa'):
    p0 = p0 * GPa_2_atm
    p1 = p1 * GPa_2_atm
    p2 = p2 * GPa_2_atm
    dp0 = dp0 * GPa_2_atm
    dp1 = dp1 * GPa_2_atm
    dp2 = dp2 * GPa_2_atm

eV_2_kcalmol = 23.0609
if (unit_E == 'eV'):
    Etot0 = Etot0 * eV_2_kcalmol
    Etot1 = Etot1 * eV_2_kcalmol
    Etot2 = Etot2 * eV_2_kcalmol
    dEtot0 = dEtot0 * eV_2_kcalmol
    dEtot1 = dEtot1 * eV_2_kcalmol
    dEtot2 = dEtot2 * eV_2_kcalmol


atm_2_GPa = 0.000101325
GPa_A3_to_kcalmol = 6.0221367 / 41.84
atm_A3_to_kcalmol = atm_2_GPa * GPa_A3_to_kcalmol

p0 = p0*atm_A3_to_kcalmol
p1 = p1*atm_A3_to_kcalmol
p2 = p2*atm_A3_to_kcalmol
dp0 = dp0*atm_A3_to_kcalmol
dp1 = dp1*atm_A3_to_kcalmol
dp2 = dp2*atm_A3_to_kcalmol


# Now the units: Energy kcal/mol, pressure kcal/mol-A3, volume A3, temperature K

dE1  = Etot1 - Etot0
dPV1 = 0.5 * (p1+p0) * (V0-V1)
The1 = dE1 - dPV1
print "Hugo1 -- should be positive"
print The1

dE2  = Etot2 - Etot0
dPV2 = 0.5 * (p2+p0) * (V0-V2)
The2 = dE2 - dPV2
print "Hugo2 -- should be negative"
print The2

if (The1*The2 < 0.0):
    print "Changed sign"

print
The_diff = The2 - The1
p_diff = p2 - p1
t_diff = T2 - T1
d_p_diff = dp2 + dp1
d_t_diff = dT2 + dT1

dpdh = p_diff/The_diff
dtdh = t_diff/The_diff

phug = p1 - The1*dpdh
thug = T1 - The1*dtdh
print "Hugoniot P(GPa)/T(K):", phug/atm_A3_to_kcalmol*atm_2_GPa, thug


#def dThe1_dE0(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -1.0
#def dThe1_dp0(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -0.5*(V0-V1)
#
#def dThe1_dE1(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return 1.0
#def dThe1_dp1(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -0.5*(V0-V1)
#
#def dThe2_dE0(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -1.0
#def dThe2_dp0(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -0.5*(V0-V2)
#
#def dThe2_dE2(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return 1.0
#def dThe2_dp2(Etot0,p0,V0,T0, Etot1,p1,V1,T1, Etot2,p2,V2,T2):
#    return -0.5*(V0-V2)

dThe1_dE0 = -1.0
dThe1_dp0 = -0.5*(V0-V1)
dThe1_dE1 = 1.0
dThe1_dp1 = -0.5*(V0-V1)
dThe1_dE2 = 0.0
dThe1_dp2 = 0.0

dThe2_dE0 = -1.0
dThe2_dp0 = -0.5*(V0-V2)
dThe2_dE1 = 0.0
dThe2_dp1 = 0.0
dThe2_dE2 = 1.0
dThe2_dp2 = -0.5*(V0-V2)

dPhug_dE0 = -(p2-p1)/((The2-The1)*(The2-The1)) * ((The2-The1)*dThe2_dE0-The2*(dThe2_dE0-dThe1_dE0))
dPhug_dE1 = -(p2-p1)/((The2-The1)*(The2-The1)) * ((The2-The1)*dThe2_dE1-The2*(dThe2_dE1-dThe1_dE1))
dPhug_dE2 = -(p2-p1)/((The2-The1)*(The2-The1)) * ((The2-The1)*dThe2_dE2-The2*(dThe2_dE2-dThe1_dE2))

dPhug_dp0 = -(p2-p1)/((The2-The1)*(The2-The1)) * ((The2-The1)*dThe2_dp0-The2*(dThe2_dp0-dThe1_dp0))
dPhug_dp1 = -(The2)/((The2-The1)*(The2-The1)) * (-1.0*(The2-The1)-(p2-p1)*(dThe2_dp1-dThe1_dp1))
dPhug_dp2 = ( (-The1+p1*dThe2_dp2)*(The2-The1)-(-p2*The1+p1*The2)*dThe2_dp2 )/((The2-The1)*(The2-The1))

dPhug = np.abs(dPhug_dE0)*dEtot0 + np.abs(dPhug_dE1)*dEtot1 + np.abs(dPhug_dE2)*dEtot2 + np.abs(dPhug_dp0)*dp0 + np.abs(dPhug_dp1)*dp1 + np.abs(dPhug_dp2)*dp2


#d_phug = dp1 + np.abs(hug1/hug_diff) * d_p_diff
#d_thug = dT1 + np.abs(hug1/hug_diff) * d_t_diff
print "std P(GPa)/T(K):", dPhug/atm_A3_to_kcalmol*atm_2_GPa

