import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate specific heat')

# Arguments supported by the code.
parser.add_argument("--file_input", default='COMPARE.dat', help='file_input')
parser.add_argument("--unit_E",     default='kcal_mol',    help='unit of energy: kcal_mol/eV')
parser.add_argument("--unit_P",     default='atm',         help='unit of energy: atm/GPa')
args        = parser.parse_args()
file_input   = args.file_input
unit_E       = args.unit_E
unit_P       = args.unit_P


# Boltzman constant, unit kcal/mol-K
kB = 0.001985875

kcal_2_kJ = 4.184
mfu = 255.14
nmolecules = 256
mass = nmolecules * mfu
g_2_kg = 0.001

print ("\nRead input file %s" % file_input)
data_T = np.loadtxt(file_input, usecols=1)
data_E = np.loadtxt(file_input, usecols=2)
data_p = np.loadtxt(file_input, usecols=3)
data_V = np.loadtxt(file_input, usecols=4)



if (unit_E=="eV"): data_E*=23.0609


if (unit_P=="GPa"): data_p*=9869.2326672


# E (kcal/mol), p (atm), V(Ang3), T(K)
## 1eV = 160.21766208 * GPa*A3 = 160.21766208 * 9869.2326672 atm*A3
## 1eV = 23.0609 kcal/mol
## Finally, atm*A3 = 23.0609 / (160.21766208 * 9869.2326672)  kcal/mol
fPV = 0.0000145841954095408
data_H = data_E + data_p * data_V * fPV

ave_p = np.mean(data_p)
ave_T = np.mean(data_T)
ave_E = np.mean(data_E)
ave_E2= np.mean(data_E*data_E)
ave_H = np.mean(data_H)
ave_H2= np.mean(data_H*data_H)
ave_V = np.mean(data_V)
ave_V2= np.mean(data_V*data_V)

sig_E = np.sqrt(ave_E2-ave_E*ave_E)
sig_H = np.sqrt(ave_H2-ave_H*ave_H)
sig_V = np.sqrt(ave_V2-ave_V*ave_V)

print ("")
#print ("factor1 = %15.9f\n" %(kcal_2_kJ / (mass * g_2_kg )  ))
#print ("factor2 = %15.9f\n" %(kcal_2_kJ / (mass * g_2_kg ) *mfu ))

#Cp = (sig_E*sig_E)/(kB*ave_T*ave_T)  * kcal_2_kJ / (mass * g_2_kg )
#print ("method-1-unit J/g-K %15.9f" % Cp)
#print ("method-1-unit J/mole-K %15.9f" % (Cp*mfu))

Cp = (sig_H*sig_H)/(kB*ave_T*ave_T)  * kcal_2_kJ / (mass * g_2_kg )
print ("method-2-unit J/g-K %15.9f" % Cp)
print ("method-2-unit J/mole-K %15.9f" % (Cp*mfu))

#Cp = (sig_E + sig_V*ave_p * fPV) * (sig_E + sig_V*ave_p * fPV) /(kB*ave_T*ave_T)  * kcal_2_kJ / (mass * g_2_kg )
#print ("method-3-unit J/g-K %15.9f" % Cp)
#print ("method-3-unit J/mole-K %15.9f" % (Cp*mfu))

