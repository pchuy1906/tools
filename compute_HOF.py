# some constants
################
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
parser.add_argument("--fileAE_Gaussian",                        default='AE_Gaussian.dat', help='AE computed from Gaussian code')
parser.add_argument("--fileHOF_atom_0K",                        default='',                help='')
parser.add_argument("--fileEnthalpy_atom_difference",           default='',                help='')
parser.add_argument("--file_freq",                              default='',                help='')
parser.add_argument("--file_numCHNO",                           default='',                help='')
parser.add_argument("--file_energy",                            default='',                help='')
parser.add_argument("--shifted_AE",                   type=int, default=0,                 help='')
parser.add_argument("--nfu",                          type=int, default=1,                 help='')


args        = parser.parse_args()
fileAE_Gaussian                = args.fileAE_Gaussian
fileHOF_atom_0K                = args.fileHOF_atom_0K
fileEnthalpy_atom_difference   = args.fileEnthalpy_atom_difference
file_freq                      = args.file_freq
file_numCHNO                   = args.file_numCHNO
file_energy                    = args.file_energy
shifted_AE                     = args.shifted_AE
nfu                            = args.nfu

AE_Gaussian = np.loadtxt( fileAE_Gaussian )
print (AE_Gaussian)

HOF_atom_0K = np.loadtxt( fileHOF_atom_0K )
print (HOF_atom_0K)

Enthalpy_atom_difference = np.loadtxt( fileEnthalpy_atom_difference )
print (Enthalpy_atom_difference)

freq = np.loadtxt( file_freq )
print (freq)

ZPE  = sum(freq) * hbar * 0.5

numCHNO = np.loadtxt( file_numCHNO )
print (numCHNO)

if (shifted_AE==1):
    AE_method = np.loadtxt("AE_method_database.dat")
    AE_ref    = np.loadtxt("AE_ref_database.dat")

Hvib = 4*kb*T + hbar* (sum(freq * (0.5 + 1.0/(np.exp(hbar*freq/(kb*T))-1.0)) ) )  - hbar*sum(freq)*0.5
print ("thermal correction to enthalpy %15.1f" %(Hvib+ZPE))

x = hbar*freq/(kb*T)
y = x/(np.exp(x)-1.0) - np.log(1.0-np.exp(-x))
S = sum(y)
print ("standard entropy %15.1f" %S)


energy = np.loadtxt( file_energy )
if (shifted_AE==1):
    energy += -sum(numCHNO*AE_method) + sum(numCHNO*AE_ref)
print (energy)

atomization_energy =  (sum(numCHNO*AE_Gaussian)- (energy + ZPE) )
enthalpy_at_0K = sum(numCHNO* HOF_atom_0K) - atomization_energy
enthalpy_at_298K = enthalpy_at_0K + Hvib - sum(numCHNO * Enthalpy_atom_difference)

ZPE                *= 1.0/float(nfu)
atomization_energy *= 1.0/float(nfu)
enthalpy_at_0K     *= 1.0/float(nfu)
enthalpy_at_298K   *= 1.0/float(nfu)


print (" zero-point energy (kcal/mol) is %15.1f" %ZPE)
print ("atomization energy (kcal/mol) is %15.1f" %atomization_energy)
print ("    enthalpy at 0K (kcal/mol) is %15.1f" %enthalpy_at_0K)
print ("  enthalpy at 298K (kcal/mol) is %15.1f" %enthalpy_at_298K)


kcal_2_kJ = 4.18400
ZPE                *= kcal_2_kJ
atomization_energy *= kcal_2_kJ
enthalpy_at_0K     *= kcal_2_kJ
enthalpy_at_298K   *= kcal_2_kJ
print (" zero-point energy (kJ/mol) is %15.1f" %ZPE)
print ("atomization energy (kJ/mol) is %15.1f" %atomization_energy)
print ("    enthalpy at 0K (kJ/mol) is %15.1f" %enthalpy_at_0K)
print ("  enthalpy at 298K (kJ/mol) is %15.1f" %enthalpy_at_298K)


