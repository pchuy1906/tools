import numpy as np

import math 

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--_lambda",       type=float,   default=1.541838, help='_lambda')
args        = parser.parse_args()
_lambda = args._lambda

pi = math.pi

print ("")
print ("Read file data")
data = np.loadtxt('AAA.dat')

twoThe = data[:,1]*pi/180.0
intensity = data[:,2]
intensity = intensity/np.max(intensity)

q = (4*pi/_lambda) * np.sin(twoThe/2.0)

print (len(q))
combined_array = np.column_stack((q,intensity))
np.savetxt('BBB.dat', combined_array, delimiter=' ')


