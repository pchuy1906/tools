import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='calculate the distribution of variables')
# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='XY data')
parser.add_argument("--nbins", type=int, default=2000, help='number of bins')
parser.add_argument("--col_use", type=int, default=0, help='column to calculate histogram')
parser.add_argument("--min_val", type=float, default=-2000.0, help='min value to calculate histogram')
parser.add_argument("--max_val", type=float, default= 2000.0, help='max value to calculate histogram')


args        = parser.parse_args()
file_input     = args.file_input
nbins          = args.nbins
col_use        = args.col_use
min_val        = args.min_val
max_val        = args.max_val


print ("")
print ("")
print ("")
print ("read file_input:", file_input)


data = loadtxt(file_input, usecols=[col_use])


bins_inp = np.linspace( min_val, max_val, num= nbins)
hist,bins = np.histogram(data, bins = bins_inp) 
res = np.vstack((bins[:-1], hist)).T

np.savetxt('hist.dat', res, fmt='%15.9f %15.9f')

from scipy.signal import argrelextrema

max_posi = argrelextrema(hist, np.greater)[0]
print (max_posi)
x_max = res[max_posi,0]
y_max = res[max_posi,1]

print (x_max)
print (y_max)

idmax = np.argmax(y_max)
density1 = x_max[idmax]
print (density1)
y_max[idmax] = -1

idmax = np.argmax(y_max)
density2 = x_max[idmax]
print (density2)

f3 = open("dens.dat", "w")
if density1 > density2:
    f3.write("%15.9f\n" %( density1))
    f3.write("%15.9f\n" %( density2))
else:
    f3.write("%15.9f\n" %( density2))
    f3.write("%15.9f\n" %( density1))
f3.close()
