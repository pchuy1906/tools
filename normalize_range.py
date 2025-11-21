import os
import math
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='normalize a function in a range')


# Arguments supported by the code.
parser.add_argument("--file_input",          default='data.dat', help='file_input')
parser.add_argument('--xranges',    nargs='+', type=float)

args      = parser.parse_args()
file_input   = args.file_input
xranges     = args.xranges

F = np.loadtxt(file_input)
x = F[:,0]
y = F[:,1]

mask = (x > xranges[0]) & (x < xranges[1])
ynew = y[mask]

max_value = np.max(ynew)

y_scale = y*1.0/max_value

res = np.vstack((x, y_scale)).T
newfile = "scaled_" + file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f')

