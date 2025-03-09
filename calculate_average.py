import os
import math
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='calculate average')


# Arguments supported by the code.
parser.add_argument("--file_input", default='external_pressure.dat', help='file_input')
parser.add_argument("--column_ID", type=int, default=1, help='column to calculate average, 1 for 2nd column')
parser.add_argument("--skiprows", type=int, default=0, help='skiprows')
args      = parser.parse_args()
file_input   = args.file_input
column_ID    = args.column_ID
skiprows     = args.skiprows

F = np.loadtxt(file_input, usecols=[column_ID], skiprows=skiprows)

try:
    Ndata = F.shape[0]
    print ("shape of the file:")
    print (F.shape)
except:
    Ndata = 1

if Ndata > 1:
    ave = sum(F)/float(Ndata)
else:
    ave = F
print ("average= %20.4f" % ave)

if Ndata > 1:
    std = np.std(F)
else:
    std = 0
print ("error_bar=%20.4f" % std)

#res = np.vstack((ave, std)).T
#newfile = "results_ave_std.dat"
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')

