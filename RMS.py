import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate RMS')

# Arguments supported by the code.
parser.add_argument("--file_input",           default='COMPARE.dat', help='file_input')
parser.add_argument("--col1",       type=int, default=0,             help='column 1')
parser.add_argument("--col2",       type=int, default=1,             help='column 2')
parser.add_argument("--skiprows",   type=int, default=0,             help='skiprows')
args        = parser.parse_args()
file_input   = args.file_input
col1         = args.col1
col2         = args.col2
skiprows     = args.skiprows

print ("Read input file ", file_input)
F = np.loadtxt(file_input, usecols=[col1,col2], skiprows=skiprows)
print (F.shape, F.shape[0])

MAE = np.mean(np.abs(F[:,0]-F[:,1]))
RMS = np.sqrt(np.mean((F[:,0]-F[:,1])**2))
ave_F_DFT = np.mean(np.abs(F[:,0]))
max_diff = np.max( np.abs(F[:,0]-F[:,1]) )

print ("RMS= %12.2f" %RMS)
print ("RMS_r= %12.2f" %(RMS/ave_F_DFT))
print ("MAE= %12.2f" %MAE)
print ("MAE_r= %12.2f" %(MAE/ave_F_DFT))
print ("maximum difference = %12.2f" %max_diff)

