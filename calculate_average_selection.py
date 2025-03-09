import os
import math
import numpy as np
import pandas as pd
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='calculate average')


# Arguments supported by the code.
parser.add_argument("--file_input", default='external_pressure.dat', help='file_input')
parser.add_argument("--column_ID", type=int, default=1, help='column to calculate average, 1 for 2nd column')
parser.add_argument("--Xmin", type=float, default=20.0,  help='X>Xmin')
parser.add_argument("--Xmax", type=float, default=180.0, help='X<Xmax')

args        = parser.parse_args()

file_input   = args.file_input
column_ID    = args.column_ID
Xmin         = args.Xmin
Xmax         = args.Xmax

print ("")
print ("Read input file ", file_input)
F1 = pd.read_csv(file_input, header=None, delimiter=r"\s+", names=["X","Y"])
#print (F1.head())

F2 = F1[ (F1["X"]>Xmin) & (F1["X"]<Xmax) ]
#print (F2.head())

F = F2.values
#print (F)
Ndata = F.shape[0]
#print ("shape of the file:")
#print (F.shape)


ave = sum(F[:,column_ID])/float(Ndata)
print "average=", ave

RMSE = 0.0
for k in range(Ndata):
    RMSE += (F[k,column_ID]-ave) * (F[k,column_ID]-ave)
RMSE = sqrt(RMSE/float(Ndata))
print "error_bar=", RMSE


f2 = open("result.dat", "w")
f2.write("%15.9f %15.9f\n" %( ave, RMSE ))
f2.close()

