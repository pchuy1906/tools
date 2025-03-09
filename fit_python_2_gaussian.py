import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='read data and fit to the function')
# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='XY data')
parser.add_argument("--xcol", type=int, default=0, help='column id of x-data')
parser.add_argument("--ycol", type=int, default=1, help='column id of y-data')
parser.add_argument("--gauwmin", type=float, default=0.01, help='min gauss-width')
parser.add_argument("--gauwmax", type=float, default=100.0, help='max gauss-width')
parser.add_argument("--minxgau1", type=float, default=0.0, help='min-peak-1')
parser.add_argument("--minxgau2", type=float, default=1.0, help='min-peak-2')
parser.add_argument("--maxxgau1", type=float, default=0.0, help='max-peak-1')
parser.add_argument("--maxxgau2", type=float, default=1.0, help='max-peak-2')


args        = parser.parse_args()
file_input     = args.file_input
xcol           = args.xcol
ycol           = args.ycol
gauwmin        = args.gauwmin
gauwmax        = args.gauwmax
minxgau1       = args.minxgau1
minxgau2       = args.minxgau2
maxxgau1       = args.maxxgau1
maxxgau2       = args.maxxgau2


print ("")
print ("")
print ("read file_input:", file_input)
print ("")
print ("")

import os
data = np.loadtxt(file_input, usecols=(xcol,ycol))

x = data[:,0]
y = data[:,1]

# guest for the parameters, can be arbitrary
gau1_1 = 1.0
gau2_1 = 1.0

gau1_2 = 1.0
gau2_2 = 2.0

gau1_3 = 1.0
gau2_3 = 1.0

print ("guess xvalues of gaussian:")
print (gau1_2, gau2_2)


print ("")
print ("")

def func(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3):
    return gau1_1 * np.exp(-(x - gau1_2)**2 / (2 * gau1_3**2)) + \
           gau2_1 * np.exp(-(x - gau2_2)**2 / (2 * gau2_3**2))

param_bounds=([0.1,0.01,1.0e-5,   0.1,0.01,1.0e-5],[1000.0,10.0,100.0,   1000.0,10.0,100.0])
popt,pcov = curve_fit(func,x,y,p0=[gau1_1,gau1_2,gau1_3,gau2_1,gau2_2,gau2_3],maxfev=5000000, bounds=param_bounds )

#popt,pcov = curve_fit(func,x,y,p0=[gau1_1,gau1_2,gau1_3,gau2_1,gau2_2,gau2_3],maxfev=5000000 )


print ("FITING RESULT: ", popt )

newx = np.linspace(x[0],x[-1],num=500)
yres = func(newx, *popt)

res = np.vstack((newx, yres)).T
newfile = "fit_" + file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f')

xgaus = np.array([popt[1],popt[4]])
xgaus = np.sort(xgaus)
print (xgaus)
print ("xgau1 = %15.9f" % xgaus[0])
print ("xgau2 = %15.9f" % xgaus[1])

