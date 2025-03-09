import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='read data and fit to the function')
# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='XY data')
parser.add_argument("--xcol", type=int, default=0, help='column id of x-data')
parser.add_argument("--ycol", type=int, default=1, help='column id of y-data')

args        = parser.parse_args()
file_input     = args.file_input
xcol           = args.xcol
ycol           = args.ycol

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
gau1_1 = 100.0
gau2_1 = 75.0

gau1_2 = 1.3
gau2_2 = 1.75

gau1_3 = 10.0
gau2_3 = 5.0

print ("guess xvalues of gaussian:")
print (gau1_2, gau2_2)


print ("")
print ("")

def func(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3):
    return gau1_1 * np.exp(-(x - gau1_2)**2 / (2 * gau1_3**2)) + \
           gau2_1 * np.exp(-(x - gau2_2)**2 / (2 * gau2_3**2))

#def func(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3):
#    return gau1_1 / ((x - gau1_2)**2 + gau1_3) + \
#           gau2_1 / ((x - gau2_2)**2 + gau2_3)

#param_bounds=([0.1,0.01,1.0e-5,   0.1,0.01,1.0e-5],[1000.0,10.0,100.0,   1000.0,10.0,100.0])
#popt,pcov = curve_fit(func,x,y,p0=[gau1_1,gau1_2,gau1_3,gau2_1,gau2_2,gau2_3],maxfev=5000000, bounds=param_bounds )
popt,pcov = curve_fit(func,x,y,p0=[gau1_1,gau1_2,gau1_3,gau2_1,gau2_2,gau2_3] )


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

