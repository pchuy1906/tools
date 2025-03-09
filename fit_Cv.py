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
para1 = 80.00
para2 = 1000.0

print ("")
print ("")

def func(x, A, the):
    xx = the/x
    return A * xx*xx* np.exp(xx) / ( (np.exp(xx)-1) * (np.exp(xx)-1))

popt,pcov = curve_fit(func, x, y, p0=[para1, para2] )

print ("FITING RESULT: ", popt )

newx = np.linspace(x[0],x[-1],num=500)
yres = func(newx, *popt)
res = np.vstack((newx, yres)).T

newfile = "fit_" + file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f')
#with open(newfile, 'r') as original: data = original.read()
#with open(newfile, 'w') as modified: modified.write("# " + str(popt[0]) + " " + str(popt[1]) + "\n" + data)

#plt.plot(x, y, 'b+:', label='data')
#plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
#plt.legend()
#plt.show()
