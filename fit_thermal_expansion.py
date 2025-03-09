import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='read data and fit to the function')
# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='XY data')

args        = parser.parse_args()
file_input     = args.file_input

print ("")
print ("")
print ("read file_input:", file_input)
print ("")
print ("")

import os
data = np.loadtxt(file_input)

x = data[:,0]
y = data[:,1]

# guest for the parameters, can be arbitrary
para1 = 1000.0
para2 = 10.0
para3 = 0.3

print ("")
print ("")

def Gauss(x, para1, para2, para3):
    return para1*x*x + para2*x+para3
def thermal_expansion(x, para1, para2, para3):
    return (2.*para1*x+para2)/(para1*x*x + para2*x+para3)

popt,pcov = curve_fit(Gauss, x, y, p0=[para1, para2, para3] )

print ("FITING RESULT: ", popt )

newx = np.linspace(x[0],x[-1],num=500)
yres = Gauss(newx, *popt)
qcal = thermal_expansion(newx, *popt)
res = np.vstack((newx, yres, qcal)).T

newfile = "fit_" + file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f %15.9f')

#with open(newfile, 'r') as original: data = original.read()
#with open(newfile, 'w') as modified: modified.write("# " + str(popt[0]) + " " + str(popt[1]) + "\n" + data)



#plt.plot(x, y, 'b+:', label='data')
#plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
#plt.legend()
#plt.show()
