import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
import numpy as np

fs_2_ps = 0.001
import argparse
parser = argparse.ArgumentParser(description='read data and fit to the function')
# Arguments supported by the code.
parser.add_argument("--file_input",             default='data.dat', help='XY data')
parser.add_argument("--xcol",       type=int,   default=0,          help='column id of x-data')
parser.add_argument("--ycol",       type=int,   default=1,          help='column id of y-data')
parser.add_argument("--tstep_fs",   type=float, default=0.1,        help='timestep in femtosecond')


args        = parser.parse_args()
file_input     = args.file_input
xcol           = args.xcol
ycol           = args.ycol
tstep_fs       = args.tstep_fs

print ("")
print ("read file_input:", file_input)
print ("")

data = np.loadtxt(file_input, usecols=(xcol,ycol))

x = data[:,0]
y = data[:,1]

x = x-x[0]
print (x)

x = x*tstep_fs*fs_2_ps
print (x)

# guest for the parameters, can be arbitrary
para1 = 1
para2 = 1
para3 = 1

def func(x, para1, para2, para3):
    return -para1 * np.exp(-x/para2) + para3

popt,pcov = curve_fit(func, x, y, p0=[para1, para2, para3] )
print ("FITING RESULT: ", popt )

#newx = np.linspace(x[0],x[-1],num=500)
#yres = func(newx, *popt)
#res = np.vstack((newx, yres)).T
#
#newfile = "fit_" + file_input
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')
#with open(newfile, 'r') as original: data = original.read()
#with open(newfile, 'w') as modified: modified.write("# " + str(popt[0]) + " " + str(popt[1]) + "\n" + data)

plt.plot(x, y, 'b+:', label='data')
plt.plot(x, func(x, *popt), 'r-', label='fit')
plt.legend()
plt.show()
