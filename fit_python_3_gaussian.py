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
parser.add_argument("--minxgau3", type=float, default=2.0, help='min-peak-3')
parser.add_argument("--maxxgau1", type=float, default=0.0, help='max-peak-1')
parser.add_argument("--maxxgau2", type=float, default=1.0, help='max-peak-2')
parser.add_argument("--maxxgau3", type=float, default=2.0, help='max-peak-3')


args        = parser.parse_args()
file_input     = args.file_input
xcol           = args.xcol
ycol           = args.ycol
gauwmin        = args.gauwmin
gauwmax        = args.gauwmax
minxgau1       = args.minxgau1
minxgau2       = args.minxgau2
minxgau3       = args.minxgau3
maxxgau1       = args.maxxgau1
maxxgau2       = args.maxxgau2
maxxgau3       = args.maxxgau3


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
gau3_1 = 1.0

gau1_3 = 0.01
gau2_3 = gauwmin
gau3_3 = gauwmin

print ("AAA")
from scipy.signal import argrelextrema
loc_max = argrelextrema(y, np.greater)[0]
val_max = y[loc_max]
print (loc_max, type(loc_max))
print (val_max, type(val_max))

id_sort_val_max = np.argsort(val_max)
print (id_sort_val_max)
print (val_max[id_sort_val_max])

if (len(id_sort_val_max)>2):
    gid = id_sort_val_max[-3:]
    print (gid, type(gid))
    gid = np.sort(gid)
    
    loc_gid = loc_max[gid]
    print (loc_gid)
    
    gau1_2 = x[loc_gid[0]]
    gau2_2 = x[loc_gid[1]]
    gau3_2 = x[loc_gid[2]]
elif (len(id_sort_val_max)==2):
    gid = id_sort_val_max[-2:]
    print (gid, type(gid))
    gid = np.sort(gid)

    loc_gid = loc_max[gid]
    print (loc_gid)

    gau1_2 = x[loc_gid[0]]
    gau2_2 = x[loc_gid[1]]
    gau3_2 = (minxgau3+maxxgau3)*0.5
else:
    print ("ERROR STOP")

print ("guess xvalues of gaussian:")
print (gau1_2, gau2_2, gau3_2)


print ("")
print ("")

#def func(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3):
#    return gau1_1 * np.exp(-(x - gau1_2)**2 / (2 * gau1_3**2)) + \
#           gau2_1 * np.exp(-(x - gau2_2)**2 / (2 * gau2_3**2)) + \
#           gau3_1 * np.exp(-(x - gau3_2)**2 / (2 * gau3_3**2))

def func(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3):
    return gau1_1 / ( (x - gau1_2)**2 / (gau1_3**2) + 1.0 ) + \
           gau2_1 / ( (x - gau2_2)**2 / (gau2_3**2) + 1.0 ) + \
           gau3_1 / ( (x - gau3_2)**2 / (gau3_3**2) + 1.0 )

def func1(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3):
    return gau1_1 / ( (x - gau1_2)**2 / (gau1_3**2) + 1.0 )

def func2(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3):
    return gau2_1 / ( (x - gau2_2)**2 / (gau2_3**2) + 1.0 )

def func3(x, gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3):
    return gau3_1 / ( (x - gau3_2)**2 / (gau3_3**2) + 1.0 )


param_bounds=([0.1,minxgau1,0.001,   0.1,minxgau2,gauwmin,   0.1,minxgau3,gauwmin],[100000,maxxgau1,gauwmax,  100000,maxxgau2,gauwmax,  100000,maxxgau3, gauwmax])
popt,pcov = curve_fit(func, x, y, p0=[gau1_1, gau1_2, gau1_3, gau2_1, gau2_2, gau2_3, gau3_1, gau3_2, gau3_3], maxfev=5000000,bounds=param_bounds )

print ("FITING RESULT: ", popt )

newx = np.linspace(x[0],x[-1],num=500)
yres = func(newx, *popt)
#y1 = func1(newx, *popt)
#y2 = func2(newx, *popt)
#y3 = func3(newx, *popt)
#res = np.vstack((newx, y1)).T
#newfile = "fit_1.dat"
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')
#res = np.vstack((newx, y2)).T
#newfile = "fit_2.dat" 
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')
#res = np.vstack((newx, y3)).T
#newfile = "fit_3.dat" 
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')

res = np.vstack((newx, yres)).T
newfile = "fit_" + file_input
np.savetxt(newfile, res, fmt='%15.9f %15.9f')
#with open(newfile, 'r') as original: data = original.read()
#with open(newfile, 'w') as modified: modified.write("# " + str(popt[0]) + " " + str(popt[1]) + "\n" + data)

#plt.plot(x, y, 'b+:', label='data')
#plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
#plt.legend()
#plt.show()

xgaus = np.array([popt[1],popt[4],popt[7]])
xgaus = np.sort(xgaus)
print (xgaus)

if (len(id_sort_val_max)>2):
    print ("xgau1 = %15.9f" % xgaus[1])
    print ("xgau2 = %15.9f" % xgaus[2])
else:
    print ("xgau1 = %15.9f" % xgaus[0])
    print ("xgau2 = %15.9f" % xgaus[1])

