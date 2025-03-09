import numpy as np
 
rmin = np.loadtxt('rmin.dat')

numcol = np.shape(rmin)[1]
for i in range(numcol):
    arr = rmin[:,i]
    print (len(arr))
    hist = np.histogram(arr, bins=100, range=(0.5,1.5))
    X = hist[1]
    Y = hist[0]
    n = min(len(X), len(Y))
    fname = "hist_"+str(i)+".dat"
    np.savetxt( fname, np.c_[X[:n], Y[:n]])

