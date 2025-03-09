import numpy as np
 
rmin = np.loadtxt('rmin.dat')
k = 6

Amatrix = np.array([])
numcol = np.shape(rmin)[1]
for i in range(numcol):
    arr = rmin[:,i]
    #print (len(arr))
    #idx = np.argpartition(arr, k)
    idx = np.argsort(arr)[:k]
    #print (idx)
    #print (arr[idx])
    Amatrix = np.append(Amatrix, idx)

#Amatrix = Amatrix.reshape((nconf, npair))
#np.savetxt('rmin.dat', Amatrix, fmt = '%.6f')
Amatrix = np.array(Amatrix, dtype=int)
#print (Amatrix)

Amatrix = np.sort(Amatrix)
#print (Amatrix)

Amatrix = np.unique(Amatrix)
#print (Amatrix)


np.savetxt('id_rmin.dat', Amatrix, fmt = '%d')

