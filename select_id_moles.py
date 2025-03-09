import numpy as np
 
data = np.loadtxt('TMP_natom.dat').astype(int)
X = data[:,0]
Y = data[:,1]
#print (Y)
nmin=4
nmax=8
ndata=5
id = np.where(np.logical_and(Y>=nmin, Y<=nmax))
if (len(id)<ndata):
    print ("the number of selected molecules is too small ", len(id))
#print (id)

X1 = X[id]
Y1 = Y[id]

id1 = []
for natom in range(nmin, nmax+1):
    #print (natom)
    tmp1 = np.where(Y1==natom)[0]
    if len(tmp1)>0: id1.append(tmp1[0])


import random

if (len(id1)<ndata):
    #print (id1)
    #print (X1[id1])
    #print (Y1[id1])

    sub_id1 = [i for i in range(len(X1)) if i not in id1]
    #print (sub_id1)
    random.shuffle(sub_id1)
    #print (sub_id1)

    nremain = ndata-len(id1)
    for i in range(nremain):
        id1.append(sub_id1[i])


print (id1)
print (X1[id1])
print (Y1[id1])


np.savetxt('selected.dat', X1[id1], fmt='%d')

