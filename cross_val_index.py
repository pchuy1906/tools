import numpy as np


import argparse
parser = argparse.ArgumentParser(description='split index for cross-validation')
# Arguments supported by the code.
parser.add_argument("--file_input",           default='file.txt', help='file index')
parser.add_argument("--nfold",      type=int, default=4,          help='fold')
parser.add_argument("--prefix",               default='prefix',   help='prefix')

args = parser.parse_args()
file_input = args.file_input
nfold      = args.nfold
prefix     = args.prefix

with open (file_input, "r") as f:
    flist = f.read().splitlines()
#print (flist)

ndata = len(flist)
#print (ndata)

neach = np.round(float(ndata)/float(nfold), 0)
neach = int(neach)
print (neach)

rand_perm = np.random.permutation(ndata)
print (rand_perm)

for i in range(nfold-1):
    fname = "set_"+str(i)+"_"+prefix+".txt"
    ibegin = i*neach
    iend   = i*neach + neach
    f2  = open(fname ,"w")
    for j in range(ibegin, iend):
        idx = rand_perm[j]
        f2.write("%1s\n" %( flist[idx] ))
    f2.close()


i = nfold-1
fname = "set_"+str(i)+"_"+prefix+".txt"
ibegin = i*neach
iend   = ndata
f2  = open(fname ,"w")
for j in range(ibegin, iend):
    idx = rand_perm[j]
    f2.write("%1s\n" %( flist[idx] ))
f2.close()



