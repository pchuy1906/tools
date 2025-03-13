import numpy as np

import argparse
parser = argparse.ArgumentParser(description='Compute mean and std for multiple files')
# Arguments supported by the code.
parser.add_argument('--file_input', nargs='+')
args        = parser.parse_args()
file_input   = args.file_input

nfile = len(file_input)
print (nfile)
ndata = 3
data = np.zeros(shape=(ndata,nfile))

i = 0
for file in file_input:
    #print (file)
    F = np.loadtxt(file)
    data[:,i] = F
    i += 1

f2 = open('final_data.dat', "w")
for idata in range(ndata):
    mean = np.mean(data[idata,:])
    std = np.std(data[idata,:])
    f2.write("%12.4f %12.5f\n" %(mean, std ))
f2.close()

