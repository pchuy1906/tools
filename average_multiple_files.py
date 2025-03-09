import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='average multiple files')
# Arguments supported by the code.
parser.add_argument('--input_files', nargs='+')

args    = parser.parse_args()
input_files          = args.input_files

#print (input_files)

ndatas = []
for ifile in input_files:
    #print (ifile)
    tmp = np.loadtxt(ifile)
    ndata = tmp.shape[0]
    ndatas.append(ndata)

#print (ndatas)
ntake = min(ndatas)
print (ntake)



bmatrix = np.array([])
for ifile in input_files:
    print (ifile)
    tmp = np.loadtxt(ifile)
    bmatrix = np.append(bmatrix, tmp[:ntake,1])
    
#print (bmatrix)
bmatrix = bmatrix.reshape(( len(input_files), ntake))
#print (bmatrix)

f2 = open("average.dat", "w")
for i in range(ntake):
    yvalues = bmatrix[:,i]
    #print (yvalues)
    #print ( tmp[i,0], np.mean(yvalues), np.std(yvalues) )
    f2.write("%15.9f %15.9f %15.9f\n" %( tmp[i,0], np.mean(yvalues), np.std(yvalues) ))
