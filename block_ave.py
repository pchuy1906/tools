import numpy as np

import argparse
parser = argparse.ArgumentParser(description='calculate average')
# Arguments supported by the code.
parser.add_argument("--file_input",          default='external_pressure.dat', help='file_input')
parser.add_argument("--column_ID", type=int, default=1,                       help='column to calculate average, 1 for 2nd column')
parser.add_argument("--nblock",    type=int, default=1,                       help='number of block to compute average and std')

args      = parser.parse_args()
file_input   = args.file_input
column_ID    = args.column_ID
nblock       = args.nblock


f2 = open("_OOO", "w")
f2.write("%s %s\n" %( "Read input file ", file_input ))

F = np.loadtxt(file_input, usecols=[column_ID])

Ndata = F.shape[0]

f2.write("%s %d\n" %( "No. of data ", Ndata ))

n0 = int(Ndata/nblock)

b_ave = []
b_std = []
for i in range(nblock):
    ibegin = n0*i+0
    iend   = n0*i+n0
    F0 = F[ibegin:iend]
    #print (len(F0))
    ave0 = np.mean(F0)
    std0 = np.std(F0)
    b_ave.append(ave0)
    b_std.append(std0)

b_ave = np.array(b_ave)
b_std = np.array(b_std)

ave = np.mean(b_ave)
print ("average= %15.2f" % ave)
#f2.write("%s %15.2f\n" %( "average=", ave ))

std = np.std(b_ave)
print ("error_bar=%15.9f" % std)
#f2.write("%s %15.2f\n" %( "error_bar=", std ))

f2.write("%s %15.2f %15.2f\n" %( "average=,std=", ave, std ))


res = np.vstack((ave, std)).T
newfile = "results_ave_std.dat"
np.savetxt(newfile, res, fmt='%15.9f %15.9f')

