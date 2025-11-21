import numpy as np

import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument('--file_data',   nargs='+')
parser.add_argument('--time_interval_fs', type=float, default=2.5)


args    = parser.parse_args()
file_data     = args.file_data
time_interval_fs = args.time_interval_fs

nfile = len(file_data)

## create an array with all -100.0
nrow = 999999
ncol = nfile
valueToFill = -100.0
array_C2t = np.full((nrow, ncol), valueToFill)

ntrunk = 8

real_MDtime = []

for i in range(nfile):

    file1 = file_data[i]
    tmp = np.loadtxt(file1)

    # cut the tail
    ndata = len(tmp[:,0])
    nbegin = 0
    nend   = ndata * (ntrunk-1)// ntrunk

    C2t = tmp[nbegin:nend,1]
    # fill the array
    array_C2t[nbegin:nend,i] = C2t

    MDtime = tmp[nbegin:nend,0]
    if len(real_MDtime) < len(MDtime):
        real_MDtime = MDtime

def compute_average(array_C2t):
    num_rows, num_cols = array_C2t.shape
    arr_ave_C2t = []
    arr_std_C2t = []
    for i in range(num_rows):
        C2t = array_C2t[i,:]
        positive_C2t = C2t[C2t > -0.5]
        if ( len(positive_C2t) > 0 ):
            ave_C2t = np.mean(positive_C2t)
            std_C2t = np.std(positive_C2t)
        else:
            ave_C2t = -100.0
            std_C2t = -100.0
        arr_ave_C2t.append(ave_C2t)
        arr_std_C2t.append(std_C2t)
    return arr_ave_C2t, arr_std_C2t

print (len(real_MDtime))
real_MDtime = real_MDtime * time_interval_fs * 0.001
ave_C2t, std_C2t = compute_average(array_C2t)

f2 = open("average_C2t.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_C2t[i] ))
f2.close()

f2 = open("average_C2t_low.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_C2t[i]-std_C2t[i] ))
f2.close()

f2 = open("average_C2t_hi.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_C2t[i]+std_C2t[i] ))
f2.close()

