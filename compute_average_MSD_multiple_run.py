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
array_MSD = np.full((nrow, ncol), valueToFill)

ntrunk = 8

real_MDtime = []

array_D = []

for i in range(nfile):

    file1 = file_data[i]
    tmp = np.loadtxt(file1)

    # cut the tail
    ndata = len(tmp[:,0])
    nbegin = 0
    nend   = ndata * (ntrunk-1)// ntrunk

    MSD = tmp[nbegin:nend,2]
    # fill the array
    array_MSD[nbegin:nend,i] = MSD

    MDtime = tmp[nbegin:nend,0]
    if len(real_MDtime) < len(MDtime):
        real_MDtime = MDtime

    # cut the head and the tail
    nbegin = ndata * 1// ntrunk
    nend   = ndata * (ntrunk-1)// ntrunk
    array_D.extend(tmp[nbegin:nend,1])

ave_D = np.mean(array_D)/time_interval_fs*1e5
std_D =  np.std(array_D)/time_interval_fs*1e5

print (f"mean and std of diffusions (unit 1e-10 * m2/s) {ave_D} {std_D}")

def compute_average(array_MSD):
    num_rows, num_cols = array_MSD.shape
    arr_ave_MSD = []
    arr_std_MSD = []
    for i in range(num_rows):
        MSD = array_MSD[i,:]
        positive_MSD = MSD[MSD > -0.5]
        if ( len(positive_MSD) > 0 ):
            ave_MSD = np.mean(positive_MSD)
            std_MSD = np.std(positive_MSD)
        else:
            ave_MSD = -100.0
            std_MSD = -100.0
        arr_ave_MSD.append(ave_MSD)
        arr_std_MSD.append(std_MSD)
    return arr_ave_MSD, arr_std_MSD

print (len(real_MDtime))
ave_MSD, std_MSD = compute_average(array_MSD)

f2 = open("average_MSD.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_MSD[i] ))
f2.close()

f2 = open("average_MSD_low.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_MSD[i]-std_MSD[i] ))
f2.close()

f2 = open("average_MSD_hi.dat", "w")
for i in range(len(real_MDtime)):
    f2.write("%15.9f %15.9f\n" %( real_MDtime[i], ave_MSD[i]+std_MSD[i] ))
f2.close()

