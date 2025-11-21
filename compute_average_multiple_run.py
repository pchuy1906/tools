import numpy as np

import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument('--file_data',   nargs='+')
parser.add_argument('--file_weight', nargs='+')
parser.add_argument('--data_type',  type=int,   default=1, help='1/1d-array 2/2d-array')

args    = parser.parse_args()
file_data     = args.file_data
file_weight   = args.file_weight
data_type     = args.data_type

ndata = len(file_data)
print (ndata)

if (data_type==1):
    ydatas = np.array([])
    weights = []
    for i in range(ndata):
        file1 = file_data[i]
        file2 = file_weight[i]
        tmp1_data = np.loadtxt(file1)
        ydatas = np.append(ydatas, tmp1_data[:,1])
        tmp2_data = np.loadtxt(file2)
        weights.append(tmp2_data)
    
    xdata = tmp1_data[:,0]
    nx = len(xdata)
    ydatas = ydatas.reshape((ndata, nx))
    
    
    f2 = open("average_1d.dat", "w")
    for i in range(nx):
        tmp_y = ydatas[:,i]*weights
        yave = sum(tmp_y)/sum(weights)
        f2.write("%15.9f %15.9f\n" %( xdata[i], yave ))


if (data_type==2):
    ydatas = np.array([])
    weights = []
    for i in range(ndata):
        file1 = file_data[i]
        file2 = file_weight[i]
        # the first line is comment, skip it
        tmp1_data = np.loadtxt(file1,skiprows=1)
        ydatas = np.append(ydatas, tmp1_data[:,2])
        tmp2_data = np.loadtxt(file2)
        weights.append(tmp2_data)

    xdata = tmp1_data[:,0:2]
    nx = len(xdata)
    ydatas = ydatas.reshape((ndata, nx))


    fline = open(file1).readline().rstrip()
    #print (fline)
    f2 = open("average_contour.dat", "w")
    f2.write("%s\n" %( fline ))
    for i in range(nx):
        tmp_y = ydatas[:,i]*weights
        yave = sum(tmp_y)/sum(weights)
        f2.write("%15.9f %15.9f %15.9f\n" %( xdata[i,0], xdata[i,1], yave ))

