import numpy as np

import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument('--file_data',  default='aaa.dat')

args    = parser.parse_args()
file_data     = args.file_data

f  = open(file_data ,"rt")

while True:
    #tmp  = f.readline()
    #line = tmp.strip()
    #if line == '': break

    x = []
    y = []
    for i in range(5):
        tmp = f.readline().split()
        #print (tmp)
        if len(tmp) == 0:
            exit()
        x.append(float(tmp[0]))
        y.append(float(tmp[1]))

    x = np.array(x)
    y = np.array(y)

    #print (x, y)
    x0 = np.mean(x)
    dx = 0.5*(x[2]-x[1])
    y0 = np.mean(y)
    dy = 0.5*(y[3]-y[4])
    print (x0, y0, dx, dy)


