import sys
import numpy as np
import scipy.linalg
import math as m
import subprocess
import os

import argparse

from numpy import *
from numpy.linalg import lstsq
from datetime import *
from subprocess import call

import glob

pairXX = ['N_N','N_H','H_H']

for XX in pairXX:
    sym = 'range-*/rdf*/gr_' + XX + '.dat'
    print (sym)
    #listA = glob.glob('range-*/rdf*/gr_H_H.dat')
    listA = glob.glob(sym)
    listA = sorted(listA)
    #print (listA)
    
    ncount = 0
    for Afile in listA:
        print (Afile)
        A = np.genfromtxt(Afile, dtype='float')
        ncount += 1
        if (ncount == 1):
            res = A[:,1]
            x = A[:,0]
        else:
            res = np.c_[res, A[:,1]]
    print (res.shape)
    
    y_lo = np.zeros(res.shape[0])
    y_hi = np.zeros(res.shape[0])
    yave = np.zeros(res.shape[0])
    for k in range(res.shape[0]):
        y = res[k,:]
        y_std = np.std(y)
        y_ave = sum(y)/len(y)
        y_lo[k]  = (y_ave-y_std)
        y_hi[k]  = (y_ave+y_std)
        yave[k]  = (y_ave)
    
    data = np.transpose([x, yave, y_lo, y_hi])
    print (data.shape)
    np.savetxt( "gr_" +XX + ".dat", data)
    
    
    x_r = np.flipud(x)
    y_hi_r = np.flipud(y_hi)
    newx = np.concatenate((x,x_r))
    newy = np.concatenate((y_lo,y_hi_r))
    print (newx.shape, newy.shape)
    data = np.transpose([newx, newy])
    np.savetxt( "error_bar_" + XX + ".dat", data)

