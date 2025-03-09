import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Compute Ax')

# Arguments supported by the code.
parser.add_argument("--fileA", default='A.txt', help='fileA')
parser.add_argument("--fileX", default='x.txt', help='fileX')
args        = parser.parse_args()
fileA   = args.fileA
fileX   = args.fileX

A = np.loadtxt(fileA)
print A.shape, A.shape[0], A.shape[1]

x = np.loadtxt(fileX)
print x.shape, x.shape[0]

Ax = np.dot(A,x)
np.savetxt('iAx.dat', Ax, fmt='%15.9f')

