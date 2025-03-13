import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate the ZPE from list frequencies')

# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='file_input')

args        = parser.parse_args()

file_input   = args.file_input

print ("")
print ("Read input file ", file_input)


f = open(file_input, 'r')
hbar_w = []

while True:
    line = f.readline().strip()
    if line == '': break
    tmp = line.split()
    hbar_w.append(float(tmp[-2]))
f.close()
hbar_w = np.array(hbar_w)

N = len(hbar_w)

kbT = 0.0862 * 300.0

ZPE = 1.0/float(N-3) * ( hbar_w/2.0 + kbT * log(1.0-exp(-hbar_w/kbT)) )

print ("ZPE in meV")
print (sum(ZPE[3:]))
