import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate RMS')

# Arguments supported by the code.
parser.add_argument("--file_input", default='COMPARE.dat', help='file_input')

args        = parser.parse_args()

file_input   = args.file_input

print "input file is: ", file_input
print "output files:"

fh = open(file_input)
#the first 3 lines are comments
fh.readline()
fh.readline()
fh.readline()

ncount = 0
while True:
    line = fh.readline()
    tmp = line.split()
    #print(len(tmp))
    if (len(tmp)==3):
        ncount += 1
        #file_output = file_input + "." + str(ncount)
        prefix = file_input.split('.')[0]
        file_output = "splitted_" + prefix + "." + tmp[0]
        print (file_output)
        OUTFILE = open(file_output,'w') 
    if (len(tmp)==4):
        aa =line.split()
        data_line = (aa[1]+" "+aa[3]+"\n")
        OUTFILE.write(data_line)
    if not line:
        break
fh.close()
