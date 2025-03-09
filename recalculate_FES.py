import sys
import numpy
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

xfile = 'x.txt'

listA = glob.glob('A.*.txt')
listA = sorted(listA)
#print (listA)

x = numpy.genfromtxt(xfile, dtype='float')

f = open('new_Ax.txt','ab')
for Afile in listA:
    print (Afile)
    A = numpy.genfromtxt(Afile, dtype='float')
    y = dot(A,x)
    numpy.savetxt( f, y)
