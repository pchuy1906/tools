import os
import math
import numpy as np
from numpy import *

from numpy import linalg as LA

import argparse

parser = argparse.ArgumentParser(description='Calculate max abs')


# Arguments supported by the code.
parser.add_argument("--file_input", default='AAA.dat', help='file_input')

args        = parser.parse_args()

file_input   = args.file_input

F = np.loadtxt(file_input)
F = F[:len(F)-2]

F = np.absolute(F)

print F.max(), LA.norm(F), sum(F)

