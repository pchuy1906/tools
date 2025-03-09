import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate weights')
# Arguments supported by the code.
parser.add_argument("--nw",       type=float,   default=1000.0)

args    = parser.parse_args()
nw       = args.nw



file_b = "b.txt"
print ("read file_b:", file_b)
b = np.loadtxt(file_b)
abs_b = np.absolute(b)
min_abs_b = min(abs_b)
max_abs_b = max(abs_b)
print (min_abs_b, max_abs_b)

print ("the number of force components:", len(abs_b))
print ("min-b", min_abs_b, " kcal/mol-A")
print ("max-b", max_abs_b, " kcal/mol-A")
print 


w = nw/(abs_b+0.001)

np.savetxt('weight_b.dat', w, fmt='%15.9f')

