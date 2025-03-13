import numpy as np

import argparse
parser = argparse.ArgumentParser(description='find the H(p) from the fit')
# Arguments supported by the code.
parser.add_argument("--A",   type=float, default=100, help='A')
parser.add_argument("--The", type=float, default=100, help='Theta')

parser.add_argument("--vmin", type=float, default=100, help='Volume min')
parser.add_argument("--vmax", type=float, default=400, help='Volume max')

args        = parser.parse_args()
vmin     = args.vmin
vmax     = args.vmax
A        = args.A
The      = args.The

# equation
def exp_cp(params, T):
    'From Joel Christenson"s slide'
    A, The = params 
    x = The/T
    vexp = np.exp(x)
    cp = A*x*x*vexp/((vexp-1.0)*(vexp-1.0))
    return cp

npoints = 500
vfit = np.linspace(vmin, vmax, npoints)

params = [A, The]
yfit = exp_cp(params, vfit)
res = np.vstack((vfit, yfit)).T
newfile = "gen_data.dat"
np.savetxt(newfile, res, fmt='%15.9f %15.9f')




