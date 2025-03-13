import sys, numpy, math
from scipy.optimize import leastsq


import argparse
parser = argparse.ArgumentParser(description='read data and fit to the function')
# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='XY data')

args        = parser.parse_args()
file_input     = args.file_input

# Murnaghan equation of state
def eos_murnaghan(params, vol):
    'From Phys. Rev. B 28, 5480 (1983)'
    B0, Bp, V0 = params
    eta = vol/V0
    P = (B0/Bp) * (eta**(-Bp) - 1.0)
    return P

def v_p_eos_murnaghan(params,p):
    B0, Bp, V0 = params
    return V0*numpy.power(p*Bp/B0+1,-1.0/Bp)

def isothermal_compress(params,p):
    B0, Bp, V0 = params
    return (1.0/B0)/(p*Bp/B0+1.0)

# Vinet equation of state
def eos_vinet(params, vol):
    'From Phys. Rev. B 70, 224107'
    B0, Bp, V0 = params 
    eta = (vol/V0)**(1.0/3.0)
    P = 3.0*B0 * (1.0-eta)/eta * numpy.exp(1.5*(Bp-1)*(1.0-eta))
    return P


print "Welcome to eos-PV-fit.py"
print
fname = "volume_P.dat"

f = open(fname, 'rt')

# read data from file
print
print "Data read from file:"
vol = []
press = []
while True:
    line = f.readline().strip()
    if line == '': break
    if line[0] == '#' or line[0] == '!': continue
    p, v = [float(x) for x in line.split()[:2]]
    vol.append(v)
    press.append(p)
    print "v,p=", v, p
print
f.close()

# transform to numpy arrays
vol = numpy.array(vol)
press = numpy.array(press)

# further questions
is_volume = True
vol_unit = "ang3"
press_unit = "GPa"
print

# fit a line to the data and get inital guess for equilibirum volume
# and bulk modulus
a, b = numpy.polyfit(vol, press, 1)
V0 = -b/a
B0 = -a
Bp = 4.0

# initial guesses in the same order used in the Murnaghan function
x0 = [B0, Bp, V0]

def print_params(label, params):
    B0, Bp, V0 = params
    print label, ": V0 = %f angstrom^3" % (V0)
    print label, ": B0 = %f GPa" % (B0)
    print label, ": Bp = %f" % (Bp)
    print

# fit the equations of state
target = lambda params, y, x: y - eos_murnaghan(params, x)
murn, ier = leastsq(target, x0, args=(press,vol))
print_params("Murnaghan", murn)


newx = numpy.linspace(press[0],press[-1],num=500)
newy = v_p_eos_murnaghan( murn , newx)
qcal = isothermal_compress( murn, newx)
res = numpy.vstack((newx, newy, qcal)).T
#
newfile = "fit_" + file_input
numpy.savetxt(newfile, res, fmt='%15.9f %15.9f %15.9f')



