import os
import numpy as np

file_b = "b.txt"
print ("read file_b:", file_b)
b = np.loadtxt(file_b)
abs_b = np.absolute(b)
min_abs_b = min(abs_b)
max_abs_b = max(abs_b)
print (min_abs_b, max_abs_b)

#f2 = open("weight_b.dat", "w")

print ("the number of force components:", len(abs_b))
print ("min-b", min_abs_b, " kcal/mol-A")
print ("max-b", max_abs_b, " kcal/mol-A")
print 


x_mean = abs_b.mean()
x_var  = abs_b.var()
x_std  = abs_b.std()

print ("x_mean, x_var, x_std =", x_mean, x_var, x_std)

bins_inp = np.linspace(min_abs_b , max_abs_b, num=100)
hist,bins = np.histogram(abs_b, bins = bins_inp) 
res = np.vstack((bins[:-1], hist)).T
np.savetxt('hist.dat', res, fmt='%15.9f %15.9f')

def weight(E, dE):
    return np.square(dE/(dE+E))

w = weight(abs_b, 2.0*x_std)

np.savetxt('weight_b.dat', w, fmt='%15.9f')

res = np.vstack((abs_b, w)).T
np.savetxt('plot_abs_weight.dat', res, fmt='%15.9f %15.9f')

