#!/usr/bin/env python3
import numpy as np
from scipy import signal
import math, os, sys, time
import matplotlib.pyplot as plt


##### PLEASE READ THE FOLLOWING INSTRUCTIONs BEFORE RUNNING SCRIPT #####
####                                                                ####
####  The Format for Running This Script:                           ####
####  python IR_total_KW.py INPUT_FILE DELTA_T WINDOW OUTPUT_FILE   ####
####                                                                ####
####  The values need to input manually when runing this script     ####
####                                                                #### 
####  (1) INPUT_FILE_NAME: The Total_Dipole_Moment_*.Diople file    ####
####           (NOTE: do NOT need to re-split the Dipole file)      ####
####                                                                #### 
####  (2) DELTA_T: The Time_step set in simulation, in unit of fs   ####
####                                                                ####
####  (3) WINDOW: The Name of the Window Function                   ####
####                                                                ####  
####  (4) OUTPUT_FILE_NAME: The Name of the Output File.            ####
####           (NOTE: do NOT need to type > sign!)                  ####
####                                                                ####
#############################  Let's Try It! ###########################






#### The values need to input manually when running this script ####
#path = input("\nPlease Enter the Directory Contained the Dipole Moment File\n")
fname = sys.argv[1]                    # The name of the input file
delta_t = float(sys.argv[2]) * 1.0e-15 # The time step in unit of femtoseconds
window = sys.argv[3]                   # The name of the window function
fout = sys.argv[4]                     # The name of the output file



#### The constants will be used in this script ####
EPSILON0=8.854187817*1.0e-12       # Vacuum permittivity.  C^2 N^-1 m^-2
kB = 1.38064852*1.0e-23            # m^-1 N^-1 K-1
T = 300                            # K
beta = 1.0/(kB * T)                # m^-1 N^-1
hbar = 1.0545718*1.0e-34           # m^2 kg / s = N * s * m
betahbar = beta*hbar               # s
#V = np.load('box.npy')
#V = np.mean(V[:,0]*V[:,4]*V[:,8]) * 1.0e-30     #  m^3
V=1000. #FIXME
c = 299792458     # m/s
PI = np.pi
# 1 e = 1.602*1.0e-19 C
# change unit to C*m for M(0)
unit_basic = 1.602176565*1.0e-19*1.0e-10;       # C m
# change unit to s for dM(0)/dt
unit = unit_basic                               # C m s^-1
unit2 = unit * unit                             # C^2 m^2 s^-2
unit_all = 2.0*PI*beta/3.0/c/V*unit2            # m^-3 N^-1 * s^-1 * C^2

# 1F = A^2 * s^4 * kg^-1 * m^-2
pref = unit_all/(4*PI*EPSILON0);                # m^-1 * s^-1, s^-1 taken by dt.
pref /= 100     # change m^-1 to cm^-1
pref /= 1000    # change cm^-1 to 1000cm^-1



#### Functions will used in this script ####

def read_data(fname,dt):
  vac=np.loadtxt(fname)
  vac=vac.T[1:].T
  return vac

def zero_padding(sample_data):
    '''
      A series of Zeros will be padded to the end of the dipole moment array 
    (before FFT performed), in order to obtain a array with the length which
    is the "next power of two" of numbers.
    #### Next power of two is calculated as: 2**np.ceil(log2(x))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    '''
    N = 2**int(math.log(len(sample_data)*2-1, 2))
    return N

def choose_window(data, kind='string',delta_t = 1e-15):
    ii = int(1e-12 / delta_t) * 0.5
    if kind == 'Gaussian':
        sigma = 2 * math.sqrt(2 * math.log(2))
        window = signal.gaussian(len(data)*2, std=ii/sigma, sym=False)[len(data):]
    elif kind == 'BH':
        window = signal.blackmanharris(len(data), sym=False)
    elif kind == 'Hamming':
        window = signal.hamming(len(data), sym=False)
    elif kind == 'Hann':
        window = signal.hann(len(data), sym=False)
    return window

def calc_ACF(array, window,delta_t):
    '''
    This function deals with the auto-correlation function (ACF) of the total
    dipole moment derivatives.

    With the Wiener-Khintchine theorem, the autocorrelation function is
    http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem

####
####  http://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
####
####  for fast convolution 
####  http://sebug.net/paper/books/scipydoc/frequency_process.html#id5
    '''
    # normalization
    #yunbiased = array - np.mean(array, axis=0)
    #ynorm = np.sum(np.power(yunbiased,2), axis=0)
    #print("the average value of input data array:")
    #print(ynorm)
    array -= np.mean(array)
    autocor = np.zeros(np.shape(array))

    for i in range(3):
        #autocor[:,i] = signal.fftconvolve(array[:,i],
        #                                  array[:,i][::-1],
        autocor[:,i] = signal.correlate(array[:,i],
                                        array[:,i],
                                        mode='full')[len(array)-1:] / np.shape(array)[0] #/ynorm[i]
    print("shape of the result3 from signal.FFTcorrelate()", np.shape(autocor))
    window = choose_window(autocor, kind=window, delta_t = delta_t)
    #WE = sum(window) / len(autocor)
    #print(window)
    #wf = window / WE
    #print(wf)
    # convolve the window function. 
    autocor = autocor * window[None,:].T
    return autocor

def calc_ACF2(array, window,delta_t):
    '''
    This function deals with the auto-correlation function (ACF) of the total
    dipole moment derivatives.

    With the Wiener-Khintchine theorem, the autocorrelation function is
    http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem

####
####  http://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
####
####  for fast convolution 
####  http://sebug.net/paper/books/scipydoc/frequency_process.html#id5
    '''
    # normalization
    incor = int(1e-12/delta_t)
    yunbiased = array - np.mean(array, axis=0)
    nframe = np.shape(array)[0]
    autocor = np.zeros((incor,3))
    nstat = np.shape(yunbiased[incor:,:])[0]

    for i in range(incor):
        autocor[i,:] = np.sum(yunbiased[:nstat,:] * yunbiased[i:i+nstat,:], axis = 0) / nstat
        

    window = choose_window(autocor, kind=window, delta_t = delta_t)
    WE = sum(window) / len(autocor)
    wf = window / WE
    # convolve the window function. 
    autocor = autocor * wf[None,:].T
    return autocor


def calc_FFT(sig, window,delta_t):
    '''
    This function is for calculating the "intensity" of the ACF at each 
    frequency by using the discrete fast Fourier transform.
    
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    '''
    # A series of number of zeros will be padded to the end of the DACF \
    # array before FFT.
    N = zero_padding(sig)
    window = choose_window(sig, kind=window, delta_t = delta_t)
    WE = sum(window) / len(sig)
    wf = window / WE
    # convolve the window function. 
    sig = sig * window[None,:].T
	
    yfft = np.fft.fft(sig, N, axis=0) # / len(sig)
    print(np.shape(yfft))
    print(int(len(yfft)/2))
    l=yfft.shape[0]
    wavenumber = np.fft.fftfreq(l, delta_t*c*100)[0:int(l/2)]
    return np.real(yfft[0:int(l/2),:]), wavenumber
    #yfft = np.fft.fft(sig, N) # / len(sig)
    #print(np.shape(yfft))
    #print(int(len(yfft)/2))
    #l=yfft.shape[0]
    #wavenumber = np.fft.fftfreq(l, delta_t)/c/100
    ##return np.real(yfft[0:int(l/2)]), wavenumber[0:int(l/2)]
    #return yfft.real, wavenumber

def calc_FFT2(sig, window, delta_t):
    '''
    This function is for calculating the "intensity" of the ACF at each 
    frequency by using the discrete fast Fourier transform.
    
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    '''
    # A series of number of zeros will be padded to the end of the DACF \
    # array before FFT.
    LL = 5000
    c = 299792458 
    wavenumber = np.array(range(LL))
    yfft = np.array(range(LL)) * 0.0
    tt = np.array(range(len(sig))) * delta_t
    for ii in range(LL):
        yfft[ii] = np.sum(2 * np.cos(wavenumber[ii] * 2*np.pi / 0.01 * c * tt) * sig[ii])

    return yfft, wavenumber

######## Save The Results to A TEXT File ########
def save_results(fout, wavenumber, intensity):
    title = ("Wavenumber", "Intensity")
    np.savetxt(fout,np.c_[wavenumber, intensity])

######## Plot The Spectrum by Using Matplotlib module ########
def visualization(D_p, DACF, wavenumber, intensity, delta_t):
    plt.subplot(3,1,1)
    L1 = np.arange(len(D_p))
    plt.plot(L1*delta_t*1e+12, D_p[:,0], color='red', linewidth=1.5)
    plt.plot(L1*delta_t*1e+12, D_p[:,1], color='green', linewidth=1.5)
    plt.plot(L1*delta_t*1e+12, D_p[:,2], color='blue', linewidth=1.5)
    plt.axis([0, len(D_p)*delta_t*1e+12, 1.1*np.min(D_p), 1.1*np.max(D_p)])
    plt.xlabel("time (ps)", fontsize=15)
    plt.ylabel("Derivative of Dipole (a.u.)", fontsize=15)

    plt.subplot(3,1,2)
    L2 = np.arange(len(DACF))
    plt.plot(L2*delta_t*1e+12, DACF[:], color='red', linewidth=1.5)
    plt.axis([0, 2, 1.1*np.min(DACF), 1.1*np.max(DACF)])
    plt.xlabel("time (ps)", fontsize=15)
    plt.ylabel("DACF (a.u.)", fontsize=15)

    plt.subplot(3,1,3)
    plt.plot(wavenumber, intensity[:,0], color='black',  linewidth=1.5, label='x')
    plt.plot(wavenumber, intensity[:,1], color='orange', linewidth=1.5, label='y')
    plt.plot(wavenumber, intensity[:,2], color='green',  linewidth=1.5, label='z')
    plt.axis([0, 6000,
             -1.1*np.min(intensity), 1.1*np.max(intensity)])
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
    plt.subplots_adjust(hspace = 0.5)
    plt.legend()
    plt.show()

######## The main program ########
def main(fname, delta_t, window, fout):
    #start = time.clock()
    vac = read_data(fname,delta_t)
    vac[0,:] /= 2
    yfft,wavenumber = calc_FFT(vac, window, delta_t)
    yfft /= delta_t 
    intensity = pref * yfft / 100
    save_results(fout, wavenumber, 2 * intensity)
    fmt = "\n Work Completed! Used Time: {:.3f} seconds"

if __name__ == '__main__':
    main(fname, delta_t, window, fout)
