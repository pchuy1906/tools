#!/bin/sh

#MSUB -N cp2k
#MSUB -l nodes=4
#MSUB -l walltime=00:00:29:00
#MSUB -A pbronze
#MSUB -q pdebug
#MSUB -V

input="cp2k.inp"
output="cp2k.out"
exe="/g/g92/pham20/codes/CP2K/cp2k-8.2/exe/Linux-x86-64-gfortran/cp2k.popt"
rm -rf ${output}
srun -n 144  ${exe}   -o  ${output}    ${input}

