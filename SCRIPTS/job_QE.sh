#!/bin/sh

#SBATCH -N 1
#SBATCH -J QE
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH -A guests
#SBATCH --exclusive

nnodes=1
nMPI=112

exe="/usr/workspace/pham20/codes/q-e-qe-7.0/bin/pw.x"

srun -N $nnodes -n $nMPI $exe -in qe.in >& output
