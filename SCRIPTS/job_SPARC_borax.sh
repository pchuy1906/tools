#!/bin/sh

#SBATCH -N 1
#SBATCH -J sparc
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH -A pbronze
#SBATCH --exclusive

nnodes=1
nMPI=36

exe="/usr/workspace/pham20/codes/SPARC/SPARC_2025_08_19/lib/sparc"

srun -N $nnodes -n $nMPI $exe -name SPARC
