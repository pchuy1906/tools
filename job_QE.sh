#!/bin/sh

#SBATCH -N 2
#SBATCH -J QE
#SBATCH -t 01:00:00
#SBATCH -p pdebug
#SBATCH -A pbronze
#SBATCH --exclusive

nnodes=2
nMPI=112

exe="/p/lustre2/pham20/SiO2/QE_from_Marcos/pw.x"
date
srun -N $nnodes -n $nMPI $exe -in qe.in >& output
date

