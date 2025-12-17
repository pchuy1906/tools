#!/bin/sh

#SBATCH -N 1
#SBATCH -J DFT
#SBATCH -t 01:00:00
#SBATCH -p pbatch
#SBATCH -A enrgmat
#SBATCH --exclusive

nnodes=1
nMPI=112

POTCARpool="/usr/gapps/emc-vasp/pseudopotentials/potpaw_PBE"

rm -rf POTCAR
while read type ; do
    echo $type
    cat ${POTCARpool}/${type}/POTCAR >> POTCAR
done < ntype.dat

VASPexe="/usr/gapps/emc-vasp/vasp.6.3.0_vtst/bin/vasp_gam"
if [ ! -f out ]; then
    srun -N $nnodes -n $nMPI $VASPexe > out
fi
