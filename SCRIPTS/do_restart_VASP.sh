nrun=`ls -d run-* | wc | awk '{print $1}'`
nnext=$(($nrun+1))

echo $nrun $nnext

old_fold="run-$nrun"
new_fold="run-$nnext"

mkdir ${new_fold}


cat > job_${nnext}.sh << EOF
#!/bin/sh

#SBATCH -N 1
#SBATCH -J DFT-3
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH -A guests
#SBATCH --exclusive

nnodes=1
nMPI=112

POTCARpool="/usr/gapps/emc-vasp/pseudopotentials/potpaw_PBE.54"

cd ${old_fold}
    cp INCAR KPOINTS POTCAR ../${new_fold}
    cp CONTCAR ../${new_fold}/POSCAR
cd ../${new_fold}
    #rm -rf POTCAR
    #while read type ; do
    #    echo \$type
    #    cat \${POTCARpool}/\${type}/POTCAR >> POTCAR
    #done < ntype.dat
    
    srun -N \$nnodes -n \$nMPI /usr/gapps/emc-vasp/vasp.6.3.0_vtst/bin/vasp_gam > out_SCAN

EOF

JOBID=`tail -1 ${old_fold}/JOB_ID | awk '{print $NF}'`
sbatch --dependency=afterany:$JOBID job_${nnext}.sh &> ${new_fold}/JOB_ID

