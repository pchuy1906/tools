nrun=`ls -d run-* | wc | awk '{print $1}'`
nnext=$(($nrun+1))

echo $nrun $nnext

old_fold="run-$nrun"
new_fold="run-$nnext"

mkdir ${new_fold}
cd ${old_fold}
    cp INCAR KPOINTS job* *.dat ../${new_fold}
    cp CONTCAR ../${new_fold}/POSCAR
cd ../${new_fold}
sbatch job.sh &> JOB_ID
