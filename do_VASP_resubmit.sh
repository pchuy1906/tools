file_job="job.sh"

nrun=`ls -d run-* | wc | awk '{print $1}'`
echo ${nrun}

nextrun=$(($nrun+1))

mkdir run-${nextrun}

cp  run-${nrun}/${file_job}    run-${nextrun}
cp  run-${nrun}/INCAR      run-${nextrun}
sed -i 's/ISTART  = 0/ISTART  = 1/g'  run-${nextrun}/INCAR
cp  run-${nrun}/KPOINTS    run-${nextrun}
cp  run-${nrun}/POTCAR     run-${nextrun}
cp  run-${nrun}/CONTCAR    run-${nextrun}/POSCAR
cp  run-${nrun}/WAVECAR    run-${nextrun}
cp  run-${nrun}/ntype.dat  run-${nextrun}

cd run-${nextrun}
    sbatch ${file_job} &> JOB_ID
cd ..
