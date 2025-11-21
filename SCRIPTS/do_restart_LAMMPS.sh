cwd=`pwd`

nrun=`ls -d run-* | wc | awk '{print $1}'`
echo $nrun
oldrun="run-${nrun}"
newrun="run-$(($nrun+1))"
echo $oldrun $newrun

mkdir $newrun
cd run-1
    cp job.sh in.lammps ../$newrun
cd ../$newrun
    sed -i 's|.*read_restart.*|read_restart     ../'${oldrun}'/restart.a|g'  in.lammps
    sed -i 's|.*read_data.*|read_restart     ../'${oldrun}'/restart.a|g'     in.lammps

    sed -i 's|read_dump|#read_dump|g'                                        in.lammps
    sed -i 's|velocity|#velocity|g'                                          in.lammps


    JOBID=`tail -1 ../${oldrun}/JOB_ID | awk '{print $NF}'`
    sbatch --dependency=afterany:$JOBID job.sh &> JOB_ID
cd ..
