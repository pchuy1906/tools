filetmp="job.sh"

for bank in pbronze eci4hps ams enrgmat mxene xanes  ; do
    file="job_${bank}.sh"
    cp $filetmp $file

    sed -i 's/.*#SBATCH -A.*/#SBATCH -A '"$bank"'/' $file
    sed -i 's/.*nMPI=.*/nMPI=224/'                  $file
    sed -i 's/.*#SBATCH -N.*/#SBATCH -N 2/'         $file
    sed -i 's/.*nnodes=.*/nnodes=2/'                $file

sbatch $file &> JOB_ID_${bank}

done

