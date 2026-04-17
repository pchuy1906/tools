rm -rf JOB_IDs

for i in $(seq 1 6) ; do
    for file in $(ls flux*.sh) ; do
        sbatch $file >> JOB_IDs
    done
done
