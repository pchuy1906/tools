for fold in $(ls -1vd 00*) ; do
    cd $fold
        if [ -f "vasp.out" ] && [ ! -f "job_done" ]; then
            pwd
            rm -rf vasp.out
        fi
    cd ../
done
