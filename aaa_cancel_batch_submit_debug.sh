for fold in $(ls -1vd VASP_*) ; do
    cd $fold
        pwd
        ~/tools/others/job_cancel_batch_submit_debug.sh
    cd ..
    echo
    echo
    echo
done
