for fold in $(ls -1vd VASP_*) ; do
    cd $fold
        pwd
        ~/tools/others/job_cancel_debug_submit_batch.sh
    cd ..
    echo
    echo
    echo
done
