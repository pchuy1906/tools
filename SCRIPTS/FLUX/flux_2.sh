#!/bin/bash
#flux: -N 5
#flux: -q pbatch
#flux: -B fuse
#flux: -t 60m
#flux: --exclusive



for fold in $(ls -1vd 002*) ; do
    cd $fold
        if [ ! -f "vasp.out" ]; then
            flux run -N5 -n480 /usr/gapps/qsg/codes/VASP/tuolumne/v5.4.4/vasp_gam  > vasp.out
            touch job_done
        fi
    cd ..
done
