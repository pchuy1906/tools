cwd=`pwd`
for file in `find -name run_md.in` ; do
    fold=${file/run_md.in/}
    echo  $fold
    cd $fold
        rm -rf output.xyz output.xyz.bak restart.bak restart.xyzv traj_bad_r.lt.rin+dp.xyz traj_bad_r.lt.rin.xyz vmd.vmd slu*
    cd $cwd
done
