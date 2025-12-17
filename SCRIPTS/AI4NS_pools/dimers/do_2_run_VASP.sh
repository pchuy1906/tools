pools="/g/g92/pham20/tools/others/SCRIPTS/AI4NS_pools"

ids=`ls POSCAR_* | sed 's/POSCAR_//' | sort -n`
for i in $ids ; do
    filePOSCAR="POSCAR_$i"
    nfold="VASP_$i"
    if [ ! -d "$nfold" ];  then
        cp -rf $pools/model_VASP $nfold
        mv $filePOSCAR $nfold/POSCAR
        cp ntype.dat $nfold
        cd $nfold
            pwd
            sbatch job.sh &> JOB_ID
        cd ..
    else
        echo "$nfold exists! Do nothing"
    fi
done
