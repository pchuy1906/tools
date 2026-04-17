cwd=`pwd`
for fold in $(ls -1vd LAMMPSSimulator/condensed_*_*_* ) ; do
    cd $fold
        pwd
        diff 00000/LMP_condensed_pair.lammps $cwd/LAMMPSSimulator/dimer_*/00000/LMP_dimer_pair.lammps
    cd $cwd
done

