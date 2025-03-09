cwd=`pwd`
for fileINCAR in `find -name INCAR` ; do
    fold=${fileINCAR/INCAR/}
    echo  $fold
    cd $fold
        rm -rf CHG CHGCAR DFT.FORCE_ENERGY_STRESS DOSCAR EIGENVAL IBZKPT ID_CHECK JOB_ID NON_ORTHO OSZICAR PCDAT REPORT STRESS_TENSOR.dat WAVECAR XDATCAR _KPOINTS _OUTPUT _OUTPUT_KP atomSym.dat box out slurm-* thisSTRESS vasprun.xml quart*
    cd $cwd
done
