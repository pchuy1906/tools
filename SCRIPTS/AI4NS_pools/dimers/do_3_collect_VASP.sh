cwd=`pwd`

rm -rf DFT.xyzf
for fold in $(ls -1vd VASP*/) ; do
    cd $fold
        pwd
        rm -rf VASP_2_ChIMES.xyzf
        python ~/tools/others/AI4NS_VASP_2_ChIMES.py \
          --vasp_outcar OUTCAR \
          --vasp_poscar POSCAR \
          --cell_type NON_ORTHO \
          --lammps_ids ../LAMMPS_ID_order_type_mol.dat
        cat VASP_2_ChIMES.xyzf >> $cwd/DFT.xyzf
    cd $cwd
done
