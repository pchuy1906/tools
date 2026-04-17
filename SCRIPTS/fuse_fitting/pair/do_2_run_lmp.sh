
cwd=`pwd`
for fold in $(ls -1vd LAMMPSSimulator/condensed_*_*_*_* ) ; do
    cd $fold

        pwd
        cp LMP_condensed_* 00000
        cp LMP_monomer_1_* 00001
        cp LMP_monomer_2_* 00002

        cd 00000
            full_path=`grep -e "initial" lammps.in | head -1 | awk '{print $2}'`
            file_name="initial.dump"
            sed -i "s|$full_path|$file_name|g" lammps.in
            /usr/workspace/fuse/lammps/lammps-2Apr2025/build_jas_tuo_11_18_2025/lmp -in lammps.in -log lammps.out
        cd ..

        cd 00001
            full_path=`grep -e "initial" lammps.in | head -1 | awk '{print $2}'`
            file_name="initial.dump"
            sed -i "s|$full_path|$file_name|g" lammps.in
            /usr/workspace/fuse/lammps/lammps-2Apr2025/build_jas_tuo_11_18_2025/lmp -in lammps.in -log lammps.out
        cd ..

        cd 00002
            full_path=`grep -e "initial" lammps.in | head -1 | awk '{print $2}'`
            file_name="initial.dump"
            sed -i "s|$full_path|$file_name|g" lammps.in
            /usr/workspace/fuse/lammps/lammps-2Apr2025/build_jas_tuo_11_18_2025/lmp -in lammps.in -log lammps.out
        cd ..

        ls -ltr */lammps.out

    cd $cwd
done

