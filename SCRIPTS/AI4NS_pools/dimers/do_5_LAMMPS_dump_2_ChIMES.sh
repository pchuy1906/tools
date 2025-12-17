rm -rf VASP_2_ChIMES.xyzf
python ~/tools/others/AI4NS_from_dump.py \
  --file_LAMMPS_dump dump_xyz_fxyz  \
  --file_LAMMPS_data lammps_data_2molecule.lmp \
  --file_LAMMPS_log log.lammps
mv VASP_2_ChIMES.xyzf LAMMPS.xyzf
