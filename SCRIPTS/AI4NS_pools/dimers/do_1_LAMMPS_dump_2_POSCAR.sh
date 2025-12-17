rm -rf POSCAR_*
python ~/tools/others/AI4NS_from_dump.py \
  --file_LAMMPS_dump dump.scan.lammpstrj  \
  --file_LAMMPS_data lammps_data_2molecule.lmp 

