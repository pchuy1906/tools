cp LAMMPS.xyzf file.xyzf
python ~/tools/others/AI4NS_fit_LJ.py \
  --file_xyzf file.xyzf  \
  --n_type_max 21 \
  --rcut 10 \
  --train_forces
