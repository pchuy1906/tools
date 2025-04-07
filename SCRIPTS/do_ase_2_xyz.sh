#source ~/.BASHRC_activate_conda
#conda activate lmp_python

aaa=`head -5 POSCAR | tail -3 | xargs`

python ~/tools/others/ase_2_xyz.py \
  --traj_ase traj.ase \
  --cell_9 $aaa

mv result_traj.xyz  traj.xyz
