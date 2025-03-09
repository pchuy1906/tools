aaa=`head -5 POSCAR | tail -3 | xargs`
python ~/tools/others/ase_2_xyz.py --cell_9 $aaa --traj_ase opt_traj.ase
