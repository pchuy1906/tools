filegen="geo_end.gen"
ntype=1

natom=`head -1 $filegen | awk '{print $1}'`
~/tools/others/genS_2_xyz_box9/a.out  $natom $ntype  geo_end.gen 
~/tools/others/do_xyz_2_qe.sh  traj.xyz 
