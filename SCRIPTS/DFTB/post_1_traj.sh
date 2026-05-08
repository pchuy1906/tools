~/tools/DFTB/do_xyz_2_xyz.sh 
a=`head -2 file.xyz | tail -1 | awk '{printf "%15.9f\n", $1}'`
b=`head -2 file.xyz | tail -1 | awk '{printf "%15.9f\n", $5}'`
c=`head -2 file.xyz | tail -1 | awk '{printf "%15.9f\n", $9}'`

echo $a $b $c
~/tools/DFTB/DFTB-xyz_2_goodxyz/a.out  file.xyz  $a $b $c
mv newfile.xyz  traj.xyz
#rm file.xyz 
