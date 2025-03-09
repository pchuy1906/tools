echo
echo "usage: \"do_script.sh  file.xyz \""
echo

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filexyz=$1

natom=`head -1 $filexyz | awk '{print $1}'`
zmin=`tail -$natom $filexyz | sort -nk4 | head -1 | awk '{print $4}'`
zmax=`tail -$natom $filexyz | sort -nk4 | tail -1 | awk '{print $4}'`
echo $zmin $zmax
zcenter=`echo "scale=9; 0.5*($zmin + $zmax)" | bc -l`
echo $zcenter

head -2 $filexyz > NEW_OUTPUT.xyz
tail -$natom $filexyz | awk '{printf "%4s %15.9f %15.9f %15.9f\n", $1, $2, $3, $4-1.0*'"$zcenter"'}' >> NEW_OUTPUT.xyz
