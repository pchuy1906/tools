echo
echo "usage: \"do_script.sh  file.xyz \""
echo

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filexyz=$1

fileout="AAA.xyz"

natom=`head -1 $filexyz | awk '{print $1}'`
tail -$natom $filexyz > TMP_XYZ

Xmin=`sort -nk2 TMP_XYZ | head -1 | awk '{print $2}'`
Ymin=`sort -nk3 TMP_XYZ | head -1 | awk '{print $3}'`
Zmin=`sort -nk4 TMP_XYZ | head -1 | awk '{print $4}'`

echo $natom > $fileout
echo "30 20 20" >> $fileout
awk '{printf "%4s %15.9f %15.9f %15.9f\n", $1, $2-1.0*'"$Xmin"'+5.0, $3-1.0*'"$Ymin"'+5.0, $4-1.0*'"$Zmin"'+5.0}' TMP_XYZ  >> $fileout

rm -rf TMP_XYZ
