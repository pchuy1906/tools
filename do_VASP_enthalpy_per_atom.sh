eV2kcalmol="23.0609"

nrun=`ls -d run-* | wc | awk '{print $1}'`

#echo ${nrun}
output="run-${nrun}/OUTCAR"
NION=`grep -e "NIONS"         ${output} | tail -1 | awk '{print $NF}'`
ETOT=`grep -e "free  energy"  ${output} | tail -1 | awk '{printf "%15.9f\n", $5*'"${eV2kcalmol}"'/'"${NION}"'}'`
P=`grep -e "external pressure"   ${output} | tail -1 | awk '{printf "%15.9f\n", $(NF-1)*0.1}'`

echo $P $ETOT
