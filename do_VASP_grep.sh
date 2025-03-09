nrun=`ls -d run-* | wc | awk '{print $1}'`
echo ${nrun}

output="run-${nrun}/OUTCAR"

output2="run-${nrun}/out"
#grep -e "Iteration"  ${output}
grep Pul  ${output} | tail -4
grep enthalpy   ${output} | tail -4
grep -e "reached required accuracy"  ${output}
grep -e "F="  $output2 | tail -4
