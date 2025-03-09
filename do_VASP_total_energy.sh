nrun=`ls -d run-* | wc | awk '{print $1}'`
#echo ${nrun}

output="run-${nrun}/OUTCAR"
#grep -e "free  energy   TOTEN"  ${output}
grep -e "free  energy   TOTEN"  ${output} | tail -1 | awk '{print $5}'

