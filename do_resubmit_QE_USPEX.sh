nrun=`ls output* | wc | awk '{print $1}'`
echo ${nrun}

mv output output_${nrun}
mv qe.in qe.in_${nrun}

nhead=`awk '/CELL_P/{ print NR-1; exit }' qe.in_${nrun}`
head -${nhead} qe.in_${nrun} > qe.in


n1=`awk '/CELL_P/{ print NR }' output_${nrun} | tail -1`
nat=`grep nat qe.in_${nrun}  | tail -1 | awk '{print $3}'`
n2=$(($n1+$nat+5))

sed -n ${n1},${n2}p  output_${nrun} >> qe.in

tail -2 qe.in_1  >> qe.in
sed -n ${n1},${n2}p  output_${nrun} > aaa

msub myrun
rm -rf aaa
