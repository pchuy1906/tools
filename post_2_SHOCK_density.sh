prefix="splitted_denz"
rm -rf _TMP.dat
for file in $(ls ${prefix}*) ; do
    cat $file >> _TMP.dat
done

min_denz=1.2
max_denz=2.5

file="_TMP.dat"

echo
echo "Calculate the distribution:"
cmd="python $HOME/tools/others/calc_distribution.py --file_input ${file} --nbins 150 --col_use 1 --min_val ${min_denz} --max_val ${max_denz} "
echo $cmd
$cmd &> _OUTPUT_1
echo

echo "Fit the histogram data"
cmd="python $HOME/tools/others/fit_python_2_gaussian.py  --file_input hist.dat   --xcol 0 --ycol 1 "
echo $cmd
$cmd &> _OUTPUT_2
echo
