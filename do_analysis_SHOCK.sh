echo "usage: \"do_script.sh  (1)splitted_denz (2)min_denz (3)max_denz (4)gauwmax (5)minxgau1 (6)maxxgau1 (7)minxgau2 (8)maxxgau2 (9)minxgau3 (10)maxxgau3 \""

if [ "$#" -ne 10 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

key=$1
file="_${key}.dat"
min_denz=$2
max_denz=$3
gauwmax=$4
minxgau1="$5"
maxxgau1="$6"
minxgau2="$7"
maxxgau2="$8"
minxgau3="$9"
maxxgau3="${10}"


echo
echo "Calculate the distribution:"
cmd="python $HOME/tools/others/calc_distribution.py --file_input ${file} --nbins 150 --col_use 0 --min_val ${min_denz} --max_val ${max_denz} "
echo $cmd
$cmd &> _OUTPUT_1
echo

echo "Fit the histogram data"
cmd="python $HOME/tools/others/fit_python_3_gaussian.py  --file_input hist.dat   --xcol 0 --ycol 1 --gauwmax ${gauwmax} --minxgau1 ${minxgau1} --maxxgau1 ${maxxgau1} --minxgau2 ${minxgau2} --maxxgau2 ${maxxgau2} --minxgau3 ${minxgau3} --maxxgau3 ${maxxgau3} "
echo $cmd
$cmd &> _OUTPUT_2
echo

xgau1=`grep -e "xgau1" _OUTPUT_2 | awk '{print $NF}'`
xgau2=`grep -e "xgau2" _OUTPUT_2 | awk '{print $NF}'`
file1=`ls -1vd ${key}* | head -1`
awk '{print $1, '"$xgau1"'}' $file1 > _TMP_1.dat
awk '{print $1, '"$xgau2"'}' $file1 > _TMP_2.dat

echo
echo "us=up*f where f="
echo "scale=5; 1/(1-$xgau1/$xgau2)" | bc -l

xmgrace  ${key}*  _TMP_1.dat  _TMP_2.dat
