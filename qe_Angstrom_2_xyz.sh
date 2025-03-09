echo "usage: \"do_script.sh    qe.in \""
if [ "$#" -ne 2 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

file=$1
fileout="file.xyz"

natom=`grep -e "nat" $file | awk '{print $3}'`

echo ${natom} > ${fileout}
if [ "$2" == "box33" ]; then
    grep -A3 "CELL_PARAMETERS" $file | tail -3 >> ${fileout}
elif [ "$2" == "box9" ]; then
    grep -A3 "CELL_PARAMETERS" $file | tail -3 | xargs >> ${fileout}
fi
grep -A${natom} "ATOMIC_POSITIONS" $file | tail -${natom} >> ${fileout}


