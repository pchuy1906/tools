echo "usage: \"do_script.sh    arg1_fileXYZ  arg2_structure_ID \""
if [ "$#" -ne 2 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

file=$1
ID=$2

natom=`head -1 ${file} | awk '{print $1}'`
neach=$(($natom+2))

nbegin=$(($neach*($ID-1)+1))
nend=$(($neach*$ID))

sed -n ${nbegin},${nend}p  $file  > ${ID}.xyz

echo "Unit-cell:"
sed -n ${ID},${ID}p AAA_Lxyz.dat


echo "Density:"
sed -n ${ID},${ID}p AAA_Den.dat
