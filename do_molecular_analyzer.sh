echo "usage: ./script.sh  file_molecule_list  file_output_molanal"

if [ "$#" -ne 2 ]; then
    echo "error. exit now !!!!"
    exit 1
fi

file_molecule_list=$1
file_output_molanal=$2
echo "file_molecule_list is: " ${file_molecule_list}
echo "file_output_molanal: " ${file_output_molanal}
echo 

keys="Concentration history for"

neach=`awk '!/'"$keys"'/{count++}/'"$keys"'/{print count; count = 0}' ${file_output_molanal} | tail -1 | awk '{print $1-1}'`
echo "neach=" $neach

ncount=0
while read line; do
    sym_curr=`echo $keys $line`
    ncount=$(($ncount+1))
    echo ${ncount} $line
    grep -A${neach} "${sym_curr}"  ${file_output_molanal} | tail -n +7  > molecule_${ncount}.dat
done < ${file_molecule_list}
