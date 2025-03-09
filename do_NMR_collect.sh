echo "usage: \"do_script.sh arg1_filename  arg2_atom_name agr3_ref\""
if [ "$#" -ne 3 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
atomname=$2
ref=$3

grep -e "Total sigma" $file | grep -e " $atomname " > TMP_1
awk '{printf "%15.9f\n", '"$ref"'-$NF}'  TMP_1 > sigma.dat
