echo
echo "usage: \"do_script.sh  arg1_POSCAR/CONTCAR  \""
echo
echo "calculate the Zero point energy"

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1

echo
echo "The zero point energy (meV) is:"
awk '{ total += $(NF-1); count++ } END { print total/count }'  $file | awk '{print $1*0.5}'
