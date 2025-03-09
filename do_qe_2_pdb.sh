echo "usage: \"do_script.sh arg1_filename\""
if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1

nat=`grep nat ${file} | awk '{print $3}'`

grep -A3 "CELL_PARA" ${file} | tail -3 > box
grep -A${nat} "ATOMIC_POSITIONS" ${file} | tail -${nat} > xyz

python ~/tools/others/Python_tools/xyz2pdb.py 
