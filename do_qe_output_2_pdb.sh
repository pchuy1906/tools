echo "usage: \"do_script.sh arg1_filename arg2_option(1-Crystal 2-Angstrom)\""
if [ "$#" -ne 2 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

file=$1
option=$2

nat=`grep -e "number of atoms" ${file}  | tail -1 | awk '{print $NF}'`

n1=`awk '/CELL_P/{ print NR}' ${file} | tail -1`
awk "NR >= $n1"  ${file} > aa

grep -A3 "CELL_PARA" aa | tail -3 > box

if [ "$option" = "1" ]; then
    grep -A${nat} "ATOMIC_POSITIONS" aa | tail -${nat} > xyz
    python  ~/tools/others/Python_tools/C2A.py > AAA
    mv Axyz.out xyz
elif [ "$option" = "2" ]; then
    grep -A${nat} "ATOMIC_POSITIONS" aa | tail -${nat} > xyz
else
    echo "unknown option2. EXIT THE PROGRAM"
    exit 1
fi
python ~/tools/others/Python_tools/xyz2pdb.py 

rm -rf aa box xyz AAA
