echo "usage: \"do_script.sh arg1_filename arg2_option(1-Angstrom  2-Crystal)\""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
option=$2

nat=`grep nat ${file} | awk '{print $3}'`

grep -A3 "CELL_PARA" ${file} | tail -3 > box
grep -A${nat} "ATOMIC_POSITIONS" ${file} | tail -${nat} > xyz


if [ "$option" = "2" ]; then
    python  ~/tools/others/Python_tools/C2A.py > AAA
    mv Axyz.xyz xyz
elif [ "$option" = "1" ]; then
    echo "good to go"
    #python ~/tools/others/Python_tools/A2goodA.py &> AAA
    #mv Axyz.out xyz
else
    echo "unknown option2. EXIT THE PROGRAM"
    exit 1
fi

aa1=`sed -n 1,1p box`
aa2=`sed -n 2,2p box`
aa3=`sed -n 3,3p box`

echo ${nat} > file.xyz
echo $aa1 $aa2 $aa3 >> file.xyz
tail -$nat xyz >> file.xyz
rm -rf box xyz AAA
