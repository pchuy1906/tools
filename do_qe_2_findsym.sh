echo "usage: \"do_script.sh  arg1_filename  arg2 (1-Crystal 2-Angstrom)\""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
ioption=$2

echo "Input is QE input file: $file"


nat=`grep nat ${file} | awk '{print $3}'`

grep -A3 "CELL_PARA" ${file} | tail -3 > box
grep -A${nat} "ATOMIC_POSITIONS" ${file} | tail -${nat} > xyz



echo "testing" > find_sym.in
echo "0.001 accuracy" >> find_sym.in
echo "1  form of lattice parameters: to be entered as box" >> find_sym.in
grep -A3 "CELL_PARA" ${file} | tail -3 >> find_sym.in


echo "1  form of primitive lattice vectors" >> find_sym.in
echo "1 0 0 primitive lattice vectors" >> find_sym.in
echo "0 1 0" >> find_sym.in
echo "0 0 1" >> find_sym.in
echo ${nat} " number of atoms in the primitive unit cell" >> find_sym.in

awk '{print $1}' xyz > xyz.1
aa=`~/tools/others/transpose_matrix.sh   xyz.1`
echo ${aa} " type of each atom" >> find_sym.in

if [ $ioption -eq 2 ] ; then
    python ~/tools/others/Python_tools/A2C.py > AAA
    awk '{print $2, $3, $4}' Cxyz.out >> find_sym.in
else
    awk '{print $2, $3, $4}' xyz      >> find_sym.in
fi

rm -rf box xyz AAA xyz.1 Cxyz.out

#sed -i 's/C/1/g'  find_sym.in
#sed -i 's/H/1/g'  find_sym.in
#sed -i 's/O/2/g'  find_sym.in
#sed -i 's/N/4/g'  find_sym.in
#sed -i 's/Li/1/g'  find_sym.in
#sed -i 's/Cs/2/g'  find_sym.in

sed -i 's/C/1/g'  find_sym.in
sed -i 's/H/2/g'  find_sym.in
sed -i 's/N/3/g'  find_sym.in
sed -i 's/Zn/4/g'  find_sym.in

