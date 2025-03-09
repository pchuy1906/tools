echo "usage: \"do_script.sh  arg1_POSCAR/CONTCAR  \""
echo "Support only Direct option of POSCAR/CONTCAR"

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1

echo "Input is POSCAR/CONTCAR: $file"

tmp=`head -7 ${file} | tail -1`
numtype=(${tmp})
ntype=`echo ${tmp} | wc -w`


syms=""
nat=0
for itype in $(seq 1 ${ntype}) ; do
    ind=$(($itype-1))
    ntmp=${numtype[$ind]}
    nat=$(($nat+$ntmp))
    for itmp in $(seq 1 $ntmp) ; do
	syms=`echo $syms $itype`
    done
done

echo $nat
echo $syms

head -5 ${file} | tail -3 > box
nhead=$(($nat+8))
head -${nhead}  ${file} | tail -${nat}  > xyz



echo "testing" > find_sym.in
echo "0.005 accuracy" >> find_sym.in
echo "1  form of lattice parameters: to be entered as box" >> find_sym.in
cat box >> find_sym.in


echo "1  form of primitive lattice vectors" >> find_sym.in
echo "1 0 0 primitive lattice vectors" >> find_sym.in
echo "0 1 0" >> find_sym.in
echo "0 0 1" >> find_sym.in
echo ${nat} " number of atoms in the primitive unit cell" >> find_sym.in

echo ${syms} " type of each atom" >> find_sym.in

cat xyz      >> find_sym.in

rm -rf box xyz
