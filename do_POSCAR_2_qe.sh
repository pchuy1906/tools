echo
echo "usage: \"do_script.sh  arg1_POSCAR/CONTCAR  \""
echo
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
echo
echo "Atom numbers:"
echo ${tmp}

tmp=`head -6 ${file} | tail -1`
symtype=(${tmp})
echo "Atom syms:"
echo $tmp

syms=""
nat=0
echo
echo "ntype=" $ntype
echo
for itype in $(seq 1 ${ntype}) ; do
    inde=$(($itype-1))
    ntmp=${numtype[${inde}]}
    nsym=${symtype[${inde}]}
    nat=$(($nat + $ntmp))
    for itmp in $(seq 1 $ntmp) ; do
       syms=`echo $syms $nsym`
    done
done


echo
echo "The number of atom is " $nat
echo $syms > TMP_SYM
$HOME/tools/others/transpose_matrix.sh  TMP_SYM > TMP_SYM_1


nDirect=`awk '/Direct/{ print NR; exit }'  $file`
nbegin=$(($nDirect+1))
nend=$(($nbegin+$nat-1))

sed -n ${nbegin},${nend}p  $file > TMP_XYZ_2
#paste TMP_SYM_1  TMP_XYZ_2


cat > qe.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = './',
    outdir='./'
    prefix='pwscf'
    tstress = .true.
    tprnfor = .true.
 /
 &system
    ibrav = 0 
    nat = $nat
    ntyp = 1,
    ecutwfc = 25.0, ecutrho = 300.0
 /
 &electrons
    diagonalization='cg'
    conv_thr = 1.0e-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Cu 63.55 Cu.pz-d-rrkjus.UPF

K_POINTS (automatic)
 8 8 8 0 0 0
EOF


echo "CELL_PARAMETERS (angstrom)" >> qe.in
head -5 $file | tail -3 >> qe.in
echo "" >> qe.in

echo "ATOMIC_POSITIONS (crystal)" >> qe.in
paste TMP_SYM_1  TMP_XYZ_2 | awk '{printf "%4s %15.9f %15.9f %15.9f\n", $1, $2, $3, $4}' >> qe.in

rm -rf TMP_SYM_1  TMP_XYZ_2 TMP_SYM
