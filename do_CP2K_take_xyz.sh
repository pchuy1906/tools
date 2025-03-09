echo "usage: \"do_script.sh    DFT.xyz  cp2k.out cp2k.inp\""
if [ "$#" -ne 3 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

fileoutxyz=$1
fileout=$2
fileinp=$3

output="__file.xyz"

natom=`head -1 ${fileoutxyz} | awk '{print $1}'`
echo "The number of atom is:" $natom

ntail=$(($natom+2))
tail -${ntail} ${fileoutxyz} > ${output}

grep -e "CELL" ${fileout} > __TMP_CELL
a=`grep -e " a " __TMP_CELL | tail -1 | awk '{print $NF}'`
b=`grep -e " b " __TMP_CELL | tail -1 | awk '{print $NF}'`
c=`grep -e " c " __TMP_CELL | tail -1 | awk '{print $NF}'`
newABC="ABC $a $b $c"
alp=`grep -e " alpha " __TMP_CELL | tail -1 | awk '{print $NF}'`
bet=`grep -e " beta "  __TMP_CELL | tail -1 | awk '{print $NF}'`
gam=`grep -e " gamma " __TMP_CELL | tail -1 | awk '{print $NF}'`
newANGLEs="ALPHA_BETA_GAMMA $alp $bet $gam"



cp ${fileinp} __${fileinp}
sed -i '/ABC/c\'"${newABC}"'' __${fileinp}
sed -i '/ALPHA_BETA_GAMMA/c\'"${newANGLEs}"'' __${fileinp}




