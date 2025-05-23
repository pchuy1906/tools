echo "usage: \"do_script.sh    file.xyz    option-cell(1-cell_3, 2-cell_9, 3-non-orthor)  1.5  C H N O  \""
if [ "$#" -lt 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filexyz="$1"
option="$2"
rcut="$3"

ncount=0
for element in "$@" ; do
    ncount=$(($ncount + 1))
    if [ "$ncount" -gt 3 ]; then
        echo "element : $element";
        types=("${types[@]}" $element)
    fi
done


RDFcode="$HOME/tools/bonds_count/"
cwd=`pwd`
cd $RDFcode
pwd
./do_compile.sh
cd $cwd

nlen=${#types[@]}

for (( i=0; i<$nlen; i++ )) ; do
    for (( j=i; j<$nlen; j++ )) ; do
        echo $nlen ${types[i]} ${types[j]}
        $RDFcode/a.out ${types[i]} ${types[j]} $rcut $filexyz $option &> _OUTPUT_RDF
        mv bond.dat bond_${types[i]}_${types[j]}.dat
    done
done

