echo "usage: \"do_script.sh    file.xyz    option-cell(1-cell_3, 2-cell_9, 3-non-orthor)  1.5  C H N O  \""
if [ "$#" -ne 5 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filexyz="$1"
option="$2"
rcut="$3"
type1="$4"
type2="$5"

RDFcode="$HOME/tools/bonds_count/"
cwd=`pwd`
cd $RDFcode
pwd
./do_compile.sh
cd $cwd

$RDFcode/a.out ${type1} ${type2} $rcut $filexyz $option &> _OUTPUT_RDF
mv bond.dat bond_${type1}_${type2}.dat

