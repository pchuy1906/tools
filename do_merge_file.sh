echo "usage: \"do_merge_file.sh arg1_filename\""
if [ "$#" -eq 0 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

OUTCAR="$1"
nfile=`ls */$1 | wc -l`

rm -rf ${OUTCAR}

for ifile in $(seq 1 $nfile); do
  echo run-${ifile}
  cd run-${ifile}/
    cat ${OUTCAR} >> ../${OUTCAR}
  cd ..
done

cp run-1/POSCAR .

