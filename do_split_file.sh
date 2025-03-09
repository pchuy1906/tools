echo "usage: \"do_split_file.sh arg1_filename arg2_nline_each\""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filename="$1"
nline=$2

nfile=`wc $1 | awk '{print $1/'"$nline"'}'`

echo "Spliting file $1 into $nfile files"


for ifile in $(seq 1 $nfile); do
  line_begin=$((1+$nline*($ifile-1)))
  line_end=$(($nline+$nline*($ifile-1)))
  sed -n ${line_begin},${line_end}p  ${filename} > split_file.${ifile}
  echo $ifile ${line_begin} ${line_end}
done

