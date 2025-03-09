echo "usage: \"script.sh    file.xyz  10000 \""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filename="$1"
nline=`head -1 $filename | awk '{print $1+2}'`
nframe=`wc ${filename} | awk '{print $1/'"$nline"'}'`

nfeach="$2"

echo
echo "File $filename has $nframe frames"

nfolder=$(($nframe/$nfeach))
echo "the number of folder is $nfolder"
echo "each folder has $nfeach configs"
echo

nline_each_folder=$(( $nfeach * $nline ))


for ifile in $(seq 1 $nfolder) ; do

    line0=$(( ${nline_each_folder} * ($ifile-1) + 1 ))
    line1=$(( ${nline_each_folder} * ($ifile-1) + ${nline_each_folder} ))

    fname="RDF-${ifile}"
    rm -rf $fname; mkdir $fname

    sed -n ${line0},${line1}p  ${filename} > $fname/traj.xyz
    cd $fname
        ~/tools/others/compute_RDF_CHNO.sh  traj.xyz 1  C H N O &> OUTPUT
    cd ..
done
