echo "usage: \"do_split_file.sh arg1_filename arg2_nline_each arg3_100\""
if [ "$#" -ne 3 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filename="$1"
nline=$2
nrange=$3

nframe=`wc ${filename} | awk '{print $1/'"$nline"'}'`

echo
echo "Spliting file $1 which has $nframe frames"

nfolder=$(($nframe/$nrange))

echo $nfolder

for ifile in $(seq 1 $nfolder) ; do
    line_begin=1
    line_end=$(($nline*$nrange+$nline*$nrange*($ifile-1)))
    rm -rf range-${ifile}; mkdir range-${ifile}
    sed -n ${line_begin},${line_end}p  ${filename} > range-${ifile}/traj.xyz
    cd range-${ifile}
        cp -rf /g/g92/pham20/tools/rdf_from_xyzf_general/   .
        cd rdf_from_xyzf_general/
#            sed -i 's/512/272/g' ./do_run.sh
            ./do_run.sh
        cd ..
    cd ..
    echo $ifile ${line_begin} ${line_end}
done
