rm -rf bfile

cat ~/tools/XMGRACE/tmp_xmgrace > bfile

declare -a array=("25" "38" "28" "32" "1")

ncount=0
for ddamp in 0.001 0.01 0.1 ; do

    fold="T_20000_l_${ddamp}"
    file="$fold/AAA_md_statistics.out"

    echo "READ BLOCK \"${file}\"" >> bfile
    echo "BLOCK xy \"2:5\"" >> bfile

    color=${array[$ncount]}
    echo "s${ncount} line color ${color}" >> bfile
    echo "s${ncount} legend \"SMOOTH ${ddamp}\"" >> bfile
    echo "s${ncount} linewidth 2" >> bfile


    ncount=$((ncount+1))

done


xmgrace -batch bfile
