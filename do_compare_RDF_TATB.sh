echo "usage: \"do_script.sh     high  C O  \""
if [ "$#" -ne 3 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

den="$1"
atom1="$2"
atom2="$3"

if [ "$den" == "high" ]; then
    path="/p/lustre2/pham20/doing/TATB/99-MD-SCAN/MD-2500-3000K-1/1-8/"
else
    path="/p/lustre2/pham20/doing/TATB/99-MD-SCAN/MD-1000-3000K/1-8/"
fi


xmgrace ~/tools/XMGRACE/tmp_RDF.agr  ${path}/gr_${atom1}_${atom2}.dat gr_${atom1}_${atom2}.dat
