echo "usage: \"do_script.sh  density(g/cm^3) \""
if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

nden=$1

prefix="splitted_denz"
file1=`ls ${prefix}* | head -1`
awk '{print $1, '"$nden"'}' $file1 > _TMP_nden.dat

xmgrace _TMP_nden.dat ${prefix}* 
