echo "usage: \"do_script.sh    file.xyz\""
if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi


file="$1"
nall=`wc $file | awk '{print $1}'`
natom=`head -1 $file | awk '{print $1}'`
neach=$(($natom+2))
nstruc=$(($nall/$neach))
echo "# configurations is " ${nstruc}
ntake=$(($nstruc/2))

nbegin=$(($neach*($ntake-1)+1))
nend=$(($neach*($ntake-1)+$neach))

awk 'NR >= '"$nbegin"' && NR <= '"$nend"'' $file > middle.xyz
