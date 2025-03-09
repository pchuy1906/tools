echo "usage: \"do_script.sh     file.xyz  1 (1:sort x;   2:sort y;   3:sort z)\""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
xyz=$2

fileout="sort_${xyz}_${file}"
head -2 $file > $fileout
natom=`head -1 $file | awk '{print $1}'`
if [ $xyz -eq 1 ] ; then
    tail -$natom $file | sort -nk2 >> $fileout
elif [ $xyz -eq 2 ] ; then
    tail -$natom $file | sort -nk3 >> $fileout
elif [ $xyz -eq 3 ] ; then
    tail -$natom $file | sort -nk4 >> $fileout
else
    echo "unknown options 2. EXIT THE PROGRAM"
    exit 1
fi
