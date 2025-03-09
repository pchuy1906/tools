while read line; do
  file=`echo $line`
  echo $file
  cp do_2_collect.sh $file
done < $1
