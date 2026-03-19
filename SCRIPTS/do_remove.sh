while read line
do
  #JOBID=`echo $line | awk '{print $1}'`
  JOBID=`echo $line`
  echo $JOBID
  rm -rf $JOBID
done < $1
