while read line
do
  JOBID=`echo $line | awk '{print $1}'`
  echo $JOBID
  find -name "*$JOBID*"
done < $1
