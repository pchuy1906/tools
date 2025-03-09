while read line
do
  JOBID=`echo $line`
  du -sh $JOBID
done < $1
