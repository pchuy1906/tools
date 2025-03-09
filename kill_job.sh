echo "usage: \"./kill_job.sh arg1_fileQSTAT\""

if [ "$#" -eq 0 ]; then
  echo "No filename found. EXIT PROGRAM NOW "
  exit 1
fi


while read line
do
JOBID=`echo $line | awk '{print $1}'`
scancel $JOBID
done < $1
