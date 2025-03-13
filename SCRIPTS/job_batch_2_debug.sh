echo
echo "usage: \"do_script.sh  job.sh  JOB_ID  \""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file="$1"
fJOBID="$2"

JOBID=`tail -1 ${fJOBID} | awk '{print $4}'` ; scancel  $JOBID ; sed -i 's/pbatch/pdebug/g' $file ; sbatch  $file  &> $fJOBID
