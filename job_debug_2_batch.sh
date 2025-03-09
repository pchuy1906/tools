echo
echo "usage: \"do_script.sh  job.sh  JOB_ID\""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file="$1"
fileJOB="$2"

JOBID=`tail -1 ${fileJOB} | awk '{print $4}'` ; scancel  $JOBID ; sed -i 's/pdebug/pbatch/g' $file ; sbatch  $file  &> JOB_ID
