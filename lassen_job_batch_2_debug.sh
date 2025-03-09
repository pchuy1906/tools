echo
echo "usage: \"do_script.sh  job.sh  JOB_ID  \""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file="$1"
fJOBID="$2"

JOBID=`tail -1 ${fJOBID} | awk '{print $2}' | sed 's/<//g' | sed 's/>//g'` 
echo $JOBID
bkill  $JOBID ; sed -i 's/pbatch/pdebug/g' $file ; bsub  $file  &> JOB_ID
