echo
echo "usage: \"do_script.sh  job.sh  JOB_ID\""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file="$1"
filejob="$2"

JOBID=`tail -1 $filejob | awk '{print $4}'` ; scancel  $JOBID ; sed -i 's/pbronze/ams/g' $file ; sbatch  $file  &> JOB_ID
