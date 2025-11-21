echo
echo "usage: \"do_script.sh  pbronze job.sh  JOB_ID  \""
echo

if [ "$#" -ne 3 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

bank="$1"
file="$2"
fJOBID="$3"

JOBID=`tail -1 ${fJOBID} | awk '{print $4}'` ; scancel  $JOBID ; sed -i '/#SBATCH -A/c\#SBATCH -A '"$bank"''  $file ; sbatch  $file  &> $fJOBID
