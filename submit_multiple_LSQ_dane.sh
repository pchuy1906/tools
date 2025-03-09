echo
echo "usage: \"do_script.sh  job2.01.sh \""
echo

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

fjob="$1"

for bank in pbronze mxene iap enrgmat ; do
    file="job_${bank}.sh"
    cp $fjob $file
            sed -i 's/.*#SBATCH -A.*/#SBATCH -A '"$bank"'/' $file
            sed -i 's/.*nMPI=.*/nMPI=224/'                  $file
            sed -i 's/.*#SBATCH -N.*/#SBATCH -N 2/'         $file
            sed -i 's/.*nnodes=.*/nnodes=2/'                $file
    sbatch $file &> JOB_ID_${bank}
done

