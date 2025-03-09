echo
echo "usage: \"do_script.sh   bank \""
echo

if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

grep Submitted JOB_ID* | awk '{print $NF}' > AAA
~/tools/others/do_cancel.sh  AAA 
rm -rf AAA

bank="$1"

sed -i 's/.*#SBATCH -A.*/#SBATCH -A '"$bank"'/' job*.sh

./job_all.sh
