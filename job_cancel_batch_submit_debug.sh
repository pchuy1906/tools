grep Submitted JOB_ID* | awk '{print $NF}' > AAA
~/tools/others/do_cancel.sh  AAA 
rm -rf AAA

sed -i 's/pbatch/pdebug/g' job*.sh
./job_all.sh
