grep Submitted JOB_ID* | awk '{print $NF}' > AAA
~/tools/others/do_cancel.sh  AAA 
rm -rf AAA
