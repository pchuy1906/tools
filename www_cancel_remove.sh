cat VASP_*/JOB* > AAA
sed -i 's/Submitted batch job//g' AAA
~/tools/others/do_cancel.sh  AAA 
rm -rf AAA VASP_*
