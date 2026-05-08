nnext=`ls -d run-* | wc -l`
ncurrent=$(($nnext-1))

cd run-$ncurrent
    pwd
    ~/tools/DFTB/do_grep_DFTB+.sh  DFTB_output 0
cd ..


cd run-$nnext
    pwd
    n0=`head -2 ini_restart.xyz | tail -1 | awk '{print $NF}'`
    ~/tools/DFTB/do_grep_DFTB+.sh  DFTB_output $n0
cd ..


xmgrace run-$ncurrent/Temperature.dat  run-$nnext/Temperature.dat
