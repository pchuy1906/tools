echo
echo "usage: \"do_script.sh  file.txt \""
echo

if [ "$#" -ne 1 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi


while read line ; do
    JOBID=`echo $line | awk '{print $1}'`
    echo $JOBID
    scancel $JOBID
    bkill $JOBID

    JOBID=`echo $line | awk '{print $4}'`
    echo $JOBID
    scancel $JOBID
    bkill $JOBID

done < $1
