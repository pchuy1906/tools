for file in $(ls -1vd restart.*.bak) ; do
    file2="${file%.bak}"
    echo $file $file2
    mv $file $file2
done
