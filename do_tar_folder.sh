for fold in $(ls -d */); do
    goodname="${fold///}"
    echo $goodname
    tar -cvzf ${goodname}.tgz  $goodname &> TMP_${goodname}
done
