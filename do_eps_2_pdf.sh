for file in $(ls *.eps) ; do
    echo $file
    epstopdf $file
done
