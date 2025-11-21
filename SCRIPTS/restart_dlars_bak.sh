for file in $(ls restart* | grep -v bak) ; do
    ls ${file}.bak 
    ls $file
    mv ${file}.bak $file
done 
