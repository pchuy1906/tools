for fold in $(ls -d */); do
    echo $fold
    cd $fold
        $HOME/tools/others/do_number_of_file.sh
    cd ..
    echo
done
