cwd=`pwd`
for file in `find -name fm_setup.in` ; do
    fold=${file/fm_setup.in/}
    echo  $fold
    cd $fold
        mv frames.all.log _frames.all.log
        rm -rf A.* b-label* dim.* restart* frame* Ax.0*.txt x.0*.txt
    cd $cwd
done
