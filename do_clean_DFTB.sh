cwd=`pwd`
for file in `find -name dftb_in.hsd` ; do
    fold=${file/dftb_in.hsd/}
    echo  $fold
    cd $fold
        rm -rf *.skf Energy.dat Pressure.dat Temperature.dat autotest.tag band.out charges.bin detailed.out dftb_pin.hsd last_velocity md.out
    cd $cwd
done
