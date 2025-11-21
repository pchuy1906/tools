echo
echo "usage: \"do_script.sh  path1 path2 \""
echo

if [ "$#" -ne 2 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

path1=$1
path2=$2

cwd=`pwd`
cd $path1
allfiles=`ls -1vd bond*.dat`
cd $cwd


for file in $allfiles ; do
    echo $file
    xmgrace $path1/$file $path2/$file
done
