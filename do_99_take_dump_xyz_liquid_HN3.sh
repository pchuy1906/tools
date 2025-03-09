echo "usage: \"do_script.sh  fileDump  nstep  nwrite\""
if [ "$#" -ne 3 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
nstep=$2
ndump=$3

echo "Write the dump file every $ndump steps"
if [ $(( $nstep % $ndump )) -ne 0 ] ; then
    echo "ERROR: the nstep % ${ndump}  should be 0"
    exit 1
fi

natom=`head -4 ${file} | tail -1 | awk '{print $1}'`
neach=$(($natom+9))

nstart=`head -2 ${file} | tail -1 | awk '{print $1}'`
if [ ${nstep} -lt ${nstart}  ] ; then
    echo "ERROR: nstep (${nstep}) should be larger than nstart (${nstart}) "
    exit 1
fi

nid=$((($nstep - $nstart)/$ndump))
nbegin=$(($neach*$nid+1))
nend=$(($neach*$nid+$neach))

sed -n ${nbegin},${nend}p  ${file}  > aaa_${file}


~/tools/LAMMPS/dump_2_xyz/a.out   aaa_dump_xyz_vxyz  2 N H
mv file.xyz ${nstep}.xyz
