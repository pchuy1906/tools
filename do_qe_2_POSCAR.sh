echo "usage: \"do_script.sh arg1_filename  arg2 (1-Crystal 2-Angstrom)\""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1
ioption=$2

nat=`grep nat ${file} | awk '{print $3}'`

grep -A3 "CELL_PARA" ${file} | tail -3 > box
grep -A${nat} "ATOMIC_POSITIONS" ${file} | tail -${nat} > xyz

if [ $ioption -eq 1 ] ; then
    python ~/tools/others/Python_tools/C2goodA.py
    mv Axyz.out xyz
fi

curr_dir=`pwd`
work_dir="$HOME/tools/VASP/box_and_xyz_2_POSCAR/"
cp box  xyz ${work_dir}

cd ${work_dir}
./run.sh
cp POSCAR ${curr_dir}
cd ${curr_dir}

