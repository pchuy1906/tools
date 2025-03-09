echo "usage: \"do_script.sh arg1_filename\""
if [ "$#" -ne 1 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file=$1

find -name "${file}*"  > _${file}_1
~/tools/others/do_remove.sh _${file}_1 &> _${file}_2
rm -rf AAA.dat
