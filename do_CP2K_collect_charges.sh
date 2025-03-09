echo
echo "usage: \"do_script.sh  cp2k.out natom\""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT PROGRAM"
  exit 1
fi

fileinput="$1"
natom="$2"

ntail=$(($natom+2))

grep -A$ntail "Hirshfeld Charges" $fileinput | tail -$natom | awk '{print $NF}' > Hirshfeld.dat
grep -A$ntail "Mulliken" $fileinput | tail -$natom | awk '{print $NF}' > Mulliken.dat

