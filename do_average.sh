echo "usage: \"do_script.sh arg1_filename  arg2_column_number \""
if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filename=$1
cnum=$2

cmd="awk '{ total += $"${cnum}"; count++ } END { print total/count }' ${filename}"

echo $cmd > TMP
chmod +x TMP
./TMP
