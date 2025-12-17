echo
echo "usage: \"do_script.sh  fileFIT filePAIR \""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

fileFIT=$1
filePAIR=$2
echo $fileFIT
grep pair_coeff $fileFIT > fitted_parameters.dat


awk '
    # Read file2 and store its lines using the first 4 columns as the key
    FNR==NR {
        key = $1 FS $2 FS $3 FS $4
        lines[key] = $0
        next
    }
    # For each line in file1, check if a matching key exists in file2
    {
        key = $1 FS $2 FS $3 FS $4
        if (key in lines) {
            print lines[key]
        } else {
            print $0
        }
    }
' fitted_parameters.dat $filePAIR > NEW_PAIR
