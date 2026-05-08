# NOTE: LIFETIMES ARE IN THE *find_molecs.out FILE

fileinput="file.gen"
FILE_INPUT_PATH="../"
MOLANAL_PATH="/g/g92/pham20/codes/molanal_2018_07_09/trunk/"

${MOLANAL_PATH}/molanal.new    ${FILE_INPUT_PATH}/${fileinput}  > ${fileinput}-molanal.out

${MOLANAL_PATH}/findmolecules.pl  ${fileinput}-molanal.out  >  ${fileinput}-find_molecs.out
