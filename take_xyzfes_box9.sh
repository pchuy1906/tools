OUTCAR="OUTCAR"
POSCAR="POSCAR"
OUTPUT="input.xyzfes"
OUTPUTF="DFT.FORCE_ENERGY_STRESS"

eV2Ha=0.0367493
Ang2Bohr=1.889725989
eV2kcal=23.061
kbar2GPa=0.1

echo "usage: \"./take_xyz_f.sh\" or \"./take_xyz_f.sh arg1_nneed arg2_nneglect\""

grep -e "free  energy   TOTEN"  ${OUTCAR} | awk '{printf "%15.9f\n", $5*'"$eV2kcal"'}' > ENERGIES.dat

grep -e "in kB"  ${OUTCAR} | awk '{printf "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", $3*'"$kbar2GPa"', $4*'"$kbar2GPa"', $5*'"$kbar2GPa"', $6*'"$kbar2GPa"', $7*'"$kbar2GPa"', $8*'"$kbar2GPa"'  }' > STRESS_TENSOR.dat

if [ -f ${OUTCAR} ]; then
  echo "Found the file ${OUTCAR}"
else
  echo "File ${OUTCAR} is not found. EXIT THE PROGRAM"
  exit 1
fi
grep -e "NIONS = "  ${OUTCAR}
natoms=`grep -e "NIONS = "  ${OUTCAR} | tail -1 | awk '{print $NF}'`
echo "The number of atoms is: " ${natoms}
echo

nprint=$(($natoms+1))
echo "--" > file.tmp1
grep -A${nprint} "TOTAL-FORCE" ${OUTCAR} >> file.tmp1
echo "export POSITIONS and FORCES --> DONE"
echo

echo "This script works only if the cell does NOT change. Modify if needed!!"
grep -A3 "direct lattice vectors" ${OUTCAR} | tail -3 | awk '{print $1, $2, $3}' | xargs > boxABC
echo "NON_ORTHO" > NON_ORTHO

nconfig=`grep -e "TOTAL-FORCE" ${OUTCAR} | wc -l`
echo "the number of configurations is " ${nconfig}
echo

if [ "$#" -eq 0 ]; then
  echo "export ALL ${nconfig} configurations "
  NSTRUC_NEED=`seq 1 $nconfig`
  #echo ${NSTRUC_NEED}
elif [ "$#" -eq 2 ]; then
  nneglect=$2
  nneed=$1
  ndel=$((($nconfig-$nneglect-1)/($nneed-1)))
  echo "export ${nneed} configurations, ignore the first ${nneglect} configurations"
  nstart=$(($nneglect+1))
  NSTRUC_NEED=`seq ${nstart} $ndel ${nconfig}`
  echo ${NSTRUC_NEED}
else
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi
echo

neach=$(($natoms+3))

if [ -f ${POSCAR} ]; then
  echo "Found the file ${POSCAR}"
else
  echo "File ${POSCAR} is not found. EXIT THE PROGRAM"
  exit 1
fi
atomSym=`sed -n 6,6p ${POSCAR}`
atomNum=`sed -n 7,7p ${POSCAR}`
numElement=`sed -n 7,7p ${POSCAR} | awk '{print NF}'`
echo "System details:" ${numElement} "elements"
echo ${atomSym}
echo ${numElement} ${atomSym} > atomSym.dat
echo ${atomNum}
echo

declare -A SymArray
ncount=0
for tmp in ${atomSym}; do
  ncount=$(($ncount+1))
  SymArray[$ncount]=$tmp
done

declare -A NumArray
ncount=0
for tmp in ${atomNum}; do
  ncount=$(($ncount+1))
  NumArray[$ncount]=$tmp
done

rm -rf ATOM.NAME
for itype in $(seq 1 $numElement); do
  echo  ${SymArray[$itype]} ${NumArray[$itype]}
  for irun in $(seq 1 ${NumArray[$itype]}); do
    echo ${SymArray[$itype]} >> ATOM.NAME
  done
done
echo

rm -rf ${OUTPUT} ${OUTPUTF}
rm -rf ID_CHECK DFT.STRESS_TENSOR

for iconfig in ${NSTRUC_NEED}; do
  echo ${natoms} >> ${OUTPUT}
  sed -n ${iconfig},${iconfig}p  ENERGIES.dat >  thisENERGY

  sed -n ${iconfig},${iconfig}p  STRESS_TENSOR.dat > thisSTRESS

  paste NON_ORTHO boxABC thisSTRESS thisENERGY atomSym.dat >> ${OUTPUT}
  #paste NON_ORTHO boxABC thisSTRESS thisENERGY  >> ${OUTPUT}

  ibegin=$(($neach*($iconfig-1)+4))
  iend=$(($neach*$iconfig))
  echo ${iconfig} ${ibegin} ${iend} >> ID_CHECK
  sed -n ${ibegin},${iend}p file.tmp1 > ATOM.XYZ
  # VASP has force (eV/A), while ChIMES used (Ha/Bohr) for input; and (kcal/mol/A) for output
  paste ATOM.NAME ATOM.XYZ | awk '{printf "%4s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", $1, $2, $3, $4, $5*'"$eV2Ha"'/'"$Ang2Bohr"', $6*'"$eV2Ha"'/'"$Ang2Bohr"', $7*'"$eV2Ha"'/'"$Ang2Bohr"'}' >> ${OUTPUT}
  while read line; do
    echo $line | awk '{printf "%15.9f %15.9f %15.9f\n", $4*'"$eV2kcal"', $5*'"$eV2kcal"', $6*'"$eV2kcal"'}' | fmt -1 | awk '{print $1, "FORCE"}' >> ${OUTPUTF}
  done < ATOM.XYZ

while IFS=" " read -r sig_xx sig_yy sig_zz sig_xy sig_yz sig_zx remainder
do
    echo $sig_xx  "DIA_STRESS" >> ${OUTPUTF} 
    echo $sig_xy  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_zx  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_xy  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_yy  "DIA_STRESS" >> ${OUTPUTF} 
    echo $sig_yz  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_zx  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_yz  "OFF_STRESS" >> ${OUTPUTF} 
    echo $sig_zz  "DIA_STRESS" >> ${OUTPUTF} 
done < "thisSTRESS"


  cat thisENERGY | awk '{print $1, "ENERGY"}' >> ${OUTPUTF}
  cat thisENERGY | awk '{print $1, "ENERGY"}' >> ${OUTPUTF}
  cat thisENERGY | awk '{print $1, "ENERGY"}' >> ${OUTPUTF}

done

rm -rf ATOM.NAME  ATOM.XYZ   boxABC  file.tmp1 thisENERGY  ENERGIES.dat
