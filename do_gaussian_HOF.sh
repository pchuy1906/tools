H2kcalmol=627.503

echo "need files: AE_QM.dat  enthalpy_0K.dat enthalpy_difference.dat"


echo "usage: \"do_script.sh     output \""
if [ "$#" -ne 1 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

output=$1

echo
grep -e "NAtoms"  ${output} | head -1
natoms=`grep -e "NAtoms"  ${output} | head -1 | awk '{print $2}'`
echo "the number of atoms" ${natoms}

echo
grep -e "SCF Done"  ${output}
grep -e "SCF Done"  ${output} | awk '{printf "%15.9f", $5*'"${H2kcalmol}"'}' > Energy.dat

echo
grep -A${natoms} "Multiplicity"  ${output} | tail -${natoms} | awk '{print $1}' | sort | uniq -c
grep -A${natoms} "Multiplicity"  ${output} | tail -${natoms} | awk '{print $1}' | sort | uniq -c | awk '{print $1}' > ntypesCHNO.dat

E0=`paste AE_QM.dat  ntypesCHNO.dat | awk '{printf "%15.9f\n", $1*$2}' | awk '{s+=$1} END {printf "%15.9f", s}'`


echo
grep -e "Zero-point correction=" ${output}
ZPE=`grep -e "Zero-point correction=" ${output} | awk '{printf "%15.9f", $3*'"${H2kcalmol}"'}'`
echo "ZERO-POINT ENERGY = " ${ZPE}



echo
grep -e "Sum of electronic and zero-point Energies="  ${output} 
E1=`grep -e "Sum of electronic and zero-point Energies="  ${output} | awk '{printf "%15.9f", $NF*'"${H2kcalmol}"'}'`

echo
atomization_energy=`echo "scale=5; $E0-1.0*${E1}" | bc -l`
echo "ATOMIZATION ENERGY = ${atomization_energy}"

echo
E2=`paste enthalpy_0K.dat ntypesCHNO.dat | awk '{printf "%15.9f\n", $1*$2}' | awk '{s+=$1} END {printf "%15.9f", s}'`
enthalpy_at_0K=`echo "scale=5; $E2-1.0*${atomization_energy}" | bc -l`
echo "ENTHALPY 0K = ${enthalpy_at_0K}"


echo 
E3=`paste enthalpy_difference.dat ntypesCHNO.dat | awk '{printf "%15.9f\n", $1*$2}' | awk '{s+=$1} END {printf "%15.9f", s}'`

echo
grep -e "Thermal correction to Enthalpy=" ${output}
E4=`grep -e "Thermal correction to Enthalpy=" ${output} | awk '{printf "%15.9f", $NF*'"${H2kcalmol}"'}'`
echo "Thermal correction to Enthalpy (kcal/mol) = $E4" 

echo
enthalpy_at_298K=`echo "scale=5; ${enthalpy_at_0K}+1.0*${E4}-1.0*${E3}-${ZPE}" | bc -l`
echo "ENTHALPY 298K = ${enthalpy_at_298K}"



