curr=`pwd`

Eall=`~/tools/others/do_VASP_total_energy.sh`
echo "Total energy (AA) is " $Eall

cd /p/lscratchh/pham20/Ferrihydrite_Magnetite/isolated_H2
E_H2=`~/tools/others/do_VASP_total_energy.sh`
echo "Total energy (H2) is " ${E_H2}
cd ${curr}


cd /p/lscratchh/pham20/Ferrihydrite_Magnetite/isolated_O2
E_O2=`~/tools/others/do_VASP_total_energy.sh`
echo "Total energy (O2) is " ${E_O2}
cd ${curr}


head -7 run-1/POSCAR | tail -2
