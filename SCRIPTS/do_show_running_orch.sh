last_output=`ls -ltr LMP_VASP_calcs/*/*/*/vasp.out LMP_VASP_calcs/*/*/*/lammps.out | tail -1 | awk '{print $NF}'`
echo ${last_output}
tail -f ${last_output}
