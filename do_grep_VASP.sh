grep -e "T="  conf0.out  | awk '{print $3}' >  Temperature.dat
grep -e "T="  conf0.out  | awk '{print $7}' >  Energy.dat 
grep -e "external pressure"  OUTCAR  | awk '{print $4}'  >  Pressure.dat

grep -e "RMM: 200" conf0.out
