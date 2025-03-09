$HOME/tools/others/do_collect_LAMMPS_MSST.sh

Volume=`python  ~/tools/others/calculate_average.py --file_input  AAA_Vol.dat | grep -e "average=" | awk '{print $NF}'`
Temp=`python  ~/tools/others/calculate_average.py --file_input  AAA_Temp.dat | grep -e "average=" | awk '{print $NF}'`
Press=`python  ~/tools/others/calculate_average.py --file_input  AAA_Press.dat | grep -e "average=" | awk '{print $NF}'`
Etot=`python  ~/tools/others/calculate_average.py --file_input  AAA_Etot.dat | grep -e "average=" | awk '{print $NF}'`

file="average_data.dat"
rm -rf ${file}
echo $Volume >  ${file}
echo $Temp   >> ${file}
echo $Press  >> ${file}
echo $Etot   >> ${file}


