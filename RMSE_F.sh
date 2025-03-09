echo "********************************************"
echo "FORCE"
paste DFT.FORCE_ENERGY_STRESS  force.txt | grep FORCE | awk '{print $1, $3}' > COMPARE.dat
python   ~/tools/others/RMS.py 

echo "********************************************"
echo "ENERGY"
paste DFT.FORCE_ENERGY_STRESS  force.txt | grep ENERGY | awk '{print $1, $3}' > COMPARE.dat
paste COMPARE.dat ref.natom | awk '{printf "%15.9f %15.9f\n", $1/$3, $2/$3}' > tmp_file
mv tmp_file  COMPARE.dat
python   ~/tools/others/RMS.py

echo "********************************************"
echo "STRESS"
paste DFT.FORCE_ENERGY_STRESS  force.txt | grep DIA_STRESS | awk '{print $1, $3*6.9479}' > COMPARE.dat
#paste DFT.FORCE_ENERGY_STRESS  force.txt | grep STRESS | awk '{print $1, $3*6.9479}' > COMPARE.dat
python   ~/tools/others/RMS.py

