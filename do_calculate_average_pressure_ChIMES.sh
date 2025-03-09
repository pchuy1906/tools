file="md_statistics.out"
awk '{print $7}' ${file} | sed -e '1,2d' > pressure.dat
python ~/tools/others/calculate_average.py   --file_input  pressure.dat 

