fileXYZ="geo_end.xyz"
fileGEN="geo_end.gen"

natom=`head -1 $fileXYZ | awk '{print $1}'`
ntail=$(($natom+2))

tail -$ntail $fileXYZ > last.xyz
cell=`tail -3 $fileGEN | xargs`
sed -i '2s/./'"$cell"'  &/'  last.xyz

python ~/tools/others/xyz_2_gen.py --file_input   last.xyz 

tail -$natom last.xyz | awk '{print $6, $7, $8}' > last_vel.dat

fname="4_RESTART_0"
rm -rf $fname
mkdir $fname
mv last.xyz last_vel.dat $fname
mv file.gen $fname/last.gen
cp charges.bin $fname
