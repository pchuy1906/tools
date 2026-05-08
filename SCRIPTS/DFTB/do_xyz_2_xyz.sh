fileXYZ="geo_end.xyz"
fileGEN="geo_end.gen"

natom=`head -1 $fileXYZ | awk '{print $1}'`
ntail=$(($natom+2))

fileout="traj.xyz"
cp $fileXYZ $fileout
cell=`tail -3 $fileGEN | xargs`
sed -i '/MD/c '"$cell"' ' $fileout
