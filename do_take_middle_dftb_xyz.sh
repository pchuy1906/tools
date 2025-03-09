file="geo_end.xyz"
nall=`wc $file | awk '{print $1}'`
natom=`head -1 $file | awk '{print $1}'`
neach=$(($natom+2))
nstruc=$(($nall/$neach))
echo "# configurations is " ${nstruc}
ntake=$(($nstruc/2))

nbegin=$(($neach*($ntake-1)+1))
nend=$(($neach*($ntake-1)+$neach))

awk 'NR >= '"$nbegin"' && NR <= '"$nend"'' $file > middle.xyz
cell=`tail -3 geo_end.gen | xargs`
sed -i '2s/.*/'"$cell"'/'  middle.xyz
