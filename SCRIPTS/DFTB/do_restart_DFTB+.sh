prefix="geo_end"
natom=`head -1 ${prefix}.xyz | awk '{print $1}'`
neach=$(($natom+2))

tail -$neach geo_end.xyz > last_xyz.xyz
~/tools/others/xyzf_2_gen/a.out last_xyz.xyz 4 C H N O

filegen=`head -2 dftb_in.hsd | tail -1 | awk '{print $2}'`
echo $filegen
mv file.gen  $filegen

tail -$natom geo_end.xyz | awk '{print $6, $7, $8}' > last_velocity

#sed -i 's/MovedAtoms = 1:-1/MovedAtoms = 1:-1\nVelocities [AA\/ps] = \{\n<<< last_velocity\n\}/g'   dftb_in.hsd
#sed -i 's/ReadInitialCharges = No/ReadInitialCharges = Yes/g'  dftb_in.hsd

