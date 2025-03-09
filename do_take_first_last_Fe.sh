~/tools/LAMMPS/dump_2_xyz/a.out  dump_xyz_vxyz 1 Fe
ntail=`head -1 file.xyz | awk '{print $1+2}'`

head -$ntail file.xyz > first.xyz
tail -$ntail file.xyz >  last.xyz

python  ~/tools/others/xyz_2_gen.py  --file_input last.xyz  --num_cell_input  3

