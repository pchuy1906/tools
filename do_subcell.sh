file='qe.in'
nat=`grep nat ${file} | awk '{print $3}'`

grep -A3 "CELL_PARA" ${file} | tail -3 > box
grep -A${nat} "ATOMIC_POSITIONS" ${file} | tail -${nat} > xyz

python ~/tools/others/subcell.py --option_xyz 0 --Nx 2 --Ny 2 --Nz 2

