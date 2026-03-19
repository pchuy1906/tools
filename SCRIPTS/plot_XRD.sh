module load python/3.12.2

tail -250 XRD.dat > AAA.dat

lambda=`grep xrd in.lammps | head -1 | awk '{print $5}'`

python ~/tools/others/tool_XRD.py --_lambda $lambda

xmgrace  BBB.dat 
