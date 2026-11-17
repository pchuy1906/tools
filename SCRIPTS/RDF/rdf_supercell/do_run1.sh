gfortran  main.f90 -o a.out  -L/usr/local/lib -llapack -lblas -O

for atom1 in N H; do
  for atom2 in N H; do
    ./a.out ${atom1} ${atom2}  ../traj.xyz  1 > OUTPUT_${atom1}_${atom2}
    mv gr.dat gr_${atom1}_${atom2}.dat
  done
done

rm -rf a.out

