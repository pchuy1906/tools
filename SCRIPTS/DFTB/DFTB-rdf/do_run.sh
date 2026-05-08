gfortran  main.f90 -o a.out  -L/usr/local/lib -llapack -lblas -O

for atom1 in C H O N ; do
  for atom2 in C H O N ; do
    ./a.out ${atom1} ${atom2} 286  13.5 13.5 13.5  ../geo_end.xyz  > OUTPUT_${atom1}_${atom2}
    mv gr.dat gr_${atom1}_${atom2}.dat
  done
done

rm -rf a.out

