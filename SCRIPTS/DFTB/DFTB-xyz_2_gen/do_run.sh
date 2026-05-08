gfortran  main.f90 -o a.out  -L/usr/local/lib -llapack -lblas -O

./a.out 286 13.5 13.5 13.5  ../geo_end.xyz  > OUTPUT

rm -rf a.out

