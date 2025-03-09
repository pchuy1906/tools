# SVD
python /g/g92/pham20/codes/LSQ/lsq2.py  --algorithm svd   --A A.txt --b b.txt --header params.header --map ff_groups.map --eps 1.0e-5  > params.txt

# LASSO
python /g/g92/pham20/codes/LSQ/lsq2.py  --algorithm lasso --A A.txt --b b.txt --header params.header --map ff_groups.map --eps 1.0e-5 --alpha 1.0e-5 > params.txt

