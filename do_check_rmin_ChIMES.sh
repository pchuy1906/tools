grep -A10 "Minimum distances between atoms:"  fm_setup.out | tail -10 | awk '{printf "%4s %4s %4s %12.2f\n", $1, $2, $3, $4-0.15}'
