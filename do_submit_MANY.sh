CURR_DIR=`pwd`

for k in $(seq 81 100 ); do
  cp -rf model  run-${k}
  cd run-$k
    k2=$(($k-1))
    ./do_submit.sh  ${CURR_DIR}/run-${k2}
  cd ..
done
