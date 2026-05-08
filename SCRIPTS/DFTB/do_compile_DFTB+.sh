if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
  source /usr/local/tools/dotkit/init.sh
  use ic-17.0.174
  use mvapich2-intel-2.2
else
  module load intel impi
fi

#Compile code LSQ

make WITH_MPI=1 WITH_ARPACK=0
