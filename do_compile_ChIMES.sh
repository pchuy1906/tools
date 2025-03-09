if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
  source /usr/local/tools/dotkit/init.sh
  use ic-17.0.174
  use mvapich2-intel-2.2
else
  module load intel impi
fi

#Compile code LSQ
rm -f *.o chimes_lsq
make -f Makefile-TS-LSQ chimes_lsq

#Compile code MD
rm -f *.o chimes_md
make -f Makefile-TS-MD chimes_md
