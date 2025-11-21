module load lapack/3.12.0-gcc-11.2.1
module load gcc/11.2.1
module load cuda/12.2.2
module load essl/6.3.0.2

OLCF_CUDA_ROOT=$CUDA_HOME
OLCF_ESSL_ROOT=/usr/tcetmp/packages/essl/essl-6.3.0.2/
OLCF_NETLIB_SCALAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.12.0-gcc-11.2.1/
OLCF_NETLIB_LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.12.0-gcc-11.2.1/

echo ${OLCF_NETLIB_SCALAPACK_ROOT}
echo ${OLCF_NETLIB_LAPACK_ROOT}
echo ${OLCF_ESSL_ROOT}
echo ${OLCF_CUDA_ROOT}

mkdir INSTALL_DIR build
cd build

FC=mpif90 \
CC=mpicc \
CXX=mpicxx \
../configure --prefix=$(pwd)/../INSTALL_DIR \
  CFLAGS="-O2 -mcpu=power9" \
  CPP="cpp -E" \
  LDFLAGS="-L${OLCF_ESSL_ROOT}/lib64 \
           -L${OLCF_NETLIB_SCALAPACK_ROOT}/lib \
           -L${OLCF_CUDA_ROOT}/lib64 -lstdc++" \
  LIBS="-lcublas -lcublasLt -lcudart -lessl -lblas -lscalapack -llapack -lgfortran -lquadmath" \
  --enable-nvidia-gpu \
  --with-cuda-path=${OLCF_CUDA_ROOT} \
  --disable-sse-assembly --disable-sse --disable-avx --disable-avx2 --disable-avx512 \
  --enable-c-tests=no \
  --with-NVIDIA-GPU-compute-capability=sm_70

make clean
make -j
make install

cd ..
