cwd=`pwd`


cd tools/toolchain/

./install_cp2k_toolchain.sh \
  --with-gcc=no            \
  --with-intel=system      \
  --with-intel-classic=no  \
  --with-cmake=system      \
  --with-openmpi=system    \
  --with-mpich=no          \
  --with-mpich-device=no   \
  --with-intelmpi=no       \
  --with-libxc=no          \
  --with-libint=no         \
  --with-libgrpp=no        \
  --with-fftw=system       \
  --with-acml=no           \
  --with-mkl=system        \
  --with-openblas=system   \
  --with-scalapack=no      \
  --with-libxsmm=no        \
  --with-elpa=no           \
  --with-cusolvermp=no     \
  --with-ptscotch=no       \
  --with-superlu=no        \
  --with-pexsi=no          \
  --with-quip=no           \
  --with-plumed=no         \
  --with-sirius=no         \
  --with-gsl=no            \
  --with-libvdwxc=no       \
  --with-spglib=no         \
  --with-hdf5=no           \
  --with-spfft=no          \
  --with-spla=no           \
  --with-cosma=no          \
  --with-libvori=no        \
  --with-libtorch=no       

cd $cwd

  cp tools/toolchain/install/arch/*  arch/
  source $cwd/tools/toolchain/install/setup
  make -j ARCH=local VERSION="psmp"

