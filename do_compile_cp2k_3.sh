module load intel-classic/2021.10.0

arch="Linux-intel-x86_64-minimal"

make -j 36 ARCH=${arch} VERSION=psmp realclean
make -j 36 ARCH=${arch} VERSION=psmp

