#!/bin/bash

module load craype-accel-amd-gfx942
module load PrgEnv-cray
module load rocm/6.2.1
#module load python
#module list
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

#export MPICH_DIR=/usr/tce/packages/cray-mpich/cray-mpich-8.1.32-cce-19.0.0-magic
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/cce/default/cce/x86_64/lib

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU

### FIXME ### Need a system wide install of libfabric from SHS 11 (or newer)
export LD_LIBRARY_PATH=/usr/workspace/wsb/accept/packages-2024/SHS11_lib:${LD_LIBRARY_PATH}

#export USE_SPINDLE='-o spindle'

#export CRAYPE_MPICH_GTL_DIR=${CRAY_MPICH_ROOTDIR}/gtl/lib
#export CRAYPE_MPICH_GTL_LIBS=-lmpi_gtl_hsa

export HIP_PATH=`hipconfig -p`
export LD_LIBRARY_PATH=${HIP_PATH}/lib:${LD_LIBRARY_PATH}

### Tell libfabric to only look for the ROCm runtime, not cuda, etc.
export FI_HMEM="rocr"

# Record/setup the compiler and FLAGS used to build LAMMPS and mpidevicemem
#export MPI_CXX_COMPILER=${HIP_PATH}/bin/hipcc
#export MY_MPI_CXX_FLAGS="-fdenormal-fp-math=ieee -fgpu-flush-denormals-to-zero -munsafe-fp-atomics -I${MPICH_DIR}/include"
#export MY_MPI_LINK_FLAGS="-L${MPICH_DIR}/lib -lmpi -L${CRAYPE_MPICH_GTL_DIR} ${CRAYPE_MPICH_GTL_LIBS} -lxpmem ${HUGEPAGES_LINK_FLAGS}"

############### start Huge Pages settings ###############

  # Share the read-only parts of executables built for libhugetlbfs
  # Don't do this because it creates persistent files in the hugetlbfs!
  # If the flux epilog properly cleans these up, then this would be safe to use
  #export HUGETLB_SHARE=1

  # Have shmget() calls use huge pages... does MPI do shmget()?
  # export HUGETLB_SHM=yes

  # Have malloc() calls use huge pages
  export HUGETLB_MORECORE=yes

  # Link flags to use libhugetlbfs, and align the exe segments to best do so
  export HUGEPAGES_LINK_FLAGS="-lhugetlbfs"
  # -Wl,--hugetlbfs-align isn't supported by the default linker
  # export HUGEPAGES_LINK_FLAGS="-lhugetlbfs -Wl,--hugetlbfs-align"

  # restrict libhugetlbfs to be enabled for these executables only:
  export HUGETLB_RESTRICT_EXE="defrag:lmp_elcap_kokkos_gpu"

  # Increase verbosity of libhugetlbfs
  #export HUGETLB_VERBOSE=3

  # Debug libhugetlbfs itself
  # export HUGETLB_DEBUG=1

  ############### end Huge Pages settings ###############

  export HSA_XNACK=1

# from Aurora, see https://docs.alcf.anl.gov/aurora/known-issues/#known-issues

#export FI_CXI_DEFAULT_CQ_SIZE=131072
#export FI_CXI_OFLOW_BUF_SIZE=8388608
#export FI_CXI_CQ_FILL_PERCENT=20
