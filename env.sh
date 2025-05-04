module purge
module load PrgEnv-gnu
module load gcc-native/12.3
module load cmake/3.24.3
module load cpe-cuda/23.02
module load cudatoolkit/11.7
module load python/3.9-anaconda-2021.11
module load craype-accel-nvidia80

export CRAY_ACCEL_TARGET=nvidia80
export CXX=CC
export CC=cc

PKG_PROP_OPT=ON




export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_CPU_TARGET=x86-64

export MPI_HOME=$MPICH_DIR
export MPI_CXX_COMPILER=$(which CC)
export MPI_CXX_COMPILER_FLAGS="$(CC --cray-print-opts=all) -I${CUDATOOLKIT_HOME}/../../math_libs/11.7/include"

#export CMAKE_PREFIX_PATH=/global/common/software/nersc/shasta2105/python/3.8-anaconda-2021.05:$CMAKE_PREFIX_PATH

export TOPDIR_HIP=/global/homes/f/fche123/m2322/fche123/code/build_script
export SRCROOT=${TOPDIR_HIP}/../src-20240121
export BUILDROOT=${TOPDIR_HIP}/build
export INSTALLROOT=${TOPDIR_HIP}/install
export LD_LIBRARY_PATH=${INSTALLROOT}/qmp/lib:${INSTALLROOT}/qdpxx/lib:${INSTALLROOT}/quda/lib:${INSTALLROOT}/chroma/lib:${INSTALLROOT}/llvm-14/lib:$CUDATOOLKIT_HOME/../../math_libs/11.7/lib64:$CUDATOOLKIT_HOME/../../math_libs/11.7/lib64:$CUDATOOLKIT_HOME/../../math_libs/lib64:${LD_LIBRARY_PATH}

module list
