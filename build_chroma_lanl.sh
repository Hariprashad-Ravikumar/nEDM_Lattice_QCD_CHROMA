source ./env.sh
pushd ${BUILDROOT}

if [ -d ./build_chroma_lanl ];
then
  rm -rf ./build_chroma_lanl
fi

mkdir  ./build_chroma_lanl
cd ./build_chroma_lanl
cmake ${SRCROOT}/chroma_lanl -DCMAKE_CXX_COMPILER=CC \
		-DCMAKE_C_COMPILER=cc -DCMAKE_C_STANDARD=99 -DCMAKE_C_EXTENSIONS=OFF  \
	 	-DBUILD_SHARED_LIBS=OFF \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DQDPXX_DIR=${INSTALLROOT}/qdpxx/lib/cmake/QDPXX \
		-DQMP_DIR=${INSTALLROOT}/qmp/lib/cmake/QMP \
	        -DLLVM_DIR=${INSTALLROOT}/llvm-14/lib/cmake/llvm \
	  -DChroma_ENABLE_QUDA=ON \
		-DQUDA_DIR=${INSTALLROOT}/quda/lib/cmake/QUDA \
		-DChroma_ENABLE_OPENMP=ON \
		-DChroma_ENABLE_JIT_CLOVER=ON \
		-DCMAKE_INSTALL_PREFIX=${INSTALLROOT}/chroma_lanl \
    -DLIBXML2_LIBRARY=${INSTALLROOT}/libxml2/lib/libxml2.a \
    -DLIBXML2_INCLUDE_DIR=${INSTALLROOT}/libxml2/include/libxml2 \
    -DCMAKE_CXX_LINK_FLAGS="-g $(xml2-config --libs) -lcufftw -lcufft -DCUFFTW -lm -L$CUDATOOLKIT_HOME/../../math_libs/11.7/lib64/" \
    -DCMAKE_CXX_FLAGS="-g $(xml2-config --cflags)"
    -DDISABLE_FORTRAN=ON
cmake --build . -j 16
cmake --install .

popd
#-DQUDA_DIR=/global/cfs/cdirs/m2322/junsik/qcd_progs/package-20220930_gpu/cuda-jit/install/quda/lib/cmake/QUDA \
