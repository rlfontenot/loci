################################################################################
#!/bin/sh
# Filename:    install_GKlib.sh
# Author:      Mark A. Hunt (CFDRC)
# Date:        2024-07-04
# Description:
#
################################################################################

LOCI_SRC=$1
DEST_PREFIX=$2

ARCH=$(uname -s)
if [ "${ARCH}" == "Darwin" ]; then
  LIB_SUFFIX=dylib
else
  LIB_SUFFIX=so
fi

cd ${LOCI_SRC}/ext/GKlib

if [ -f ${DEST_PREFIX}/lib/libGKlib.${LIB_SUFFIX} ]; then
  echo "GKlib found in '${DEST_PREFIX}/lib/libGKlib.${LIB_SUFFIX}'"
  exit 0
fi

mkdir -p ${DEST_PREFIX}
rm -rf build
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DBUILD_SHARED_LIBS=on \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ > config.out 2>&1 && \
make -j install > make.out && \
rm -rf * && \
cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ >> config.out 2>&1 && \
make -j install >> make.out
cd ../
rm -rf build
git checkout .
