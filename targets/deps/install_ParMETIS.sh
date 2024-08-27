################################################################################
#!/bin/sh
# Filename:    install_ParMETIS.sh
# Author:      Mark A. Hunt (CFDRC)
# Date:        2024-07-04
# Description:
#
################################################################################

LOCI_SRC=$1
DEST_PREFIX=$2
GKLIB_BASE=$3
METIS_BASE=$4

ARCH=$(uname -s)
if [ "${ARCH}" == "Darwin" ]; then
  LIB_SUFFIX=dylib
else
  LIB_SUFFIX=so
fi

cd ${LOCI_SRC}/ext/ParMETIS

if [ -f ${DEST_PREFIX}/lib/libparmetis.${LIB_SUFFIX} ]; then
  echo "ParMETIS found in '${DEST_PREFIX}/lib/libparmetis.${LIB_SUFFIX}'"
  exit 0
fi

mkdir -p ${DEST_PREFIX}/lib
rm -rf build
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DGKLIB_PATH=${GKLIB_BASE} \
      -DMETIS_PATH=${METIS_BASE} \
      -DSHARED=on \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC -Wno-unused-command-line-argument" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ > config.out 2>&1 && \
make -j install > make.out
rm -rf *
cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DGKLIB_PATH=${GKLIB_BASE} \
      -DMETIS_PATH=${METIS_BASE} \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC -Wno-unused-command-line-argument" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ >> config.out 2>&1 && \
make -j install >> make.out
cd ../
rm -rf build
git checkout .
