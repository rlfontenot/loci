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

source ${LOCI_SRC}/targets/deps/installFunctions.sh

cd ${LOCI_SRC}/ext/ParMETIS

if [ -f ${DEST_PREFIX}/lib/libparmetis.so ]; then
  echo "ParMETIS found in '${DEST_PREFIX}/lib/libparmetis.so'"
  exit 0
fi

mkdir -p ${DEST_PREFIX}
rm -rf build
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DGKLIB_PATH=${GKLIB_BASE} \
      -DMETIS_PATH=${METIS_BASE} \
      -DSHARED=on \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ > config.out 2>&1 && \
make -j install > make.out
rm -rf *
cmake -DCMAKE_INSTALL_PREFIX=${DEST_PREFIX} \
      -DGKLIB_PATH=${GKLIB_BASE} \
      -DMETIS_PATH=${METIS_BASE} \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_C_FLAGS="-O3 -fPIC" \
      -DCMAKE_EXE_LINKER_FLAGS="-lmpi" \
      ../ >> config.out 2>&1 && \
make -j install >> make.out
cd ../
rm -rf build
git checkout .
