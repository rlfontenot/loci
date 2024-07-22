################################################################################
#!/bin/sh
# Filename:    install_METIS.sh
# Author:      Mark A. Hunt (CFDRC)
# Date:        2024-07-04
# Description:
#
################################################################################

LOCI_SRC=$1
DEST_PREFIX=$2
GKLIB_BASE=$3

source ${LOCI_SRC}/targets/deps/installFunctions.sh

cd ${LOCI_SRC}/ext/METIS


if [ -f ${DEST_PREFIX}/lib/libmetis.so ]; then
  echo "METIS found in '${DEST_PREFIX}/lib/libmetis.so'"
  exit 0
fi

mkdir -p ${DEST_PREFIX}
rm -rf build

/usr/bin/make config \
              prefix=${DEST_PREFIX} \
              gklib_path=${GKLIB_BASE} \
              shared=1 \
              cc=mpicc > config.out 2>&1 && \
make -j install > make.out && \
rm -rf build

make config \
     prefix=${DEST_PREFIX} \
     gklib_path=${GKLIB_BASE} \
     cc=mpicc >> config.out 2>&1 && \
/usr/bin/make -C build/ install >> make.out
rm -rf build
git checkout .

