#!/usr/bin/bash

REVISION_NAME='$Name:  $'

INSTALL_DIR=${LOCI_INSTALL_DIR-/usr/local}

ARCH=`arch`
REV=`echo $REVISION_NAME| sed -e 's/.*: *//' -e 's/ *\$$//'`
INSTALL_PATH=$INSTALL_DIR/Loci-$ARCH-$REV/

echo INSTALL_PATH = $INSTALL_PATH

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib

echo Installing Library Files
cp Tools/libTools.a $INSTALL_PATH/lib
cp hdf5CC/libhdf5CC.a $INSTALL_PATH/lib
cp System/libLoci.a $INSTALL_PATH/lib

echo cp Loci.conf comp.conf sys.conf $INSTALL_PATH
cp Loci.conf comp.conf sys.conf $INSTALL_PATH

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include

for i in hdf5CC Tools Config; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done

cp -rf include/g++-fixes $INSTALL_PATH/include/

