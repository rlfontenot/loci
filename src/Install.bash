#!/bin/bash

INSTALL_PATH=$INSTALL_DIR/$LOCI_INSTALL_DIR

echo INSTALL_PATH = $INSTALL_PATH

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib
mkdir -p $INSTALL_PATH/bin

echo Installing Library Files
cp Tools/libTools.so $INSTALL_PATH/lib
cp System/libLoci.so $INSTALL_PATH/lib

echo Installing Loci Tools
cp lpp/lpp $INSTALL_PATH/bin

echo cp Loci.conf comp.conf sys.conf $INSTALL_PATH
cp Loci.conf comp.conf sys.conf $INSTALL_PATH

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include
cp include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb ; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done


#echo Installing gcc 2.95 include fixes
#mkdir -p $INSTALL_PATH/include/g++-fixes
#cp -f include/g++-fixes/* $INSTALL_PATH/include/g++-fixes

chmod -R a+rX $INSTALL_PATH
