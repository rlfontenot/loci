#!/bin/bash

REVISION_NAME='$Name: rel-1-0-beta-10 $'

INSTALL_DIR=${LOCI_INSTALL_DIR-/usr/local}

COMP_NAME=`echo $CPP | sed -e 's/ .*//' -e 's/.*\///'`
SYSTEM=`uname -s`
MACHINE=`uname -p`
if [ $SYSTEM == "Linux" ]; then MACHINE=`uname -m` ; fi

REV=`echo $REVISION_NAME| sed -e 's/.*: *//' -e 's/ *\$$//'`
echo $REV
echo $REVISION_NAME
# If no revision name, set the default to be month-day-year
if [ $REV == "" ]; then REV=`date +%m.%d.%y` ; fi

INSTALL_PATH=$INSTALL_DIR/Loci-$SYSTEM-$MACHINE-$COMP_NAME-$REV/

echo INSTALL_PATH = $INSTALL_PATH

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib

echo Installing Library Files
cp Tools/libTools.a $INSTALL_PATH/lib
cp System/libLoci.a $INSTALL_PATH/lib

echo cp Loci.conf comp.conf sys.conf $INSTALL_PATH
cp Loci.conf comp.conf sys.conf $INSTALL_PATH

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include

for i in  Tools Config ; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done


echo Installing gcc 2.95 include fixes
mkdir -p $INSTALL_PATH/include/g++-fixes
cp -f include/g++-fixes/* $INSTALL_PATH/include/g++-fixes

chmod -R a+rX $INSTALL_PATH
