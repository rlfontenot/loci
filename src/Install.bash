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
cp FVMMod/fvm_m.so $INSTALL_PATH/lib

echo Installing Loci Tools
cp lpp/lpp $INSTALL_PATH/bin
cp FVMtools/cobalt2xdr $INSTALL_PATH/bin
cp FVMtools/extract $INSTALL_PATH/bin
cp FVMtools/pb $INSTALL_PATH/bin
cp FVMtools/plot3d2xdr $INSTALL_PATH/bin
cp FVMtools/solidMesh2xdr $INSTALL_PATH/bin
cp FVMtools/solidMesh2xdrMPI $INSTALL_PATH/bin
cp FVMtools/xdr2cobalt $INSTALL_PATH/bin
cp FVMtools/xdropt $INSTALL_PATH/bin

echo cp Loci.conf comp.conf sys.conf $INSTALL_PATH
cp Loci.conf comp.conf sys.conf $INSTALL_PATH

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include
cp include/*.lh $INSTALL_PATH/include
cp include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb ; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done

mkdir -p $INSTALL_PATH/docs
mkdir -p $INSTALL_PATH/docs/1D-Diffusion
mkdir -p $INSTALL_PATH/docs/heat
mkdir -p $INSTALL_PATH/docs/Datatypes

cp Tutorial/docs/tutorial.pdf $INSTALL_PATH/docs
cp Tutorial/1D-Diffusion/Makefile $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.loci $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.lh $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.cc $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.vars $INSTALL_PATH/docs/1D-Diffusion

cp Tutorial/heat/Makefile $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.loci $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.lh $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.xdr $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.vars $INSTALL_PATH/docs/heat

cp Tutorial/Datatypes/Makefile $INSTALL_PATH/docs/Datatypes
cp Tutorial/Datatypes/*.cc $INSTALL_PATH/docs/Datatypes

chmod -R a+rX $INSTALL_PATH
