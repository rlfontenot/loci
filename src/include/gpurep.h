//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef GPUREP_H
#define GPUREP_H

// This header file contains the class definition of
// gpustoreRepI, store, and const_store.
// Their corresponding template implementation is in
// store_impl.h
// This separation is necessary to resolve some dependency
// problems in the class hierarchy.
// The same design applies to all other container classes.
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <DMap.h>
#include <store_rep.h>
#include <data_traits.h>
#include <sstream>
#include <hdf5_readwrite.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>

namespace Loci {
  class gpuRep : public storeRep {
  public:
    gpuRep() { } ;
    ~gpuRep() {} ;
    // Code to copy from cpu container to gpu container
    virtual void copyFrom(const storeRepP &p, entitySet set) = 0 ;
    // code to copy from gpu container to cpu container
    virtual void copyTo(storeRepP &p, entitySet set) const = 0 ;
  } ;

  typedef NPTR<gpuRep> gpuRepP ;
}
#endif
