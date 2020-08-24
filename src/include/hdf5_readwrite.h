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
#ifndef HDF5_READWRITE_H
#define HDF5_READWRITE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <entitySet.h>
#include <Map.h>
#include <data_traits.h>

namespace Loci {
  
  void HDF5_WriteDomain(hid_t group_id, const entitySet &en ) ;
  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g) ;
  void HDF5_ReadDomain( hid_t group_id, entitySet &eset) ;
  void HDF5_WriteVecSize(hid_t group_id, const int &size ) ;
  void HDF5_ReadVecSize(hid_t group_id, int *size ) ;

}
#endif
