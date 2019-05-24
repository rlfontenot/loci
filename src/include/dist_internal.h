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
#ifndef DIST_INTERNAL_H
#define DIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <hdf5_readwrite.h>
#include <vector>
namespace Loci {

  void read_multi_vector_int(hid_t group_id, const char* name, int dim,  std::vector<int>& vint, MPI_Comm comm) ;
  void read_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size, MPI_Comm comm) ;
  void write_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, MPI_Comm comm) ;

  std::vector<int> all_collect_sizes(int size,MPI_Comm comm) ;
  inline std::vector<int> all_collect_sizes(int size) {
    return all_collect_sizes(size,MPI_COMM_WORLD) ;
  }

}

#endif
