//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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


#include "Config/conf.h"
#include <gkey_manager.h>

namespace Loci {
  gEntitySet
  distKeyManager::generate_key(gEntity size) {
    if(p>1){
      vector<gEntity> recv_buf(p) ;
      MPI_Allgather(&size, 1, MPI_GENTITY_TYPE, &recv_buf[0], 1, MPI_GENTITY_TYPE, comm) ;
      gEntity local_min = max_allocated+1 ;
      for(int i = 0; i < rank; ++i)
        local_min += recv_buf[i] ;
      gEntity total_size = 0 ;
      for(int i = 0; i < p; ++i)
        total_size += recv_buf[i] ;
      gEntity local_max = local_min + size -1;
      max_allocated += total_size;
      return gEntitySet(std::pair<gEntity, gEntity>(local_min, local_max));
    }else{
      gEntity local_min = max_allocated+1 ;
      gEntity local_max = local_min + size -1;
      max_allocated += size;
      return gEntitySet(std::pair<gEntity, gEntity>(local_min, local_max));
    }
  }



  //orthKeyManager implementation, not tested, for future use 
  orthKeyManager::orthKeyManager(gEntity global_start) {
    // first we obtain the limits on the gEntity type
    // we would like it to be slightly smaller than
    // the true max value
    gEntity global_max_key = std::numeric_limits<gEntity>::max() - 5 ;
    gEntity global_range = global_max_key - global_start ;
    // get the number of processes
    int comm_size = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size) ;
    gEntity local_range = global_range / comm_size ;

    rank = 0 ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    min_key = rank * local_range + global_start ;
    max_key = min_key + local_range - 1 ;

    int np = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD, &np) ;
    min_dist.resize(np) ;
    max_dist.resize(np) ;
    MPI_Allgather(&min_key, 1, MPI_INT,
                  &min_dist[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    MPI_Allgather(&max_key, 1, MPI_INT,
                  &max_dist[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    range_size = max_key - min_key + 1 ;

    // initially all are free keys
    freeset = interval(min_key, max_key) ;
  }

  
  gEntitySet
  orthKeyManager::generate_key(gEntity size) {
    if( (gEntity)freeset.size() < size) {
      std::stringstream err ;
      err << "orthKeyManager Error: process " << rank
          << " is running out of keys!" ;
      throw StringError(err.str()) ;
    }
    gEntitySet k ;
    for(gEntitySet::const_iterator i=freeset.begin();
        size!=0;++i,--size)
      k += *i ;

    freeset -= k ;

    return k ;
  }

  void
  orthKeyManager::recycle_key(gEntity k) {
    freeset += k ;
  }

  void
  orthKeyManager::recycle_key(const gEntitySet& keys) {
    freeset += keys ;
  }

}
