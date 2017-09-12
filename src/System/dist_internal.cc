//#############################################################################
//#
//# Copyright 2015, Mississippi State University
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
#include <vector>
using std::vector ;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include "dist_internal.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>


namespace Loci {
  std::vector<int> all_collect_sizes(int size,MPI_Comm comm) {
    int prank=0,pnum = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&pnum) ;
    std::vector<int> vset( pnum) ;
    if(pnum > 1) {
      int *recv_count = new int[ pnum ] ;
      int *send_count = new int[ pnum ] ;
      for(int i = 0; i <  pnum; ++i) 
	send_count[i] = size ;
      
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,comm) ;
      for(int i = 0; i <  pnum ; ++i)
	vset[i] = recv_count[i] ;
      
      delete [] send_count ;
      delete [] recv_count ;
    }
    else
      vset[0] = size ;
    return vset ;
  }

}
