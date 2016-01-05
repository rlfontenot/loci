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
#include <gkeyspace.h>
#include <gkey_manager.h>
#include <iostream>
#include <limits>

using std::cout ;
using std::cerr ;
using std::endl ;
using std::vector ;
using std::string ;

#include <sstream>
using std::ostringstream ;
#include <map>
using std::map ;
#include <set>
using std::set ;

#include <parSampleSort.h>      // for balanceDistribution()
#include <distribute.h>
#include <Tools/except.h>

#include <execute.h>
#include <rule.h>
#include "sched_tools.h"

namespace Loci {

  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
  extern std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm);
  //static variable
  gKeySpace::gKeySpaceManager* gKeySpace::ksm=0;
  std::map<std::string, gKeySpaceP>  gKeySpace::gKeySpaceManager::space_map;
  
  void gKeySpace::set_send_recv_ptn(const std::vector<gEntitySet>& ptn){
      send_ptn = ptn;
      if( np>1 ) recv_ptn = transposePtn(ptn, comm);
      else recv_ptn = ptn;
  }
  

   gKeySpaceP gKeySpace::get_space(const std::string& spacename, const std::string& casename){
     return ksm->get_space(spacename, casename);
   }

  std::vector<gKeySpaceP> gKeySpace::get_all_spaces(){
    return ksm->get_all_spaces();
   }

  
  void gKeySpace::register_space(const std::string& spacename, const std::string& casename){
    ksm->add_space(spacename, casename, gKeySpaceP(this));
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

 
} // end of namespace Loci



