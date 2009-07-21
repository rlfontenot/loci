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
//////////////////////////////////////////////////////////////////////////////////////////////////////
//                       set_bnodes_index.cc
// 
//
//////////////////////////////////////////////////////////////////////////////////////////////////////


#include <Loci.h>
#include <stdio.h>
#include "mpi.h"

typedef Loci::vector3d<double> vect3d;
namespace Loci{
std::vector<int> all_collect_sizes(int size);
}





class get_node_offset : public pointwise_rule {
  store<int> node_offset;
  
public:
  get_node_offset(){
    name_store("boundary_node_index", node_offset);
    output("boundary_node_index");
    constraint("boundary_nodes");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    
    if(num_procs ==  1){
      
      
      int offset = 0 ; //in vog, index of node start with 0
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        node_offset[*ei] = offset;
        offset++;
      }
      return;
    }
       
    int num_local_boundary_nodes = seq.size();
    std::vector<int> local_nodes_sizes;
    local_nodes_sizes = Loci::all_collect_sizes(num_local_boundary_nodes);

    
    int noffset = 0 ;
    
    for(int i = 0; i < my_id; i++){
      noffset += local_nodes_sizes[i];
    }
    
    for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
      node_offset[*ei] = noffset;
      noffset++;
    }
    
    
  }
} ;
register_rule<get_node_offset> register_get_node_offset;


