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
//                       set_offset.cc
// 
//
//////////////////////////////////////////////////////////////////////////////////////////////////////


#include <Loci.h>
#include <stdio.h>
#include "mpi.h"
#include "defines.h"
using Loci::storeRepP;

typedef Loci::vector3d<double> vect3d;
namespace Loci{
  std::vector<int> all_collect_sizes(int size);
// Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);

   void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                        MPI_Comm comm);
  
}





class get_node_offset : public pointwise_rule {
  const_param<int> num_original_nodes;
  const_store<int> num_inner_nodes;
  store<int> node_offset;
  
public:
  get_node_offset(){
    name_store("num_inner_nodes", num_inner_nodes);
    name_store("num_original_nodes", num_original_nodes);
    name_store("node_offset", node_offset);
    input("num_inner_nodes");
    input("num_original_nodes");
    output("node_offset");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {

    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
   
  
    Loci::constraint edges, geom_cells, faces;
    Loci::storeRepP e2n = Loci::exec_current_fact_db->get_variable("edge2node");
    *edges = e2n->domain();
    faces = Loci::exec_current_fact_db->get_variable("faces");
    geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");


 

    if(num_procs ==  1 || dist == 0){
      
      // int offset = *num_original_nodes + 1; //in cobalt, index of node start with 1
      int offset = *num_original_nodes ; //in vog, index of node start with 0
      FORALL(*edges, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      //write out cell_nodes first, then write face_nodes
      FORALL(*geom_cells, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      
      FORALL(*faces, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      return;
    }
    
    Loci::constraint my_entities ; 
    my_entities = dist->my_entities ;

    entitySet local_nodes, local_edges, local_faces, local_geom_cells;
    //  local_nodes = *my_entities & nodes ;
    //don't know if it's necessray
    local_edges = (*my_entities) & (*edges) ;
    local_faces = (*my_entities) & (*faces);
    local_geom_cells = (*my_entities)&(*geom_cells);

     //each process computes its node  offset

    int noffset = *num_original_nodes ;

    int offset = 0;
    int num_local_nodes = 0;
    std::vector<int> local_nodes_sizes;
    storeRepP localVar=node_offset.Rep();
    //put into block to make the memory for store<int> freeed after it's used 
    {
      store<int> edge_num_inner_nodes;
      edge_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
      FORALL(edge_num_inner_nodes.domain(), ee){
        num_local_nodes += edge_num_inner_nodes[ee];
      }ENDFORALL;
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      
    
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> edge_file_offset;
      edge_file_offset.allocate(edge_num_inner_nodes.domain());
      FORALL(edge_num_inner_nodes.domain(), ei){
        edge_file_offset[ei] = noffset;
        noffset += edge_num_inner_nodes[ei];
      }ENDFORALL;
      File2LocalOrder(localVar, local_edges,
                      edge_file_offset.Rep(), offset,
                      dist,
                      MPI_COMM_WORLD);
      
    }
    //finish with edge nodes
    int total_num_edge_nodes = 0;
    for(int i = 0; i < num_procs; i++){
      total_num_edge_nodes += local_nodes_sizes[i];
    } 

    noffset = *num_original_nodes + total_num_edge_nodes;

    
   
    
    //repeat for geom_cells

 
    
    
    offset= 0;
    // Create container vardist that is ordered across processors in the
    // file numbering, the domain of this container shifted by offset
    // is the actual file numbering.
    {
      store<int> cell_num_inner_nodes;
      cell_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
   
      num_local_nodes = 0;
      FORALL(cell_num_inner_nodes.domain(), ei){
        num_local_nodes += cell_num_inner_nodes[ei];
      }ENDFORALL;
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> cell_file_offset;
      cell_file_offset.allocate(cell_num_inner_nodes.domain());
      FORALL(cell_num_inner_nodes.domain(), ei){
        cell_file_offset[ei] = noffset;
        noffset += cell_num_inner_nodes[ei];
      }ENDFORALL;
      //File2Local use the offset value set by Local2File    
      File2LocalOrder(localVar, local_geom_cells,
                      cell_file_offset.Rep(), offset,
                      dist,
                      MPI_COMM_WORLD);
      //finish with cell nodes
    }
    //update noffset
    int total_num_cell_nodes = 0;
    for(int i = 0; i < num_procs; i++){
      total_num_cell_nodes += local_nodes_sizes[i];
    } 
    noffset = *num_original_nodes + total_num_edge_nodes + total_num_cell_nodes;

    
    
     
    //repeat for faces
    offset = 0; 
    // Create container vardist that is ordered across processors in the
    // file numbering, the domain of this container shifted by offset
    // is the actual file numbering.
    {
      store<int> face_num_inner_nodes;
      face_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_faces, offset,dist,MPI_COMM_WORLD) ;
          
      num_local_nodes = 0;
      FORALL(face_num_inner_nodes.domain(), ei){
        num_local_nodes += face_num_inner_nodes[ei];
      }ENDFORALL;
      
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> face_file_offset;
      face_file_offset.allocate(face_num_inner_nodes.domain());
      FORALL(face_num_inner_nodes.domain(), ei){
        face_file_offset[ei] = noffset;
        noffset += face_num_inner_nodes[ei];
      }ENDFORALL;
      
      File2LocalOrder(localVar, local_faces,
                    face_file_offset.Rep(), offset,
                      dist,
                      MPI_COMM_WORLD);
      
    } 
  }
} ;
register_rule<get_node_offset> register_get_node_offset;


class get_cell_offset : public pointwise_rule {
  const_store<int> num_fine_cells;
  store<int> cell_offset;
    
public:
  get_cell_offset(){
    name_store("num_fine_cells", num_fine_cells);
    name_store("cell_offset", cell_offset);
    input("num_fine_cells");
    output("cell_offset");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    entitySet local_geom_cells = entitySet(seq);

   
    
    if(num_procs ==  1){
      int offset = 0;
      FORALL(local_geom_cells, ei){
        cell_offset[ei] = offset;
        offset += num_fine_cells[ei];
      }ENDFORALL;
      
      return;
    }
    
  
    
    //compute num_local_fine_cells on each process
    
    int num_local_fine_cells = 0;
    FORALL(local_geom_cells, ei){
      num_local_fine_cells += num_fine_cells[ei];
    }ENDFORALL;
    
   
    std::vector<int> local_cells_sizes;
    local_cells_sizes = Loci::all_collect_sizes(num_local_fine_cells);
    
    //each process computes its cell  offset
    int coffset = 0;
    for(int i = 0; i < my_id; i++){
      coffset += local_cells_sizes[i];
    }
  
   
    
    //compute the store values
    FORALL(local_geom_cells, ei){
      cell_offset[ei] = coffset;
      coffset += num_fine_cells[ei];
    }ENDFORALL;
    
   
  }
    
} ;
register_rule<get_cell_offset> register_get_cell_offset;




class init_npnts : public unit_rule{
  param<int> npnts;
public:
  init_npnts(){
    name_store("npnts", npnts);
    output("npnts");
    constraint("UNIVERSE");
  
  }
  //parameter, no loop, 
  virtual void compute(const sequence &seq){


    *npnts = 0;
  }
}; 
register_rule<init_npnts> register_init_npnts;

class apply_npnts : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> npnts;
  const_store<int> num_inner_nodes;
public:
  apply_npnts(){
    name_store("npnts", npnts);
    name_store("num_inner_nodes", num_inner_nodes);
    input("num_inner_nodes");
    input("npnts");
    output("npnts");
   
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    *npnts += num_inner_nodes[cc];
  }
}; 
register_rule<apply_npnts> register_apply_npnts;







// class init_ncells : public unit_rule{
//   param<int> ncells;
// public:
//   init_ncells(){
//     name_store("ncells", ncells);
//     output("ncells");
//     constraint("UNIVERSE");
    
//   }
//   //parameter, no loop
//   virtual void compute(const sequence &seq){
//     *ncells = 0;
//   }
// }; 
// register_rule<init_ncells> register_init_ncells;

// class apply_ncells : public apply_rule<param<int>, Loci::Summation<int> >{
//   param<int> ncells;
//   const_store<int> num_fine_cells;
// public:
//   apply_ncells(){
//     name_store("ncells", ncells);
//     name_store("num_fine_cells", num_fine_cells);
//     input("ncells");
//     input("num_fine_cells");
//     output("ncells");
//     constraint("geom_cells");
//   }
//   virtual void compute(const sequence &seq){
//      do_loop(seq, this);
//   }
//   void calculate(Entity cc){
//     join(*ncells, num_fine_cells[cc]);
//   }
// }; 
// register_rule<apply_ncells> register_apply_ncells;


// class init_nfaces : public unit_rule{
//   param<int> nfaces;
  
// public:
//   init_nfaces(){
//     name_store("nfaces", nfaces);
//     output("nfaces");
//     constraint("UNIVERSE");
//   }
//   //parameter, no loop
//   virtual void compute(const sequence &seq){
//     *nfaces = 0;
//   }
// }; 
// register_rule<init_nfaces> register_init_nfaces;

// class apply_nfaces : public apply_rule<param<int>, Loci::Summation<int> >{
//   param<int> nfaces;
//   const_store<Loci::FineFaces> fine_faces;
// public:
//   apply_nfaces(){
//     name_store("nfaces", nfaces);
//     name_store("fine_faces", fine_faces);
//     input("nfaces");
//     input("fine_faces");
//     output("nfaces");
//   }
//   virtual void compute(const sequence &seq){
//     do_loop(seq, this);
//   }
//   void calculate(Entity cc){
//     join(*nfaces, fine_faces[cc].size());
//   }
// }; 
// register_rule<apply_nfaces> register_apply_nfaces;

