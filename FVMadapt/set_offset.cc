//////////////////////////////////////////////////////////////////////////////////////////////////////
//                       set_offset.cc
// 
//
//////////////////////////////////////////////////////////////////////////////////////////////////////


#include <Loci.h>
#include <stdio.h>
#include "mpi.h"
#include "defines.h"
typedef Loci::vector3d<double> vect3d;
namespace Loci{
std::vector<int> all_collect_sizes(int size);
}





class get_node_offset : public pointwise_rule {
  const_param<int> num_original_nodes;
  const_store<int> num_inner_nodes;
  store<int> node_offset;
  
public:
  get_node_offset(){
    name_store("num_inner_nodes_copy", num_inner_nodes);
    name_store("num_original_nodes", num_original_nodes);
    name_store("node_offset", node_offset);
    input("num_inner_nodes_copy");
    input("num_original_nodes");
    output("node_offset");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;

  
    Loci::constraint edges, geom_cells, faces;
    Loci::storeRepP e2n = Loci::exec_current_fact_db->get_variable("edge2node");
    *edges = e2n->domain();
    faces = Loci::exec_current_fact_db->get_variable("faces");
    geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");


    // store<vect3d> pos ;
//     pos = Loci::exec_current_fact_db->get_variable("pos") ;
//     entitySet nodes = pos.domain() ;

    if(num_procs ==  1){
      
      int offset = *num_original_nodes + 1; //in cobalt, index of node start with 1
      
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
    fact_db::distribute_infoP d = Loci::exec_current_fact_db->get_distribute_info() ;
    Loci::constraint my_entities ; 
    my_entities = d->my_entities ;

    entitySet local_nodes, local_edges, local_faces, local_geom_cells;
    //  local_nodes = *my_entities & nodes ;
    //don't know if it's necessray
    local_edges = (*my_entities) & (*edges) ;
    local_faces = (*my_entities) & (*faces);
    local_geom_cells = (*my_entities)&(*geom_cells);
    
    
    // std::vector<int> local_original_nodes_sizes;
    // local_original_nodes_sizes = Loci::all_collect_sizes(local_nodes.size());
    //int num_original_nodes = 0;
    // for(unsigned int i = 0; i < local_original_nodes_sizes.size(); i++) num_original_nodes += local_original_nodes_sizes[i];
    
    
    //compute num_local_inner_nodes/num_local_fine_cells on each process
    int num_local_inner_nodes = 0;
    
    FORALL(local_edges, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    //write out cell_nodes first, then write face_nodes
    FORALL(local_geom_cells, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    
    FORALL(local_faces, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    
   
    std::vector<int> local_nodes_sizes;
    local_nodes_sizes = Loci::all_collect_sizes(num_local_inner_nodes);

    //each process computes its node  offset
    int noffset = *num_original_nodes +1;
    for(int i = 0; i < my_id; i++){
      noffset += local_nodes_sizes[i];
    }
    //compute the store values
    FORALL(local_edges, ei){
      node_offset[ei] = noffset;
      noffset += num_inner_nodes[ei];
    }ENDFORALL;
    //write out cell_nodes first, then write face_nodes
    FORALL(local_geom_cells, ei){
      node_offset[ei] = noffset;
      noffset += num_inner_nodes[ei];
    }ENDFORALL;
    
    FORALL(local_faces, ei){
      node_offset[ei] = noffset;
      noffset += num_inner_nodes[ei];
    }ENDFORALL;
    
  }
} ;
register_rule<get_node_offset> register_get_node_offset;


class get_cell_offset : public pointwise_rule {
  const_store<int> num_fine_cells;
  store<int> cell_offset;
    
public:
  get_cell_offset(){
    name_store("num_fine_cells_copy", num_fine_cells);
    name_store("cell_offset", cell_offset);
    input("num_fine_cells_copy");
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
    name_store("num_inner_nodes_copy", num_inner_nodes);
    input("num_inner_nodes_copy");
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







class init_ncells : public unit_rule{
  param<int> ncells;
public:
  init_ncells(){
    name_store("ncells", ncells);
    output("ncells");
    constraint("UNIVERSE");
    
  }
  //parameter, no loop
  virtual void compute(const sequence &seq){
    *ncells = 0;
  }
}; 
register_rule<init_ncells> register_init_ncells;

class apply_ncells : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> ncells;
  const_store<int> num_fine_cells;
public:
  apply_ncells(){
    name_store("ncells", ncells);
    name_store("num_fine_cells_copy", num_fine_cells);
    input("ncells");
    input("num_fine_cells_copy");
    output("ncells");
    constraint("geom_cells");
  }
  virtual void compute(const sequence &seq){
     do_loop(seq, this);
  }
  void calculate(Entity cc){
    join(*ncells, num_fine_cells[cc]);
  }
}; 
register_rule<apply_ncells> register_apply_ncells;


class init_nfaces : public unit_rule{
  param<int> nfaces;
  
public:
  init_nfaces(){
    name_store("nfaces", nfaces);
    output("nfaces");
    constraint("UNIVERSE");
  }
  //parameter, no loop
  virtual void compute(const sequence &seq){
    *nfaces = 0;
  }
}; 
register_rule<init_nfaces> register_init_nfaces;

class apply_nfaces : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> nfaces;
  const_store<Loci::FineFaces> fine_faces;
public:
  apply_nfaces(){
    name_store("nfaces", nfaces);
    name_store("fine_faces_copy", fine_faces);
    input("nfaces");
    input("fine_faces_copy");
    output("nfaces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    join(*nfaces, fine_faces[cc].size());
  }
}; 
register_rule<apply_nfaces> register_apply_nfaces;

